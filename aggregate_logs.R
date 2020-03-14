aggregate_logs = function(dir,prefix='log',suffix='csv',stoich=14.7) {
# dir = 'C:/Users/Michael/Dropbox/Documents/!CAR STUFF/MX-5 Tuning/logs/tune21 (Jun-Aug 2019)/'

wid = 900;  ht = 600  # dimensions of PNG outputs
home = getwd()
setwd(dir)
loglist = Sys.glob(paste(prefix,'*',suffix,sep=''))
nlogs = length(loglist)

# Initialize data matrices
rpm_labels = seq(500,7500,500)  # define RPM bin labels
rpm_breaks = seq(250,7750,500)  # define RPM bin bounds
R = length(rpm_labels)
load_labels = round(seq(0,1,.0625),3)  # define load bin labels
load_breaks = c(0,seq(.03125,1.03125,.0625))  # define load bin bounds
L = length(load_labels)
mat0 = matrix(data=0,nrow=R,ncol=L,dimnames=list(as.character(rpm_labels),as.character(load_labels)))
CLcounts=mat0;  CLlambda=mat0;  CLkr=mat0;  CLtrim=mat0
OLcounts=mat0;  OLlambda=mat0;  OLkr=mat0
overview = data.frame(num_rows=vector('integer',nlogs), minutes=vector('numeric',nlogs), trim_avg=vector('numeric',nlogs),
                      ect_min=vector('numeric',nlogs), ect_max=vector('numeric',nlogs), 
                      num_wot=vector('integer',nlogs), hi_det=vector('logical',nlogs), kr_max=vector('numeric',nlogs))
row.names(overview) = loglist

# myinterp = function(xin, yin=NULL, xout) {return(approx(x=xin, y=yin, xout=xout, method='linear')$y)}

switch = TRUE
for (ilog in 1:nlogs) {                        # cycle over logs
  print(noquote(loglist[ilog]))
  data = read.csv(loglist[ilog], header=TRUE)  # load data
  nr = length(data[,1])                        # number of rows
  overview$num_rows[ilog] = nr
  overview$minutes[ilog] = round((data$time[nr]-data$time[1])/60,1)
  data = data[,3:length(data[1,])]             # remove first two columns
  data = data[2:nr,]                           # remove first row
  overview$kr_max[ilog] = round(max(data$kr),2)
  overview$ect_min[ilog] = min(data$ect)
  overview$ect_max[ilog] = max(data$ect)
  data = data[data$ect>70,]                    # eliminate cold coolant temps
  nr = length(data[,1])                        # number of rows
  if (nr == 0) next
  
  # keep = c('etc','lambda','status','ltft','stft','rpm','load','kr','hidet')
  # data = data[keep]                            # no longer need ect
  
  # # interpolate data
  # t1 = 1:nr
  # t2 = seq(1, nr, 0.5)
  # data = as.data.frame(apply(data, 2, 'myinterp', xin=t1, xout=t2))
  
  kr_event = c(0, diff(data$kr,1))
  data$kr[kr_event <= 0] = 0       # only keep KR values that have increased from the prior sample
  
  if (max(data$hidet) == 1) {    # eliminate hidet data
    print(noquote(' >>> WARNING: High detonation mode enabled'))
    ii = which(data$hidet == 1)
    ii = ii[2:length(ii)]        # keep first hidet observation
    data = data[-ii,]
    overview$hi_det[ilog] = TRUE
  } else {overview$hi_det[ilog] = FALSE}
  data = data[-which(data$status == 1),]                  # eliminate records with fuel status = 1
  data = data[-which((data$status==4) & (data$etc<10)),]  # eliminate open loop decel records
  if (length(data[,1]) == 0) next
  
  data$trim = (data$ltft + data$stft)/100                 # calculate total trim, and convert units
  data$irpm = as.integer(cut(data$rpm,rpm_breaks))        # Bin the RPM data
  data$iload = as.integer(cut(data$load,load_breaks))     # Bin the Load data
  CLdata = data[data$status==2,]                          # Closed Loop data
  
  grand_mean = mean(CLdata$trim[which(CLdata$rpm > 1500)])  # fuel trim grand mean
  overview$trim_avg[ilog] = round(grand_mean,4)
  nt = length(CLdata$trim)
  if (nt < 13000) {
    CLdata$trim = CLdata$trim - grand_mean  # convert absolute trim to relative trim
  } else {
    mean_vec = filter(CLdata$trim, filter=rep(1/12001,12001), method='convolution', sides=2)
    mean_vec[1:6000] = mean_vec[6001]
    mean_vec[(nt-5999):nt] = mean_vec[nt-6000]
    CLdata$trim = CLdata$trim - mean_vec
  }
  
  OLdata = data[((data$status==4) & (data$load>0.25)),]  # Open Loop data
  ii = which(OLdata$etc == 100)                          # WOT data
  if (length(ii) > 2) {
    jj=c(2,diff(ii,1));  kk=which(jj>1);  ii=ii[-kk]  # eliminate 1st observation of each WOT pull
    jj=c(2,diff(ii,1));  kk=which(jj>1);  ii=ii[-kk]  # eliminate 2nd observation of each WOT pull
    overview$num_wot[ilog] = length(kk)
    wot = OLdata[ii,]
    wot = wot[c('rpm','lambda','ltft')]
    if (switch) {
      WOTdata = wot
      switch = FALSE
    } else {
      WOTdata = rbind(WOTdata,wot)
    }
  } else {overview$num_wot[ilog] = 0}
  
  keep = c('irpm','iload','kr','lambda')
  rm('data');  OLdata = OLdata[keep];  CLdata = CLdata[c(keep,'trim')]  # clear memory of non-essentials
  
  ###########################################################################################
  # Calculate counts and sums by cell
  for (i in 1:R) {
    for (j in 1:L) {
      # Closed Loop
      dsub = subset(CLdata, (CLdata$irpm==i) & (CLdata$iload==j))
      n = length(dsub[,1])  # data count
      CLcounts[i,j] = CLcounts[i,j] + n
      if (n > 0) {
        CLtrim[i,j] = CLtrim[i,j] + sum(dsub$trim)
        CLlambda[i,j] = CLlambda[i,j] + sum(dsub$lambda)
        nKR = length(dsub$kr[which(dsub$kr > 0)])  # KR count
        if (nKR > 0) {
          nKR = nKR + length(dsub$kr[which(dsub$kr > 1.8)])
          nKR = nKR + 2*length(dsub$kr[which(dsub$kr > 3.6)])
        }
        CLkr[i,j] = CLkr[i,j] + nKR
      }
      
      # Open Loop
      dsub = subset(OLdata, (OLdata$irpm==i) & (OLdata$iload==j))
      n = length(dsub[,1])  # data count
      OLcounts[i,j] = OLcounts[i,j] + n
      if (n > 0) {
        nKR = length(dsub$kr[which(dsub$kr > 0)])  # KR count
        if (nKR > 0) {
          nKR = nKR + length(dsub$kr[which(dsub$kr > 1.8)])
          nKR = nKR + 2*length(dsub$kr[which(dsub$kr > 3.6)])
        }
        OLkr[i,j] = OLkr[i,j] + nKR
      }
    }
  }
}

# calculate mean values by cell
CLlambda = CLlambda/CLcounts
CLkr_pct = CLkr/CLcounts*100
OLkr_pct = OLkr/OLcounts*100
CLtrim = CLtrim/CLcounts

#########################################################################################
# Generate outputs
dir.create('output')
write.csv(CLcounts,'output/CLcounts.csv')
write.csv(OLcounts,'output/OLcounts.csv')
write.csv(round(CLlambda,2),'output/CLlambda.csv')
write.csv(round(CLkr_pct,2),'output/CLkr_pct.csv')
write.csv(round(OLkr_pct,2),'output/OLkr_pct.csv')
write.csv(round(CLtrim,3),'output/CLtrim.csv')
write.csv(overview,'output/overview.csv')

WOTdata = WOTdata[which(WOTdata$lambda < 0.925),]
wot_labels = seq(3125,7375,250)  # define RPM bin labels
rpm_breaks = seq(3000,7500,250)  # define RPM bin bounds
WOTdata$irpm = as.integer(cut(WOTdata$rpm,rpm_breaks))  # Bin the RPM data
WOTlambda = matrix(data=0, nrow=3, ncol=length(wot_labels))
for (i in 1:length(wot_labels)) {
  dsub = WOTdata$lambda[which(WOTdata$irpm==i)]
  WOTlambda[,i] = quantile(dsub, probs=c(0.05,0.5,0.95))
}

png(filename='output/wot_afr.png',width=wid,height=ht)
par(mar=c(4,5,2.5,1))
plot(WOTdata$rpm,stoich*WOTdata$lambda,xlim=c(3000,7500),ylim=c(11.6,13.6),xlab='',ylab='',yaxt='n',type='n')
mtext(side=3,line=1,'WOT AFR vs. Engine Speed',cex=1.1)
mtext(side=1,line=2.5,'Engine Speed (RPM)')
mtext(side=2,line=3.5,paste('Air Fuel Ratio (stoich = ',stoich,')',sep=''))
axis(side=2, at=seq(11.6,13.6,0.4),las=1)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4], col='grey90');  box()
abline(h=seq(11.6,13.6,0.4), col='white')
abline(v=seq(3000,7500,500), col='white')
ii = which(WOTdata$ltft < 2)
points(WOTdata$rpm[ii],stoich*WOTdata$lambda[ii],cex=1.2,col=rgb(red=0,green=0,blue=0,alpha=0.30))
ii = which(WOTdata$ltft>=2 & WOTdata$ltft<4)
points(WOTdata$rpm[ii],stoich*WOTdata$lambda[ii],cex=1.2,col=rgb(red=0,green=1,blue=0,alpha=0.40))
ii = which(WOTdata$ltft >= 4)
points(WOTdata$rpm[ii],stoich*WOTdata$lambda[ii],cex=1.2,col=rgb(red=1,green=0,blue=0,alpha=0.40))
lines(wot_labels,stoich*WOTlambda[1,],col='blue2',lwd=2,lty=2)
lines(wot_labels,stoich*WOTlambda[2,],col='blue3',lwd=3)
lines(wot_labels,stoich*WOTlambda[3,],col='blue2',lwd=2,lty=2)
legend('topright',inset=0.02,legend=c('LTFT < 2', '2 < LTFT < 4','LTFT > 4'),col=c('black','green','red'),pch=1,cex=1.2)
dev.off()

# png(filename='output/wot_lambda.png',width=wid,height=ht)
# par(mar=c(4,5,2.5,1))
# plot(WOTdata$rpm,WOTdata$lambda,xlim=c(3000,7500),ylim=c(0.8,0.925),xlab='',ylab='',yaxt='n',type='n')
# mtext(side=3,line=1,'WOT Lambda vs. Engine Speed',cex=1.1)
# mtext(side=1,line=2.5,'Engine Speed (RPM)')
# mtext(side=2,line=3.5,'Lambda')
# axis(side=2, at=seq(0.8,0.925,0.025),las=1)
# rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4], col='grey90');  box()
# abline(h=seq(0.8,0.925,0.025), col='white')
# abline(v=seq(3000,7500,500), col='white')
# ii = which(WOTdata$ltft < 2)
# points(WOTdata$rpm[ii],WOTdata$lambda[ii],cex=1.2,col=rgb(red=0,green=0,blue=0,alpha=0.25))
# ii = which(WOTdata$ltft>=2 & WOTdata$ltft<4)
# points(WOTdata$rpm[ii],WOTdata$lambda[ii],cex=1.2,col=rgb(red=0,green=1,blue=0,alpha=0.30))
# ii = which(WOTdata$ltft >= 4)
# points(WOTdata$rpm[ii],WOTdata$lambda[ii],cex=1.2,col=rgb(red=1,green=0,blue=0,alpha=0.30))
# lines(wot_labels,WOTlambda[1,],col='blue2',lwd=2,lty=2)
# lines(wot_labels,WOTlambda[2,],col='blue3',lwd=3)
# lines(wot_labels,WOTlambda[3,],col='blue2',lwd=2,lty=2)
# legend('topright',inset=0.02,legend=c('LTFT < 2', '2 < LTFT < 4','LTFT > 4'),col=c('black','green','red'),pch=1,cex=1.2)
# dev.off()

setwd(home)
print(noquote('Processing complete'))
}
