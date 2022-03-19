### This function is in development for astrochron: An R Package for Astrochronology
### Copyright (C) 2020 Stephen R. Meyers
###
###########################################################################
### function timeOptL - (SRM: April 12-14, 2017; September 23-24, 2020)
###
### modified from FUNCTION-timeOptL_v9.R, to perform timeOpt using formal
### likelihood functions
###########################################################################


timeOptL <- function (dat,sedmin=0.5,sedmax=5,numsed=100,linLog=1,limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,detrend=T,output=0,title=NULL,genplot=T,check=T,verbose=T)
{

if(verbose) cat("\n----- TimeOpt: Assessment of Amplitude Modulation & Bundling-----\n")

cormethod=1

# prepare data array
   dat = data.frame(dat)      
   npts <- length(dat[,1]) 
   dx <- dat[2,1]-dat[1,1]

if(check)
{
# initial error checking 
   if(dx<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[1], na.last = NA, decreasing = F), ]
       dx <- dat[2,1]-dat[1,1]
       npts <- length(dat[,1])
     }
   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }
}

if (verbose) 
 {
   cat(" * Number of data points in stratigraphic series:",npts,"\n")
   cat(" * Stratigraphic series length (meters):",(npts-1)*dx,"\n")
   cat(" * Sampling interval (meters):",dx,"\n\n")
 }

# detrend
if (detrend) 
  {
    lm.1 <- lm(dat[,2] ~ dat[,1])
    dat[2] <- dat[2] - (lm.1$coeff[2]*dat[1] + lm.1$coeff[1])
    if(verbose) cat(" * Linear trend subtracted. m=",lm.1$coeff[2],"b=",lm.1$coeff[1],"\n")
  }

# standardize data series
   dat[2]=dat[2]-colMeans(dat[2])
   dat[2]=dat[2]/sapply(dat[2],sd)

# convert sedmin and sedmax from cm/ka to m/ka for processing
   sedmin=sedmin/100
   sedmax=sedmax/100

# if sedmin equals sedmax, ensure numsed = 1
  if(sedmin == sedmax || numsed == 1) 
   {
     if(numsed != 1)
      {
        if(verbose) cat("\n**** WARNING: sedmin = sedmax, so numsed will be reset to 1.\n\n")
        numsed=1
      }  
   }  
   
  if(sedmin != sedmax && numsed == 1)   
   {
     cat("\n**** ERROR: sedmin is not equal to sedmax. You must specify numsed >1\n")
     stop("**** TERMINATING NOW!")
   }
  
#######################################################################################
# set up default bandpass frequencies and targets
#  first for precession
if(fit == 1)
 {
   if(is.null(flow)) 
    {
      flow = 0.035
      if(verbose) cat(" * Using default flow =",flow,"\n")
    }  
   if(is.null(fhigh)) 
    {
      fhigh = 0.065
      if(verbose) cat(" * Using default fhigh =",fhigh,"\n")
    } 
   if(is.null(roll))
    {
      roll = 10^3
      if(verbose) cat(" * Using default roll =",roll,"\n")
     }   
        
   if(is.null(targetP))
    {
# the four dominant precession peaks, based on spectral analysis of 
#   Laskar et al. (2004), 0-10 Ma
      targetP <- double(4)
      targetP[1] = 23.62069
      targetP[2] = 22.31868
      targetP[3] = 19.06768 
      targetP[4] = 18.91979 
      if(verbose) cat(" * Using default precession target periods (ka)=",targetP,"\n")
    }
  }

if(fit == 2 && !is.null(targetP)) 
  {
    if(verbose) cat("\n**** WARNING: targetP is defined but will not be used in fitting!\n")
  }

# next for short eccentricity
if(fit == 2)
 {
   if(is.null(flow))
    {
      flow = 0.007
      if(verbose) cat(" * Using default flow =",flow,"\n")
    }  
   if(is.null(fhigh))
    {
      fhigh = 0.0115
      if(verbose) cat(" * Using default fhigh =",fhigh,"\n")
     }
   if(is.null(roll))
    {
      roll = 10^5
      if(verbose) cat(" * Using default roll =",roll,"\n")
     }   
 }  

if(is.null(targetE))
 {
# the five domintant eccentricity peaks based on spectral analysis of LA10d solution 
#   (Laskar et al., 2011), 0-20 Ma
     targetE <- double(5)
     targetE[1] = 405.6795
     targetE[2] = 130.719
     targetE[3] = 123.839
     targetE[4] = 98.86307
     targetE[5] = 94.87666
     if(verbose) cat(" * Using default eccentricity target periods (ka)=",targetE,"\n")     
  }   

# targetTot is for plotting, and fitting if precession modulations assessed
if(fit == 1 && fitModPwr) targetTot = c(targetE,targetP)
if(fit == 1 && !fitModPwr) targetTot = c(targetP)
if(fit == 2 && fitModPwr) targetTot = c(targetE)
if(fit == 2 && !fitModPwr) targetTot = c(targetE[-1])

if(flow < 0.5/max(targetTot) && verbose) cat("\n**** NOTE: flow is less than half of the smallest target frequency. Did you intend this?\n")
if(fhigh > 2/min(targetTot) && verbose) cat("\n**** NOTE: fhigh is 2 times greater than the largest target frequency. Did you intend this?\n")  

# check minimum and maximum sedimentation rates. sedmin is now in m/ka, dx is in meters.
   NyqFreq=sedmin/(2*dx)
# fhigh is the upper half-power point for the filter   
   if(fhigh>NyqFreq)
    {
     if(limit) 
      {
        sedmin = 2*dx*fhigh
        if(verbose) cat("\n**** WARNING: minimum sedimentation rate is too low for full signal recovery.\n")
        if(verbose) cat("              sedmin reset to",100*sedmin,"cm/ka\n\n")
      }
     if(verbose && !limit) 
      {
# note: when sedimentation rates exceed this value, the filter starts to behave more-and-more like a high-pass filter
        cat("\n**** WARNING: minimum sedimentation rate is too low for full signal recovery, but it will still be evaluated. USE WITH CAUTION!\n")
        cat("              The high frequency half-power point of filter is exceeded when sedrate=",100*2*dx*fhigh,"cm/ka\n")
        if(fit == 1)
         {
           cat("              The shortest precession period is undetectable when sedrate <",100*2*dx/min(targetP),"cm/ka\n")
           cat("              The longest precession period is undetectable when sedrate <",100*2*dx/max(targetP),"cm/ka\n\n")
         }
        if(fit == 2)
         {   
           cat("              The shortest eccentricity period is undetectable when sedrate <",100*2*dx/min(targetE),"cm/ka\n")
           cat("              The longest eccentricity period is undetectable when sedrate <",100*2*dx/max(targetE),"cm/ka\n\n")
         }  
      }
    }

# check maximum sedimentation rate. sedmax is in m/ka. dx is in meters.
RayFreq = sedmax/(npts*dx)
# freqLow identifies the frequency of the longest period in the target.
freqLow=1/max(targetE)
if(RayFreq>freqLow)
  {
    if(limit) 
      {
        sedmax = npts*dx*freqLow
        if(verbose) cat("\n**** WARNING: maximum sedimentation rate is too high for full signal recovery.\n")
        if(verbose) cat("              sedmax reset to",100*sedmax,"cm/ka\n\n")
      }
     if(!limit) 
      {
        if(verbose) cat("\n**** WARNING: maximum sedimentation rate is too high for full signal recovery, but it will still be evaluated. USE WITH CAUTION!\n")
        cat("              The longest eccentricity period is undetectable when sedrate >",100*npts*dx*freqLow,"cm/ka\n")
        cat("              The shortest eccentricity period is undetectable when sedrate >",100*npts*dx/min(targetE),"cm/ka\n\n")
      }  
  }

if(sedmin>sedmax)
  {
    cat("\n**** ERROR: sedmin > sedmax\n")
    stop("**** TERMINATING NOW!")
  }

#######################################################################################
# Definition of FUNCTIONS: genCycles, fitIt, calcLogLH
# function to generate cos (real) and sin (imaginary) terms for each target period, 
#   and convert to spatial cycles, given a particular sed rate in m/ka
genCycles <- function(sedrate1, targetIn, n) 
  {
    x <- matrix(0, n, 2*length(targetIn))
    for (i in 1:length(targetIn)) 
      {
        x[,2*i-1] <- cos( (2*pi)/(targetIn[i]) * (dx/sedrate1) * (1:n))
        x[,2*i] <- sin( (2*pi)/(targetIn[i]) * (dx/sedrate1) * (1:n))
      }
    return(x)
  }
  
# function to perform fitting
#  dx, npts passed into function transparently
fitIt <- function(sedrate1,timeSeries,targetIn) 
  {
    xm <- genCycles(sedrate1, targetIn, npts)
    lm.0 <- lm(timeSeries[,2] ~ xm)
# calculate rho and sigma based on residuals
    rho = cor(lm.0$residuals[1:(npts-1)],lm.0$residuals[2:npts])
    sigma = sd(lm.0$residuals)
    logLH = calcLogLH(lm.0$residuals,rho,sigma)
    return(cbind(sedrate1,logLH,rho,sigma))
   } 
  
# function to perform log-likelihood calculation, including
#  assessment of correlated residuals, assuming AR1 model
#  npts passed into function transparently
calcLogLH <- function(residuals,rho,sigma)
 {
# calculate Re^-1 (ReInv), as in EQ. A-7 of Malinverno & Briggs (2004)
# set up array, with zeros
   ReInv <- double(npts*npts)
   dim(ReInv) <-c(npts,npts)
# put 1+rho^2 on diagonal
   ReInv[row(ReInv)==col(ReInv)] = 1+rho^2
# except at (1,1) and (npts,npts), which have a value of 1
   ReInv[1,1] = 1
   ReInv[npts,npts] = 1
# put -rho on subdiagonal  
   ReInv[(row(ReInv)-1)==col(ReInv)] = -1*rho
# put -rho on superdiagonal  
   ReInv[(row(ReInv)+1)==col(ReInv)] = -1*rho
# now multiple matrix by 1/(1-rho^2)
   ReInv = ReInv/(1-rho^2)
# calculate log-likelihood 
   logLH2 = ( npts*log(2*pi) ) + ( 2*npts*log(sigma) ) + ( (npts-1)*log(1-rho^2) )
   logLH2 = logLH2 + (1/(sigma^2)) * t(residuals) %*% ReInv %*% residuals 
   logLH2 = -0.5*logLH2
   return(logLH2)
 }


#######################################################################################  
# set up sedimentation rate grid array, dimension appropriately
# 'ans' will contain sedrate, envelope likelihood, power likelihood, env. rho, env. sigma, power rho, power sigma 
ans<- rep(NA,numsed*7)
dim(ans) <- c(numsed,7)
sedrate <-double(numsed)

if(numsed == 1) sedrate=sedmin

if(numsed != 1)
 {
if(linLog==0)
  {
    sedinc = (sedmax-sedmin)/(numsed-1)
    for (ii in 1:numsed)
      {
         sedrate[ii] = sedmin + ((ii-1)*sedinc)
      }
   }
   
if(linLog==1)
  {
    sedinc = (log10(sedmax)-log10(sedmin))/(numsed-1)
    for (ii in 1:numsed)
      {
         sedrate[ii] = log10(sedmin) + ((ii-1)*sedinc)
         sedrate[ii] = 10^sedrate[ii]
      }
   }
 }
 
 
#######################################################################################
# begin sedimentation rate loop


if(verbose) 
  {  
    cat("\n * PLEASE WAIT: Performing Optimization\n")
    cat("\n0%       25%       50%       75%       100%\n")
# create a progress bar
    progress = utils::txtProgressBar(min = 0, max = numsed, style = 1, width=43)
  }

i=0
for (ii in 1:numsed) 
{
    if(verbose) utils::setTxtProgressBar(progress, ii)
    i=i+1
# CALIBRATE DEPTH SERIES (m) TO TIME (ka)
    ts = dat
# create new time vector
# it is the index vector for time
    it <- seq(1,npts,by=1)
    time = (dx/sedrate[ii]) * (it-1)
    ts[1] = time

# bandpass precession or short eccentricity band
    bp = taner(ts,padfac=2,flow=flow,fhigh=fhigh,roll=roll,demean=T,detrend=F,addmean=F,genplot=F,verbose=F)

# hilbert transform for instantaneous amplitude
    hil = hilbert(bp,padfac=2,demean=T,detrend=F,addmean=F,genplot=F,verbose=F)
# standardize hil to unit variance
    hil[2]=hil[2]-colMeans(hil[2])
    hil[2]=hil[2]/sd(hil[,2])    

# execute functions
# for precession modulations
    if(fit == 1) 
      {
         res = fitIt(sedrate[ii],hil,targetE)
         pwrOut = fitIt(sedrate[ii],ts,targetTot)
       }
# for short eccentricity modulations
    if(fit == 2) 
      {
         res = fitIt(sedrate[ii],hil,targetE[1])
         pwrOut = fitIt(sedrate[ii],ts,targetE)
       }
    
    ans[i,1] <- res[1]
    ans[i,2] <- res[2]
    ans[i,3] <- pwrOut[2]
    ans[i,4] <- res[3]
    ans[i,5] <- res[4]
    ans[i,6] <- pwrOut[3]
    ans[i,7] <- pwrOut[4]

# end sedimentation rate loop
}
   if(verbose) close(progress) 
   rPwr = ans[,2] + ans[,3]
   
#######################################################################################
# find sedimentation rate with maxima in power
ptPwr=which.max(ans[,3])
# check for multiple maxima
sumMax = sum( ans[,3] == max(ans[,3]) )
if(sumMax > 1 && verbose) cat("\n**** WARNING: Multiple spectral power maxima detected.\n") 
# find sedimentation rate with maxima in correlation
ptCor=which.max(ans[,2])
# check for multiple maxima
sumMax = sum( ans[,2] == max(ans[,2]) )
if(sumMax > 1 && verbose) cat("\n**** WARNING: Multiple envelope maxima detected.\n") 
# find sedimentation rate with maxima in r*p
ptRp=which.max(rPwr)
# check for multiple maxima
sumMax = sum( rPwr == max(rPwr) )
if(sumMax > 1 && verbose) cat("\n**** WARNING: Multiple (spectral power*envelope) maxima detected.\n") 

if(verbose)
 {  
  cat("\n * Maximum (spectral power likelihood)=", ans[ptPwr,3],"at sedimentation rate of", 100*ans[ptPwr,1],"cm/ka\n")
  cat(" * Maximum (envelope likelihood)=", ans[ptCor,2],"at sedimentation rate of", 100*ans[ptCor,1],"cm/ka\n")
  cat(" * Maximum (envelope likelihood) + (spectral power likelihood) =", rPwr[ptRp],"at sedimentation rate of", 100*ans[ptRp,1],"cm/ka\n")
 }
  

#######################################################################################
if(genplot == T || output == 2)
 { 
# plotting
# recalculate bandpasses and hilbert of optimal sedimentation rate for plotting
    ts = dat
# create new time vector
# index vector for time
    it <- seq(1,npts,by=1)    
    time = (dx/ans[ptRp,1]) * (it-1)
    ts[1] = time 

# filter record    
    bp = taner(ts,padfac=2,flow=flow,fhigh=fhigh,roll=roll,demean=T,detrend=F,addmean=F,genplot=F,verbose=F)
# get filter window
    bpWin= taner(ts,padfac=2,flow=flow,fhigh=fhigh,roll=roll,demean=T,detrend=F,addmean=F,genplot=F,verbose=F,output=2)
    hil = hilbert(bp,padfac=2,demean=T,detrend=F,addmean=F,genplot=F,verbose=F)
 
# perform fitting at optimal sedimentation rate for plotting
if(fit == 1) 
 {
   xm <- genCycles(ans[ptRp,1], targetE, npts)
   xm2 <- genCycles(ans[ptRp,1], targetTot, npts)
 }  
if(fit == 2) 
 {
   xm <- genCycles(ans[ptRp,1], targetE[1], npts)
   xm2 <- genCycles(ans[ptRp,1], targetE, npts)
 }
lm.0 <- lm(hil[,2] ~ xm)
lm.2 <- lm(ts[,2] ~ xm2)
 }
 
if(genplot)
 { 
  if(is.null(title)) title = c("TimeOpt Results")
  if(numsed != 1) 
   {
     dev.new(title=title,height=7,width=7)
     par(mfrow=c(3,2))

     plot(100*ans[,1],ans[,2],cex=.75,cex.lab=1.2,cex.main=1.3,col="red",xlab="",ylab="",main="Fit: envelope (red) and power (gray)")
# plot red numbers on left axis
     axis(2,col.axis="red")
     mtext("Log Likelihood",side=2,line=2,cex=0.9)
     mtext("Sedimentation rate (cm/ka)",side=1,line=2.3,cex=0.8)
     par(new=T)
     plot(100*ans[,1],ans[,3],col="#00000064",xlab="",ylab="",type="l",axes=F,lwd=2)
     axis(4, ylim=c(0,max(ans[,3])),lwd=1,col="black")
        
     plot(100*ans[,1],rPwr,type="l",lwd=2,cex.lab=1.2,cex.main=1.3,col="black",xlab="",ylab="",main="Optimal Fit")
#    points(100*ans[,1],rPwr,cex=.75,pch=21,bg="white")
     mtext("Log Likelihood",side=2,line=1.9,cex=0.9)
     mtext("Sedimentation rate (cm/ka)",side=1,line=2.3,cex=0.8)

# plot hil and ecc estimated by least squares fitting
     ylim1=c(min(hil[,2],lm.0$fitted),max(hil[,2],lm.0$fitted))
     plot(hil,cex=.5,cex.lab=1.2,cex.main=1.2,xlab="",ylab="",main="Envelope (red); Reconstructed Ecc. Model (black)", col="red",type="l",ylim=ylim1)
     lines(hil[,1],lm.0$fitted,lwd=1.5)     
     mtext("Std. Value",side=2,line=2,cex=0.9)
     mtext("Time (ka)",side=1,line=2.3,cex=0.8)

# plot bp and hilbert of bp
     ylim1=c(min(bp[,2],-1*hil[,2]),max(bp[,2],hil[,2]))
     plot(bp,col="blue",cex=.5,cex.lab=1.2,cex.main=1.3,xlab="",ylab="",main="Envelope (red); Filtered Data (blue)",ylim=ylim1)
     lines(bp,col="blue")
     lines(hil,col="red")
     lines(hil[,1],-1*hil[,2],col="red")
     abline(h=0,col="black",lty=3)
     mtext("Std. Value",side=2,line=2,cex=0.9)
     mtext("Time (ka)",side=1,line=2.3,cex=0.8)

# cross plot of hil and ecc
     plot(hil[,2],lm.0$fitted, cex.lab=1.2, cex.main=1.3, main="Envelope vs. Reconstructed Ecc. Model", xlab="",ylab="")
     mtext("Reconst. Ecc. Model",side=2,line=2,cex=0.9)
     mtext("Data Envelope",side=1,line=2.3,cex=0.8)

# periodogram of stratigraphic series
     fft = periodogram( data.frame(cbind(ts[,1],ts[,2])), output=1, verbose=F,genplot=F)
# remove f(0) for log plot of power
     fft = subset(fft,(fft[,1] > 0))
     plot(fft[,1],fft[,3],cex.lab=1.2,cex.main=1.3,xlim=c(0,0.1),type="l",xlab="",ylab="",main="Power Spectrum (black=linear; gray=log)")
     lines(bpWin[,1],bpWin[,2]*max(fft[,3]),col="blue")
     mtext("Power",side=2,line=2,cex=0.9)
     mtext("Frequency (cycles/ka)",side=1,line=2.3,cex=0.8)
# plot a second y-axis
     par(new=TRUE)
     plot(fft[,1],log(fft[,3]),xlim=c(0,0.1),type="l",yaxt="n",col="gray",xlab="",ylab="")
     abline(v=1/targetTot, col="red",lty=3)
#     abline(v=c(flow,fhigh), col="blue",lty=2)   
   }

  if(numsed == 1) 
   {
     dev.new(title=title,height=7,width=7)
     par(mfrow=c(3,2))

# plot bp and hilbert of bp
     ylim1=c(min(bp[,2],-1*hil[,2]),max(bp[,2],hil[,2]))
     plot(bp,col="blue",cex=.5,cex.lab=1.2,cex.main=1.3,xlab="",ylab="",main="Envelope (red); Filtered Data (blue)",ylim=ylim1)
     lines(bp,col="blue")
     lines(hil,col="red")
     lines(hil[,1],-1*hil[,2],col="red")
     abline(h=0,col="black",lty=3)
     mtext("Std. Value",side=2,line=2,cex=0.9)
     mtext("Time (ka)",side=1,line=2.3,cex=0.8)

# periodogram of stratigraphic series
     fft = periodogram( data.frame(cbind(ts[,1],ts[,2])), output=1, verbose=F,genplot=F)
# remove f(0) for log plot of power
     fft = subset(fft,(fft[,1] > 0))
     plot(fft[,1],fft[,3],cex.lab=1.2,cex.main=1.3,xlim=c(0,0.1),type="l",xlab="",ylab="",main="Power Spectrum (black=linear; gray=log)")
     lines(bpWin[,1],bpWin[,2]*max(fft[,3]),col="blue")
     mtext("Power",side=2,line=2,cex=0.9)
     mtext("Frequency (cycles/ka)",side=1,line=2.3,cex=0.8)
# plot a second y-axis
     par(new=TRUE)
     plot(fft[,1],log(fft[,3]),xlim=c(0,0.1),type="l",yaxt="n",col="gray",xlab="",ylab="")
     abline(v=1/targetTot, col="red",lty=3)
#     abline(v=c(flow,fhigh), col="blue",lty=2)   

# plot hil and ecc estimated by least squares fitting
     ylim1=c(min(hil[,2],lm.0$fitted),max(hil[,2],lm.0$fitted))
     plot(hil,cex=.5,cex.lab=1.2,cex.main=1.2,xlab="",ylab="",main="Envelope (red); Reconstructed Ecc. Model (black)", col="red",type="l",ylim=ylim1)
     lines(hil[,1],lm.0$fitted,lwd=1.5)     
     mtext("Std. Value",side=2,line=2,cex=0.9)
     mtext("Time (ka)",side=1,line=2.3,cex=0.8)

# plot calibrated times series and fit estimated by least squares fitting    
     ylim2=c(min(ts[,2],lm.2$fitted),max(ts[,2],lm.2$fitted))
     plot(ts,cex=.5,cex.lab=1.2,cex.main=1.2,xlab="",ylab="",main="Data (red); Full Regression Model (black)", col="red",type="l",ylim=ylim2)
     lines(ts[,1],lm.2$fitted,lwd=1.5)     
     mtext("Std. Value",side=2,line=2,cex=0.9)
     mtext("Time (ka)",side=1,line=2.3,cex=0.8)
    
# cross plot of hil and ecc
     plot(hil[,2],lm.0$fitted, cex.lab=1.2, cex.main=1.3, main="Envelope vs. Reconstructed Ecc. Model", xlab="",ylab="")
     mtext("Reconst. Ecc. Model",side=2,line=2,cex=0.9)
     mtext("Data Envelope",side=1,line=2.3,cex=0.8)

# cross plot of data and full regression model
     plot(ts[,2],lm.2$fitted, cex.lab=1.2, cex.main=1.3, main="Data vs. Full Regression Model", xlab="",ylab="")
     mtext("Regression Model",side=2,line=2,cex=0.9)
     mtext("Data",side=1,line=2.3,cex=0.8)    
   }      
# end genplot section
 }
     
# return sedimentation rate grid, envelope, spectral power, envelope+spectral power 
     if(output == 1) 
       {
         out <- data.frame (cbind(100*ans[,1],ans[,2],ans[,3],rPwr,ans[,4],ans[,5],ans[,6],ans[,7]) )
         colnames(out) <- c("sedrate","likelihood_env","likelihood_spec","likelihood_combined","env_rho","env_sigma","pwr_rho","pwr_sigma")
         return(out)
       }
       
# return optimal time series, hilbert and fitted periods
     if(output == 2) 
      {
        out <- data.frame (cbind(ts[,1],ts[,2],bp[,2],hil[,2],lm.0$fitted) )
        colnames(out) <- c("time","value","filtered_precession","precession_envelope","reconstructed_ecc_model")
        return(out)
      }  

### END function timeOptL
}
