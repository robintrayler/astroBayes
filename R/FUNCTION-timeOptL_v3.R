### This function is in development for astrochron: An R Package for Astrochronology
### Copyright (C) 2020 Stephen R. Meyers
###
###########################################################################
### function timeOptL - (SRM: April 12-14, 2017; September 23-24, 2020)
###
### modified from FUNCTION-timeOptL_v9.R, to perform timeOpt using formal
### likelihood functions
###########################################################################


time_opt_likelihood <- function (dat,
                                 sed_rate = 1,
                                 linLog=1,
                                 limit=T,
                                 fit=1,
                                 fitModPwr=T,
                                 flow=NULL,
                                 fhigh=NULL,
                                 roll=NULL,
                                 targetE=NULL,
                                 targetP=NULL,
                                 output=0,
                                 title=NULL,
                                 genplot=T,
                                 verbose=T)
{

  # prepare data array
  dat = data.frame(dat)
  npts <- length(dat[, 1])
  dx <- dat[2, 1] - dat[1, 1]

  # standardize data series
  dat[2]=dat[2] - colMeans(dat[2])
  dat[2]=dat[2] / sapply(dat[2],sd)

  # convert sed_rate cm/ka to m/ka for processing
  sed_rate = sed_rate / 100

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
    # the five dominant eccentricity peaks based on spectral analysis of LA10d solution
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

  #######################################################################################
  # Definition of FUNCTIONS: genCycles, fitIt, calcLogLH
  # function to generate cos (real) and sin (imaginary) terms for each target period,
  #   and convert to spatial cycles, given a particular sed rate in m/ka
  genCycles <- function(sed_rate,
                        target_frequencies,
                        n) {
    # set up storage matrix
    storage <- matrix(data = 0,
                      nrow = n,
                      ncol = 2 * length(target_frequencies))

    # loop through target frequencies
    for (i in 1:length(target_frequencies)) {
      storage[,2*i-1] <- cos( (2*pi) / (target_frequencies[i]) * (dx / sed_rate) * (1:n))
      storage[,2*i] <- sin( (2*pi) / (target_frequencies[i]) * (dx / sed_rate) * (1:n))
    }
    return(storage)
  }

  # function to perform fitting
  #  dx, npts passed into function transparently
  fitIt <- function(sed_rate,
                    time_series,
                    target_frequencies) {
    xm <- genCycles(sed_rate, target_frequencies, npts)
    lm.0 <- lm(time_series[,2] ~ xm)

    # calculate rho and sigma based on residuals
    rho = cor(lm.0$residuals[1:(npts-1)],
              lm.0$residuals[2:npts])

    sigma = sd(lm.0$residuals)

    logLH = calcLogLH(lm.0$residuals, rho, sigma)

    return(cbind(sed_rate, logLH, rho, sigma))
  }

  # function to perform log-likelihood calculation, including
  #  assessment of correlated residuals, assuming AR1 model
  #  npts passed into function transparently
  calcLogLH <- function(residuals, rho, sigma)
  {
    # calculate Re^-1 (ReInv), as in EQ. A-7 of Malinverno & Briggs (2004)
    # set up array, with zeros
    ReInv <- double(npts * npts)
    dim(ReInv) <-c(npts, npts)
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
  ans <- rep(NA, 7)
  dim(ans) <- c(7)

  #######################################################################################
  # begin sedimentation rate loop

  i = 0
  # CALIBRATE DEPTH SERIES (m) TO TIME (ka)
  ts = dat
  # create new time vector

  # it is the index vector for time
  it <- seq(1, npts, by=1)
  time = (dx / sed_rate) * (it - 1)
  ts[1] = time

  # bandpass precession or short eccentricity band
  bp = taner(ts,
             padfac=2,
             flow=flow,
             fhigh=fhigh,
             roll=roll,
             demean=T,
             detrend=F,
             addmean=F,
             genplot=F,
             verbose=F)

  # hilbert transform for instantaneous amplitude
  hil = hilbert(bp,
                padfac=2,
                demean=T,
                detrend=F,
                addmean=F,
                genplot=F,
                verbose=F)
  # standardize hil to unit variance

  hil[2] = hil[2] - colMeans(hil[2])
  hil[2] = hil[2] / sd(hil[,2])

  # execute functions
  # for precession modulations
  if(fit == 1) {
    res = fitIt(sed_rate, hil, targetE)
    pwrOut = fitIt(sed_rate, ts , targetTot)
  }

  # for short eccentricity modulations
  if(fit == 2) {
    res = fitIt(sed_rate, hil, targetE[1])
    pwrOut = fitIt(sed_rate, ts, targetE)
  }

  ans[1] <- res[1]
  ans[2] <- res[2]
  ans[3] <- pwrOut[2]
  ans[4] <- res[3]
  ans[5] <- res[4]
  ans[6] <- pwrOut[3]
  ans[7] <- pwrOut[4]

  rPwr = ans[2] + ans[3]
    output <- data.frame('sed_rate' = 100 * ans[1],
                               'likelihood_env' = ans[2],
                               'likelihood_spec' = ans[3],
                               'likelihood_combined' = rPwr,
                               'env_rho' = ans[4],
                               'env_sigma' = ans[5],
                               'pwr_rho' = ans[6],
                               'pwr_sigma' = ans[7])

    return(output)
}
