### This function is in development for astrochron: An R Package for Astrochronology
### Copyright (C) 2020 Stephen R. Meyers
###
###########################################################################
### function timeOptL - (SRM: April 12-14, 2017; September 23-24, 2020)
###
### modified from FUNCTION-timeOptL_v9.R, to perform timeOpt using formal
### likelihood functions
###########################################################################


time_opt_likelihood <- function (cyclostrat_data,
                                 tuning_frequency,
                                 sed_rate = 1,
                                 fitModPwr=T,
                                 flow=NULL,
                                 fhigh=NULL,
                                 roll=NULL)
{

  # prepare data array
  npts <- length(cyclostrat_data[, 1])
  dx <- cyclostrat_data[2, 1] - cyclostrat_data[1, 1]

  # standardize data series
  cyclostrat_data[2]=cyclostrat_data[2] - colMeans(cyclostrat_data[2])
  cyclostrat_data[2]=cyclostrat_data[2] / sapply(cyclostrat_data[2],sd)

  # convert sed_rate cm/ka to m/ka for processing
  sed_rate = sed_rate / 100

  #######################################################################################
  # set up default bandpass frequencies and targets
  #  first for precession
  if(is.null(flow)) {flow = 0.035}
  if(is.null(fhigh)) {fhigh = 0.065}
  if(is.null(roll)) {roll = 10 ^ 3}

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
      storage[,2 * i-1] <- cos( (2*pi) /
                                  (target_frequencies[i]) * (dx / sed_rate) * (1:n))

      storage[,2 * i] <- sin( (2*pi) /
                                (target_frequencies[i]) * (dx / sed_rate) * (1:n))
    }
    return(storage)
  }

  # function to perform fitting
  #  dx, npts passed into function transparently
  fitIt <- function(sed_rate,
                    time_series,
                    target_frequencies) {
    xm <- genCycles(sed_rate, target_frequencies, npts)
    lm.0 <- lm(time_series[, 2] ~ xm)

    # calculate rho and sigma based on residuals
    rho = cor(lm.0$residuals[1:(npts - 1)],
              lm.0$residuals[2:npts])

    sigma = sd(lm.0$residuals)

    # calculate log likelihood
    logLH = calcLogLH(lm.0$residuals, rho, sigma)

    # return everything
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
    logLH2 = -0.5 * logLH2
    return(logLH2)
  }

  ans <- vector(mode = 'numeric', length = 7)

  # CALIBRATE DEPTH SERIES (m) TO TIME (ka)
  time_series = cyclostrat_data

  # create new time vector
  # it is the index vector for time
  it <- seq(1, npts, by=1)
  time = (dx / sed_rate) * (it - 1)
  time_series[1] = time

  # bandpass precession or short eccentricity band
  bp = taner(time_series,
             padfac = 2,
             flow = flow,
             fhigh = fhigh,
             roll = roll,
             demean = T,
             detrend = F,
             addmean = F,
             genplot = F,
             verbose = F)

  # hilbert transform for instantaneous amplitude
  hil = hilbert(bp,
                padfac = 2,
                demean = T,
                detrend = F,
                addmean = F,
                genplot = F,
                verbose = F)

  # standardize hil to unit variance
  hil[2] = hil[2] - colMeans(hil[2])
  hil[2] = hil[2] / sd(hil[, 2])

  # execute functions
  # for precession modulations
  res = fitIt(sed_rate, hil, tuning_frequency)
  pwrOut = fitIt(sed_rate, time_series , tuning_frequency)

  output <- data.frame('sed_rate' = 100 * res[1],
                       'likelihood_env' = res[2],
                       'likelihood_spec' = pwrOut[2],
                       'env_rho' = res[3],
                       'env_sigma' = res[4],
                       'pwr_rho' = pwrOut[3],
                       'pwr_sigma' = pwrOut[4]) %>%
    mutate(likelihood_combined = likelihood_env + likelihood_spec)

  return(output)
}
