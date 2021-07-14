##-----------------------------------------------------------------------------
## this function is an implementation of the adaptive proposal algorithm of
## Haario et al. (1999). As implemented it only works for 1-dimensional
## sampling.
## ROBIN B. TRAYLER APRIL 26, 2019
##-----------------------------------------------------------------------------
## Inputs: chain = the current Markov Chain for a given parameter
##             h = how many previous iterations of the chain to use to calculate
##                 a proposal sd
##        update = how often should the mhsd update
##            cd = 2.4, empirically derived scaling factor. See Gelman 1996
##  mhsd_current = current mhsd
## Outputs mhsd = new proposal standard deviation

##-----------------------------------------------------------------------------
adaptive_proposal <- function(i,
                             chain,
                             h = 200,
                             update = 200,
                             cd = 2.4,
                             mhsd_current){
  # calculate the length of the chain
  # calculate the mean of the previous h iterations in the chain
  if(i%%update == 0 & i > update){
    K <- chain[(i - h + 1): i] - mean(chain[(i - h + 1): i])
    # calculate using the haario formula
    mhsd <- sqrt(var(K) * cd ^ 2)
    # set the proposal sd to the sd of the whole chain if it goes to 0 or NA
    mhsd <- ifelse(mhsd == 0, cd * sd(chain, na.rm = T), mhsd)
    mhsd <- ifelse(is.na(mhsd), cd * sd(chain, na.rm = T), mhsd)
  }
  else{
    mhsd <- mhsd_current
  }
  # return a new proposal
  return(mhsd)
}
