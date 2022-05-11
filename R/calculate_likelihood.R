# this is a wrapper function to calculate model likelihoods
# there are two options `malinverno` and `time_opt`
#

calculate_likelihood <- function(sed_rate,
                                 segment_edges,
                                 cyclostrat_data,
                                 tuning_frequency,
                                 method = NA){

  # preallocate storage -------------------------------------------------------
  LL <- vector(length = length(sed_rate))

  for(i in 1:length(sed_rate)) {
    # loop through all segments an calculate scaled periodogram ---------------

    current_cyclostrat <- cyclostrat_data %>%
      filter(position > segment_edges[i] &
               position < segment_edges[i + 1])

    if(method == 'malinverno') {
      LL[i] <- malinverno_likelihood(cyclostrat_data = current_cyclostrat,
                                     sed_rate = sed_rate[i],
                                     tuning_frequency = tuning_frequency)
    }

    if(method == 'time_opt') {
      LL[i] <- time_opt_likelihood(cyclostrat_data = current_cyclostrat,
                                   sed_rate = sed_rate[i],
                                   tuning_frequency = tuning_frequency)$LL_total
    }
  }

  # calculate joint probability
  LL <- sum(LL)

  return(LL)
}


