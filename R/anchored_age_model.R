#' Run a modified version of the algorithm from Malinverno et al., 2010 with an additional
#' modified to anchor the resulting age model using radioisotopic (e.g., U-Pb Ar/Ar) geochronology
#' @param geochron_data data frame containing radiometric age determinations, their uncertainties and stratigraphic positions, and a identifier variable.
#' column headers *must* be named exactly as follows:
#' * age: the weighted mean age in Ma
#' * age_sd: 1 sigma uncertainty in Ma
#' * position: stratigraphic position in depth coordinates (i.e. the top is 0)
#' * id: sample name. currently only one sample per unique id is allowed.
#'
#' @param cyclostrat_data a data frame where the first column is the stratigraphic position
#' and the second column is the cyclostratigraphic proxy record.
#' column headers *must* be named exactly as follows:
#' * position: stratigraphic position in depth coordinates (i.e. the top is 0).
#' * value: the measurement value of the proxy record.
#'
#' @param tuning_frequency data frame of tuning frequencies to use. Currently frequencies must be in cycles / Ma (i.e., long eccentricity ~ 1 / 0.405).
#' column headers *must* be named exactly as follows:
#' * frequency: the tuning frequency or frequencies to use
#  * orbital_cycle: character string name of each orbital cycle. Multiple rows with the same name are allowed.
#'
#' @param segment_edges stratigraphic points where sedimentation rate changes can  occur. Must be in the same stratigraphic scheme as geochron_data and cyclostrat_data.
#'
#' @param sed_prior_range vector of 2 sedimentation rates that define the minimum
#' and maximum allowable sedimentation rates. These two rates define a uniform
#' prior distribution for sedimentation rate.
#'
#' @param sed_mhsd metropolis hasting standard deviation for sedimentation rate
#' markov chains in depth / Ma. Defaults to 0.05 m/Ma. Someday this will be adaptive
#' but for now it must be tuned by hand.
#'
#' @param anchor_mhsd metropolis hasting standard deviation for model anchor age
#' markov chains in Ma. Defaults to 0.01 Ma. Someday this will be adaptive
#' but for now it must be tuned by hand.
#'
#' @param iterations how many Markov Chain Monte Carlo iterations should the model run for
#'
#' @param burn how many initial iterations to toss when calculating credible intervals
#'
#' @import "tidyverse"
#' @import "dplyr"
#' @import "astrochron"
#' @import "tibble"
#' @importFrom magrittr "%>%"
#'
#' @return a list with a whole bunch of stuff in it.
#' @md
#' @export
#'
anchored_age_model <- function(geochron_data,
                               cyclostrat_data,
                               tuning_frequency,
                               segment_edges,
                               sed_prior_range,
                               sed_mhsd = rep(0.05, nrow(segment_edges) - 1),
                               anchor_mhsd = 0.01,
                               iterations = 100000,
                               burn = 5000) {
  # Need a whole bunch of error checking here ---------------------------------

  # check to make sure columns are named correctly
  if ((colnames(geochron_data) %in%
       c("id", "age", "age_sd", "position", "thickness")) %>%
      all()) {cat('geochron_data column names correct')} else {
        cat('ERROR: column names for geochron_data must be *exactly*:
          \n "id", "age", "age_sd", "position", "thickness")')
        stop("**** TERMINATING NOW")
      }

  # check to make sure things are in order
  geochron_data <- geochron_data %>% arrange(position)
  cyclostrat_data <- cyclostrat_data %>% arrange(position)

  # define interpolation grid -------------------------------------------------
  position_grid <- seq(from = min(c(cyclostrat_data$position, geochron_data$position)),
                       to = max(c(cyclostrat_data$position, geochron_data$position)),
                       length = 5000)

  # set up model storage ------------------------------------------------------
  # store the sedimentation rate for each segment
  sed_rate <- sed_rate <- matrix(nrow = iterations,
                                 ncol = nrow(segment_edges) - 1)

  master_edges <- segment_edges
  # store all probabilities and anchor point
  parameter_storage <- data.frame(
    anchor_point = vector(length = iterations),
    pgram_LL     = vector(length = iterations),
    sed_prior_LL = vector(length = iterations),
    radio_LL     = vector(length = iterations),
    bayes_LL     = vector(length = iterations))

  # storage for the age model
  model_storage <- matrix(nrow = length(position_grid),
                          ncol = iterations)

  # set initial conditions ----------------------------------------------------
  # set anchor point using a simple linear regression
  parameter_storage$anchor_point[1] <- geochron_data %>%
    lm(age ~ position, data = .) %>%
    predict(newdata = data.frame(position = min(position_grid))) %>%
    as.vector()

  # set the starting sed_rate around the average rate
  mean_rate <- geochron_data %>%
    lm(position ~ age, data = .) %>%
    coef() %>%
    dplyr::nth(2)

  # randomly adjust starting rates
  sed_rate[1, ] <- rnorm(nrow(segment_edges) - 1,
                         mean = mean_rate,
                         sd = sed_mhsd)
  rm(mean_rate)

  # anchor the initial model in time
  anchored_model <- anchor_sed_model(segment_edges$position,
                                     sed_rate[1, ],
                                     parameter_storage$anchor_point[1])

  # calculate initial probabilities -------------------------------------------
  # sed rate likelihood
  parameter_storage$pgram_LL[1]     <- pgram_likelihood(sed_rate[1, ],
                                                        segment_edges$position,
                                                        cyclostrat_data,
                                                        tuning_frequency$frequency)
  # sed rate prior
  parameter_storage$sed_prior_LL[1] <- sed_prior(sed_rate[1, ],
                                                 min(sed_prior_range),
                                                 max(sed_prior_range))

  # calculate radiometric probability
  parameter_storage$radio_LL[1]     <- radio_likelihood(segment_edges$position,
                                                        anchored_model,
                                                        geochron_data$age,
                                                        geochron_data$age_sd,
                                                        runif(nrow(geochron_data),
                                                              min = geochron_data$position - geochron_data$thickness / 2,
                                                              max = geochron_data$position + geochron_data$thickness / 2))
  # calculate joint probability
  parameter_storage$bayes_LL[1]     <- parameter_storage$pgram_LL[1] +
    parameter_storage$sed_prior_LL[1] +
    parameter_storage$radio_LL[1]

  # interpolate and store the model -------------------------------------------
  f <- approxfun(x = segment_edges$position,
                 y = anchored_model)
  model_storage[, 1] <- f(position_grid)
  rm(f)

  # run the MCMC model ----------------------------------------------------------
  pb <- progress::progress_bar$new(total = iterations,
                                   format = '[:bar] :percent eta: :eta')


  for(j in 2:iterations){
    # update progress bar
    pb$tick()

    # store the previous iterations in case things get rejected ---------------
    sed_rate[j, ] <- sed_rate[j - 1, ]
    model_storage[, j] <- model_storage[, j - 1]
    parameter_storage[j, ] <- parameter_storage[j - 1, ]

    # propose new sed rates. use a gamma dist. to keep it from going negative
      new_rates <- rgamma(ncol(sed_rate),
                          shape = sed_rate[j - 1, ] ^ 2 / sed_mhsd ^ 2,
                          rate  = sed_rate[j - 1, ]  / sed_mhsd ^ 2)

    # adjust anchor point. No need to worry about this one going negative
    new_anchor <- rnorm(1,
                        mean = parameter_storage$anchor_point[j - 1],
                        sd = anchor_mhsd)

    # update segment edges ----------------------------------------------------
    segment_edges$position <- runif(nrow(master_edges),
                                    min = master_edges$position - master_edges$thickness,
                                    max = master_edges$position + master_edges$thickness)

    # anchor the model in time
    anchored_model <- anchor_sed_model(segment_edges$position,
                                       new_rates,
                                       new_anchor)

    # calculate probabilities of new model ------------------------------------
    pgram_proposed <- pgram_likelihood(new_rates,
                                       segment_edges$position,
                                       cyclostrat_data,
                                       tuning_frequency$frequency)
    # sed rate prior
    sed_prior_proposed <- sed_prior(new_rates,
                                    min(sed_prior_range),
                                    max(sed_prior_range))

    # calculate radiometric probability
    radio_proposed <- radio_likelihood(segment_edges$position,
                                       anchored_model,
                                       geochron_data$age,
                                       geochron_data$age_sd,
                                       runif(nrow(geochron_data),
                                             min = geochron_data$position - geochron_data$thickness / 2,
                                             max = geochron_data$position + geochron_data$thickness / 2))

    # calculate joint probability
    bayes_LL <- pgram_proposed + sed_prior_proposed + radio_proposed

    # keep it from going to NA or Inf
    bayes_LL <- ifelse(test = is.na(bayes_LL) | is.infinite(bayes_LL),
                       yes  = -1000,
                       no   = bayes_LL)


    # use a metropolis hasting algorithm to accept or reject
    if(exp(bayes_LL - parameter_storage$bayes_LL[j - 1]) > runif(1)){
      parameter_storage$pgram_LL[j] <- pgram_proposed
      parameter_storage$sed_prior_LL[j] <- sed_prior_proposed
      parameter_storage$radio_LL[j] <- radio_proposed
      parameter_storage$bayes_LL[j] <- bayes_LL
      parameter_storage$anchor_point[j] <- new_anchor
      sed_rate[j, ]   <- new_rates
      # segment_storage[j, ]   <- segment_edges$position

      # interpolate the model onto a grid
      f <- approxfun(x = segment_edges$position,
                     y = anchored_model)
      model_storage[, j] <- f(position_grid)
    }
  }

  # clean up the results and organize for output ------------------------------

  # calculate credible interval minus burn in
  CI <- apply(X = model_storage[, burn:iterations],
              MARGIN = 1,
              FUN = quantile,
              prob = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    t() %>%
    data.frame()
  names(CI) <- c('CI_2.5', 'median', 'CI_97.5')

  CI <- CI %>% add_column(grid = position_grid)


  output = list(CI = CI,
                parameters = parameter_storage,
                iterations = iterations,
                burn = burn,
                model_iterations = model_storage,
                sed_rate = sed_rate,
                segment_edges = master_edges,
                geochron_data = geochron_data,
                cyclostrat_data = cyclostrat_data,
                sed_prior_range = sed_prior_range,
                tuning_frequency = tuning_frequency)
  return(output)
}

