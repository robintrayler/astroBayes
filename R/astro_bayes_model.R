#' Run a modified version of the algorithm from Malinverno et al., 2010 with an additional
#' modified to anchor the resulting age model using radioisotopic (e.g., U-Pb Ar/Ar) geochronology
#'
#' @name astro_bayes_model
#'
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

astro_bayes_model <- function(geochron_data,
                              cyclostrat_data,
                              tuning_frequency,
                              segment_edges,
                              sed_prior_range,
                              iterations = 100000,
                              burn = 5000) {
  # Need a whole bunch of error checking here ---------------------------------

  # check to make sure things are in order
  geochron_data   <- geochron_data %>% arrange(position)
  cyclostrat_data <- cyclostrat_data %>% arrange(position)

  # define interpolation grid -------------------------------------------------
  position_grid <- seq(from = min(segment_edges$position),
                       to   = max(segment_edges$position),
                       length = 1000)

  # set up model storage ------------------------------------------------------
  # store the sedimentation rate for each segment
  sed_rate <- sed_rate <- matrix(nrow = iterations,
                                 ncol = nrow(segment_edges) - 1)

  master_edges <- segment_edges
  master_geochron <- geochron_data

  # store all probabilities and anchor point
  anchor_point <- vector(length = iterations)

  # storage for the age model
  model_storage <- matrix(nrow = length(position_grid),
                          ncol = iterations)

  # set initial conditions ----------------------------------------------------
  # set anchor point using a simple linear regression
  anchor_point[1] <- geochron_data %>%
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
                         sd = 0.001)
  rm(mean_rate)

  # anchor the initial model in time
  anchored_model <- anchor_sed_model(segment_edges$position,
                                     sed_rate[1, ],
                                     anchor_point[1])

  # interpolate and store the model -------------------------------------------
  f <- approxfun(x = segment_edges$position,
                 y = anchored_model)
  model_storage[, 1] <- f(position_grid)
  rm(f)

  # run the MCMC model ----------------------------------------------------------
  pb <- progress::progress_bar$new(total = iterations,
                                   format = '[:bar] :percent eta: :eta')

  for(j in 2:iterations) {
    # update progress bar
    pb$tick()
    # store the previous iterations in case things get rejected ---------------
    sed_rate[j, ] <- sed_rate[j - 1, ]
    model_storage[, j] <- model_storage[, j - 1]
    anchor_point[j] <- anchor_point[j - 1]
    # update segment edges ----------------------------------------------------
    segment_edges$position <- runif(nrow(master_edges),
                                    min = master_edges$position - master_edges$thickness,
                                    max = master_edges$position + master_edges$thickness)
    # update geochron positions -----------------------------------------------
    geochron_data$position = runif(nrow(master_geochron),
                                   min = master_geochron$position - master_geochron$thickness / 2,
                                   max = master_geochron$position + master_geochron$thickness / 2)
    # update sed rates --------------------------------------------------------
    sed_rate[j, ] <- sed_rate[j - 1, ]
    for(q in 1:ncol(sed_rate)) {
      # propose a new rate ------------------------------------------
      proposed_rate <- adaptive_update(chain = sed_rate[, q],
                                       i = j,
                                       start_index = burn/2,
                                       initial_Cd = 0.001,
                                       distribution = 'gamma')

      # calculate the probability -----------------------------------
      proposed_prob <- pgram_likelihood(sed_rate = proposed_rate,
                                        segment_edges = segment_edges$position[q:(q + 1)],
                                        cyclostrat = cyclostrat_data,
                                        tuning_frequency = tuning_frequency$frequency)

      prior_proposed <- sed_prior(proposed_rate,
                                  min(sed_prior_range),
                                  max(sed_prior_range))
      current_prob  <- pgram_likelihood(sed_rate = sed_rate[j - 1, q],
                                        segment_edges = segment_edges$position[q:(q + 1)],
                                        cyclostrat = cyclostrat_data,
                                        tuning_frequency = tuning_frequency$frequency)
      prior_current <- sed_prior(sed_rate[j - 1, q],
                                 min(sed_prior_range),
                                 max(sed_prior_range))

      a <- (proposed_prob + prior_proposed) - (current_prob + prior_current)

      if(!is.na(a)) {
        if(!is.infinite(a)) {
          if(exp(a) > runif(1)) {
            sed_rate[j, q] <- proposed_rate
          }
        }
      }
    }

    # anchor the model in time --------------------------------------------------
    new_anchor <- adaptive_update(chain = anchor_point,
                                  i = j,
                                  start_index = burn/2,
                                  initial_Cd = 0.0001,
                                  distribution = 'gaussian')

    model_proposed <- anchor_sed_model(segment_edges$position,
                                       sed_rate[j, ],
                                       new_anchor)

    model_current <- anchor_sed_model(segment_edges$position,
                                      sed_rate[j, ],
                                      anchor_point[j-1])
    # calculate radiometric probability
    radio_proposed <- radio_likelihood(segment_edges = segment_edges$position,
                                       anchored_model = model_proposed,
                                       age = geochron_data$age,
                                       age_sd = geochron_data$age_sd,
                                       position = geochron_data$position)

    radio_current <- radio_likelihood(segment_edges$position,
                                      model_current,
                                      geochron_data$age,
                                      geochron_data$age_sd,
                                      position = geochron_data$position)

    if(exp(radio_proposed - radio_current) > runif(1)) {
      anchor_point[j] <- new_anchor
      f <- approxfun(x = segment_edges$position,
                     y = model_proposed)
      model_storage[, j] <- f(position_grid)
    }
  }

  # clean up the results and organize for output ------------------------------
  # calculate credible interval minus burn in
  CI <- apply(X = model_storage[, burn:iterations],
              MARGIN = 1,
              FUN = quantile,
              prob = c(0.025, 0.5, 0.975),
              na.rm = TRUE) %>%
    t() %>%
    data.frame()

  names(CI) <- c('CI_2.5', 'median', 'CI_97.5')

  CI <- CI %>% add_column(grid = position_grid)

  output = list(CI = CI,
                anchor_point = anchor_point,
                iterations = iterations,
                burn = burn,
                model_iterations = model_storage,
                sed_rate = sed_rate,
                segment_edges = master_edges,
                geochron_data = master_geochron,
                cyclostrat_data = cyclostrat_data,
                sed_prior_range = sed_prior_range,
                tuning_frequency = tuning_frequency)
  class(output) <- "astroBayesModel"
return(output)
}
