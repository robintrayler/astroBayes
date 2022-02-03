#' Run a modified version of the algorithm from Malinverno et al., 2010 with an additional
#' modified to anchor the resulting age model using radioisotopic (e.g., U-Pb Ar/Ar) geochronology
#'
#' @name astro_bayes_hiatus_model
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
#'
#' * frequency: the tuning frequency or frequencies to use
#'  * orbital_cycle: character string name of each orbital cycle.
#'   Multiple rows with the same name are allowed.
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
#' @return The function returns a list of class `astroBayesModel`
#' which contains the following objects:
#'  `CI` data frame containing the credible interval for the age model
#'  `anchor_point` a data.frame containing the posterior sample for the anchor_point
#'    parameter(s)
#'  * `sed_rate` a matrix containing the posterior sample of sedimentation rate
#'  for each model segment
#'  * `model_iterations` a matrix containing the individual model iterations used
#'  to calculate `CI`. These can be plotted against the `position` column in `CI`
#'  to visualize.
#'
#'  The output object also includes all the model inputs: shown below. These are
#'  carried through for ease of data visualization.
#'  * `segment_edges`
#'  * `geochron_data`
#'  * `cyclostrat_data`
#'  * `sed_prior_range`
#'  * `tuning_frequency`
#'
#' @md
#' @export

astro_bayes_hiatus_model <- function(geochron_data,
                                     cyclostrat_data,
                                     tuning_frequency,
                                     segment_edges,
                                     sed_prior_range = c(0, 1000),
                                     iterations = 10000,
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
  sed_rate <- matrix(nrow = iterations,
                     ncol = nrow(segment_edges) - 1)

  master_edges <- segment_edges
  master_geochron <- geochron_data

  # store all probabilities and anchor point
  anchor_point <- vector(length = iterations)

  # hiatus duration
  n_hiatus <- sum(segment_edges$hiatus_boundary)
  if(n_hiatus > 0) {
    hiatus_storage <- matrix(0,
                             nrow = iterations,
                             ncol = n_hiatus)
  } else {
    hiatus_storage <- NA
  }

  # storage for the age model
  model_storage <- matrix(nrow = length(position_grid),
                          ncol = iterations)

  tuned_cyclostrat <- matrix(nrow = nrow(cyclostrat_data),
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

  # anchor the initial model in time ------------------------------------------
  model_storage[, 1] <- anchor_hiatus_sed_model(segment_edges = segment_edges,
                                                sed_rate[1, ],
                                                anchor_point = anchor_point[1],
                                                position_grid = position_grid,
                                                hiatus_duration =
                                                  ifelse(
                                                    test = n_hiatus > 0,
                                                    yes = hiatus_storage[1, ],
                                                    no = NA)
  )

  # run the MCMC model --------------------------------------------------------
  pb <- progress::progress_bar$new(total = iterations,
                                   format = '[:bar] :percent eta: :eta')

  for(j in 2:iterations) {
    # update progress bar
    pb$tick()

    # store the previous iterations in case things get rejected ---------------
    sed_rate[j, ] <- sed_rate[j - 1, ]
    model_storage[, j] <- model_storage[, j - 1]
    anchor_point[j] <- anchor_point[j - 1]

    # if there are hiatuses
    if(n_hiatus > 0) {
      hiatus_storage[j, ] <- hiatus_storage[j - 1, ]
    }

    # randomly adjust segment edges -------------------------------------------
    segment_edges$position <- runif(nrow(master_edges),
                                    min = master_edges$position -
                                      master_edges$thickness / 2,
                                    max = master_edges$position +
                                      master_edges$thickness / 2)

    # randomly geochronology positions ----------------------------------------
    geochron_data$position = runif(nrow(master_geochron),
                                   min = master_geochron$position -
                                     master_geochron$thickness / 2,
                                   max = master_geochron$position +
                                     master_geochron$thickness / 2)

    # update sedimentation rates ----------------------------------------------
    for(q in 1:ncol(sed_rate)) {

      # propose a new rate ----------------------------------------------------
      proposed_rate <- adaptive_update_truncated(chain = sed_rate[, q],
                                                 i = j,
                                                 start_index = burn / 2,
                                                 initial_Cd = 0.001,
                                                 lower = sed_prior_range[1],
                                                 upper = sed_prior_range[2])

      # calculate the probability ---------------------------------------------
      proposed_prob <- pgram_likelihood(sed_rate = proposed_rate,
                                        segment_edges = segment_edges$position[q:(q + 1)],
                                        cyclostrat = cyclostrat_data,
                                        tuning_frequency = tuning_frequency$frequency)

      current_prob  <- pgram_likelihood(sed_rate = sed_rate[j - 1, q],
                                        segment_edges = segment_edges$position[q:(q + 1)],
                                        cyclostrat = cyclostrat_data,
                                        tuning_frequency = tuning_frequency$frequency)


      # use a Metropolis-Hastings algorithm to accept or reject
      p <- (proposed_prob) - (current_prob)

      if(!is.na(p)) {
        if(!is.infinite(p)) {
          if(exp(p) > runif(1)) {
            sed_rate[j, q] <- proposed_rate
          }
        }
      }
    }

    # Update the anchor point -------------------------------------------------
    # propose a new age for the anchor point
    new_anchor <- adaptive_update_truncated(chain = anchor_point,
                                            i = j,
                                            start_index = burn / 2,
                                            initial_Cd = 0.001,
                                            lower = -1e-10,
                                            upper = 1e10)


    # anchor the proposed and current models in absolute time -----------------
    model_proposed <- anchor_hiatus_sed_model(segment_edges = segment_edges,
                                              sed_rate = sed_rate[j, ],
                                              anchor_point = new_anchor,
                                              hiatus_duration = ifelse(
                                                test = n_hiatus > 0,
                                                yes = hiatus_storage[j - 1, ],
                                                no = NA
                                              ),
                                              position_grid = position_grid)

    model_current <- anchor_hiatus_sed_model(segment_edges = segment_edges,
                                             sed_rate = sed_rate[j, ],
                                             anchor_point = anchor_point[j - 1],
                                             hiatus_duration = ifelse(
                                               test = n_hiatus > 0,
                                               yes = hiatus_storage[j - 1, ],
                                               no = NA
                                             ),
                                             position_grid = position_grid)


    # calculate radiometric probability ---------------------------------------
    # only the anchor point changes at this point
    # hiatus duration are updated later
    radio_proposed <- hiatus_radio_likelihood(anchored_model = model_proposed,
                                              position_grid = position_grid,
                                              age = geochron_data$age,
                                              age_sd = geochron_data$age_sd,
                                              position = geochron_data$position)

    radio_current <- hiatus_radio_likelihood(anchored_model = model_current,
                                             position_grid = position_grid,
                                             age = geochron_data$age,
                                             age_sd = geochron_data$age_sd,
                                             position = geochron_data$position)

    # use a Metropolis-Hastings algorithm to accept or reject
    p <- radio_proposed - radio_current

    if(!is.na(p)) {
      if(!is.infinite(p)) {
        if(exp(p) > runif(1)) {
          anchor_point[j] <- new_anchor
          model_storage[, j] <- model_proposed
        }
      }
    }
    if(n_hiatus > 0) {
      # update hiatus duration --------------------------------------------------
      duration <- vector() # temporary vector of duration
      for(k in 1:n_hiatus) {
        duration[k] <- adaptive_update_truncated(chain = hiatus_storage[, k],
                                                 i = j,
                                                 start_index = burn/2,
                                                 initial_Cd = 0.00001,
                                                 lower = -1e-10,
                                                 upper = 1e10)
      }

      # anchor the model in absolute time ---------------------------------------
      model_proposed <- anchor_hiatus_sed_model(segment_edges = segment_edges,
                                                sed_rate = sed_rate[j, ],
                                                anchor_point = anchor_point[j],
                                                hiatus_duration = duration,
                                                position_grid = position_grid)

      model_current <- anchor_hiatus_sed_model(segment_edges = segment_edges,
                                               sed_rate = sed_rate[j, ],
                                               anchor_point = anchor_point[j],
                                               hiatus_duration = hiatus_storage[j - 1, ],
                                               position_grid = position_grid)


      # calculate radiometric probability ---------------------------------------
      # only the duration(s) change at this point
      # anchor point is constant
      duration_proposed <- hiatus_radio_likelihood(anchored_model = model_proposed,
                                                   position_grid = position_grid,
                                                   age = geochron_data$age,
                                                   age_sd = geochron_data$age_sd,
                                                   position = geochron_data$position)

      duration_current <- hiatus_radio_likelihood(anchored_model = model_current,
                                                  position_grid = position_grid,
                                                  age = geochron_data$age,
                                                  age_sd = geochron_data$age_sd,
                                                  position = geochron_data$position)

      # use a Metropolis-Hastings algorithm to accept or reject
      p <- duration_proposed - duration_current

      if(!is.na(p)) {
        if(!is.infinite(p)) {
          if(exp(p) > runif(1)) {
            hiatus_storage[j, ] <- duration
            model_storage[, j]  <- model_proposed

          }
        }
      }
    }
    # predict age of cyclostrat each iteration --------------------------------
    f <- approxfun(x = position_grid,
                   y = model_storage[, j])

    tuned_cyclostrat[, j] <- f(cyclostrat_data$position)
  } # end of the main loop

  # clean up the results and organize for output ------------------------------
  # calculate credible interval minus burn-in
  CI <- apply(X = model_storage[, burn:iterations],
              MARGIN = 1,
              FUN = quantile,
              prob = c(0.025, 0.5, 0.975),
              na.rm = TRUE) %>%
    t() %>%
    data.frame()

  cyclostrat_CI <- apply(X = tuned_cyclostrat[, burn:iterations],
                         MARGIN = 1,
                         FUN = quantile,
                         prob = c(0.025, 0.5, 0.975),
                         na.rm = TRUE) %>%
    t() %>%
    data.frame()

  # make it ggplot friendly
  names(CI) <- c('CI_2.5', 'median', 'CI_97.5')
  CI <- CI %>% add_column(position = position_grid)

  # make it ggplot friendly
  names(cyclostrat_CI) <- c('CI_2.5', 'median', 'CI_97.5')
  cyclostrat_CI <- cyclostrat_CI %>%
    add_column(value = cyclostrat_data$value)

  # gather the inputs and outputs into a list
  output = list(CI = CI, # credible interval
                cyclostrat_CI = cyclostrat_CI, # credible interval for cyclostrat data
                anchor_point = anchor_point, # anchor point age
                iterations = iterations, # number of iterations
                burn = burn, # burn-in
                model_iterations = model_storage, # individual age models
                sed_rate = sed_rate, # sedimentation rate chains
                segment_edges = master_edges, # segment_edges input
                geochron_data = master_geochron, # geochron input
                cyclostrat_data = cyclostrat_data, # cyclostrat_data input
                sed_prior_range = sed_prior_range, # sed prior input
                tuning_frequency = tuning_frequency, # tuning frequencies input
                hiatus_durations = hiatus_storage) # hiatus duration

  # assign a class and return
  class(output) <- "astroBayesModel"
  return(output)
}
