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
#' @param target_frequency data frame of target frequencies to use. Frequencies must be in cycles / Ma (i.e., long eccentricity ~ 1 / 0.405).
#' The column headers *must* be named exactly as follows:
#'
#'  * frequency: the target frequency or frequencies to use
#'  * orbital_cycle: character string name of each orbital cycle.
#'    * Multiple rows with the same name are allowed.
#'
#' @param layer_boundaries A data frame specifying where sedimentation rate changes can occur. Must be in the same stratigraphic scheme as geochron_data and cyclostrat_data. The column headers *must* be named exactly as follows
#'  * `position`: Stratigraphic position in depth coordinates (i.e. the top is 0). Note that the top and bottom of the section must be specified.
#'  * `thickness`: Uniform uncertainty of `position`.
#'  * `hiatus_boundary`: `TRUE` or `FALSE`. Should the layer boundary be treated as a hiatus?
#'  * `sed_min`: minimum sedimentation rate boundary for the layer.
#'  * `sed_max`: minimum sedimentation rate boundary for the layer.
#'       * Note: While `sed_min` and `sed_max` are specified in the same rows as the layer boundaries, they are applied to the layer between *their own* row and the row below.
#'
#' @param iterations How many Markov Chain Monte Carlo iterations should the model run for
#'
#' @param burn How many initial iterations to toss when calculating credible intervals
#'
#' @param method Probability calculation method. Currently only `malinverno` is supported.
#'
#' @import "tidyverse"
#' @import "tidyr"
#' @import "dplyr"
#' @import "astrochron"
#' @import "tibble"
#' @importFrom magrittr "%>%"
#'
#' @return The function returns a list of class `astroBayesModel`
#' which contains the following objects:
#'  * `CI` data frame containing the credible interval for the age model
#'  * `anchor_point` a data.frame containing the posterior sample for the anchor_point
#'    parameter(s)
#'  * `cyclostrat_CI` data frame containg the credible interval for the cyclostratigraphic data.
#'  * `sed_rate` a matrix containing the posterior sample of sedimentation rate
#'  for each model layer.
#'  * `model_iterations` a matrix containing the individual model iterations used
#'  to calculate `CI`. These can be plotted against the `position` column in `CI`
#'  to visualize.
#'
#'  The output object also includes all the model inputs: shown below. These are
#'  carried through for ease of data visualization.
#'  * `layer_boundaries`
#'  * `geochron_data`
#'  * `cyclostrat_data`
#'  * `target_frequency`
#'  * `hiatus_durations`
#'
#' @md
#' @export

astro_bayes_model <- function(geochron_data,
                              cyclostrat_data,
                              target_frequency,
                              layer_boundaries,
                              iterations = 10000,
                              burn = 5000,
                              method = 'malinverno') {
  # error checking ------------------------------------------------------------
  # check to see if multiple cyclostratigrapic records are in use
  message('initial error checking...')
  if(class(cyclostrat_data) == 'data.frame') {
    cycle_list <- list()
    cycle_list[[1]] <- cyclostrat_data
    cyclostrat_data <- cycle_list
    rm(cycle_list)
  }

  # geochron data
  if(!all(c('id', 'age', 'age_sd', 'position', 'thickness') %in%
          names(geochron_data))
  ) {stop('geochron_data columns are named incorrectly. \n
          Column names must be exactly: "id", "age", "age_sd", "position", "thickness"')}

  # layer_boundaries
  if(!all(c('position', 'thickness', 'hiatus_boundary', 'sed_min', 'sed_max') %in% names(layer_boundaries))
  ) {stop('layer_boundaries columns are named incorrectly. \n
          Column names must be exactly: "position", "thickness", "hiatus_boundary", "sed_min", "sed_max"')}

  # cyclostrat data
  # if(!all(names(cyclostrat_data) %in% c('position', 'value'))
  # ) {stop('cyclostrat_data columns are named incorrectly. \n
  #         Column names must be exactly: "position", "value"')}

  # target_frequency
  if(!all(c('frequency', 'orbital_cycle') %in% names(target_frequency))
  ) {stop('target_frequency columns are named incorrectly. \n
          Column names must be exactly: "frequency", "orbital_cycle"')}

  # check for missing values --------------------------------------------------
  # geochron_data
  if(!(nrow(geochron_data) > 0)) {
    stop('geochron_data has 0 rows. It must contain at least two rows')
  }

  #  layer_boundaries
  if(!(nrow(layer_boundaries) > 1)) {
    stop('layer_boundaries has less than 2 rows. It must contain at least two rows')
  }

  # cyclostrat data
  if(!any(unlist(lapply(cyclostrat_data, nrow)) > 0)) {
    stop('cyclostrat_data has 0 rows.')
  }

  # target frequency
  if(!(nrow(target_frequency) > 0)) {
    stop('target_frequency has 0 rows. It must contain at least 1 row')
  }

  # use malinverno probability if not specified
  if(is.na(method)) {
    method = 'malinverno'
    message('no method specified, defaulting to "malinverno" spectral power method.')
  }

  # check to make sure things are in order ------------------------------------
  geochron_data   <- geochron_data %>% arrange(position)

  for(i in seq_along(cyclostrat_data)) {
    cyclostrat_data[[i]] <- cyclostrat_data[[i]] %>% arrange(position)
  }
  rm(i)

  message('finished')
  # store the input data for when positions move around later
  master_edges <- layer_boundaries
  master_geochron <- geochron_data

  # define interpolation grid -------------------------------------------------
  position_grid <- seq(from = min(layer_boundaries$position),
                       to   = max(layer_boundaries$position),
                       length = 1000)

  # set up model storage ------------------------------------------------------
  # store the sedimentation rate for each layer
  sed_rate <- matrix(nrow = iterations,
                     ncol = nrow(layer_boundaries) - 1)

  layer_storage <- matrix(nrow = iterations,
                            ncol = nrow(layer_boundaries))

  layer_storage[1, ] <- layer_boundaries$position

  # store all probabilities and anchor point
  anchor_point <- vector(length = iterations)

  # hiatus duration
  n_hiatus <- sum(layer_boundaries$hiatus_boundary)
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

  tuned_cyclostrat <- list()
  for(i in seq_along(cyclostrat_data)) {
    tuned_cyclostrat[[i]] <- matrix(nrow = nrow(cyclostrat_data[[i]]),
                                    ncol = iterations)
  }

  # set initial conditions ----------------------------------------------------
  # set anchor point using a simple linear regression
  anchor_point[1] <- geochron_data %>%
    lm(age ~ position, data = .) %>%
    predict(newdata = data.frame(position = min(position_grid))) %>%
    as.vector()

  message('Setting initial sedimentation rates to maximum likelihood... \n (this can take a few minutes)')
  # set the starting rates to the maximum likelihood
  sed_rate[1, ] <- visualize_likelihood(layer_boundaries  = layer_boundaries,
                                         cyclostrat_data  = cyclostrat_data,
                                        target_frequency = target_frequency,
                                         resolution = 0.1,
                                         method = method,
                                        plot = FALSE) |>
    group_by(name) |>
    filter(probability == max(probability, na.rm = TRUE)) |>
    pull(sed_rate) |>
    rev()

  # anchor the initial model in time ------------------------------------------
  model_storage[, 1] <- anchor_sed_model(layer_boundaries = layer_boundaries,
                                         sed_rate[1, ],
                                         anchor_point = anchor_point[1],
                                         position_grid = position_grid,
                                         hiatus_duration =
                                           if(n_hiatus > 0) {
                                             hiatus_storage[1, ]
                                           } else {
                                             NA
                                           }
  )

  # run the MCMC model --------------------------------------------------------
  message('Staring the MCMC model...')
  pb <- progress::progress_bar$new(total = iterations,
                                   format = '[:bar] :percent eta: :eta')

  for(j in 2:iterations) {
    # update progress bar
    pb$tick()

    # store the previous iterations in case things get rejected ---------------
    sed_rate[j, ] <- sed_rate[j - 1, ]
    model_storage[, j] <- model_storage[, j - 1]
    anchor_point[j] <- anchor_point[j - 1]

    # if there are hiatuses ---------------------------------------------------
    if(n_hiatus > 0) {
      hiatus_storage[j, ] <- hiatus_storage[j - 1, ]
    }

    # randomly adjust layer edges ---------------------------------------------
    layer_boundaries$position <- runif(nrow(master_edges),
                                    min = master_edges$position -
                                      master_edges$thickness / 2,
                                    max = master_edges$position +
                                      master_edges$thickness / 2)

    layer_storage[j, ] <- layer_boundaries$position
    # randomly geochronology positions ----------------------------------------
    geochron_data$position = runif(nrow(master_geochron),
                                   min = master_geochron$position -
                                     master_geochron$thickness / 2,
                                   max = master_geochron$position +
                                     master_geochron$thickness / 2)

    # update sedimentation rates ----------------------------------------------
    for(q in 1:ncol(sed_rate)) {
      # get current sed rates
      current_rates <- sed_rate[j - 1, ]
      # propose a new rate ----------------------------------------------------
      proposed_rate <- adaptive_update(chain = sed_rate[, q],
                                       i = j,
                                       start_index = burn / 2,
                                       initial_Cd = 0.01,
                                       lower = layer_boundaries$sed_min[q],
                                       upper = layer_boundaries$sed_max[q])

      # calculate the probability ---------------------------------------------
      proposed_prob <- vector()
      for(i in seq_along(cyclostrat_data)) {
        proposed_prob[i] <- calculate_likelihood(sed_rate = proposed_rate,
                                                 layer_boundaries = layer_boundaries$position[q:(q + 1)],
                                                 cyclostrat = cyclostrat_data[[i]],
                                                 target_frequency = target_frequency,
                                                 method = method)

      }
      proposed_prob <- sum(proposed_prob)

      current_prob <- vector()
      for(i in seq_along(cyclostrat_data)) {
        current_prob[i] <- calculate_likelihood(sed_rate = sed_rate[j - 1, q],
                                                layer_boundaries = layer_boundaries$position[q:(q + 1)],
                                                cyclostrat = cyclostrat_data[[i]],
                                                target_frequency = target_frequency,
                                                method = method)

      }
      current_prob <- sum(current_prob)

      # calculate the radioisotopic probability -------------------------------
      current_rates[q] <- proposed_rate
      model_proposed <- anchor_sed_model(layer_boundaries = layer_boundaries,
                                         sed_rate = current_rates,
                                         anchor_point = anchor_point[j - 1],
                                         hiatus_duration =
                                           if(n_hiatus > 0) {
                                             hiatus_storage[j - 1, ]
                                           } else {
                                             NA
                                           },
                                         position_grid = position_grid)

      radio_proposed <- radio_likelihood(anchored_model = model_proposed,
                                         position_grid = position_grid,
                                         age = geochron_data$age,
                                         age_sd = geochron_data$age_sd,
                                         position = geochron_data$position)

      model_current <- anchor_sed_model(layer_boundaries = layer_boundaries,
                                        sed_rate = sed_rate[j - 1, ],
                                        anchor_point = anchor_point[j - 1],
                                        hiatus_duration =
                                          if(n_hiatus > 0) {
                                            hiatus_storage[j - 1, ]
                                          } else {
                                            NA
                                          },
                                        position_grid = position_grid)

      radio_current <- radio_likelihood(anchored_model = model_current,
                                        position_grid = position_grid,
                                        age = geochron_data$age,
                                        age_sd = geochron_data$age_sd,
                                        position = geochron_data$position)

      # use a Metropolis-Hastings algorithm to accept or reject
      p <- (proposed_prob + radio_proposed) - (current_prob + radio_current)

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
    new_anchor <- adaptive_update(chain = anchor_point,
                                  i = j,
                                  start_index = burn / 2,
                                  initial_Cd = 0.001,
                                  lower = -1e10,
                                  upper = 1e10)

    # anchor the proposed and current models in absolute time -----------------
    model_proposed <- anchor_sed_model(layer_boundaries = layer_boundaries,
                                       sed_rate = sed_rate[j, ],
                                       anchor_point = new_anchor,
                                       hiatus_duration =
                                         if(n_hiatus > 0) {
                                           hiatus_storage[j - 1, ]
                                         } else {
                                           NA
                                         },
                                       position_grid = position_grid)

    model_current <- anchor_sed_model(layer_boundaries = layer_boundaries,
                                      sed_rate = sed_rate[j, ],
                                      anchor_point = anchor_point[j - 1],
                                      hiatus_duration =
                                        if(n_hiatus > 0) {
                                          hiatus_storage[j - 1, ]
                                        } else {
                                          NA
                                        },
                                      position_grid = position_grid)

    # calculate radiometric probability ---------------------------------------
    # only the anchor point changes at this point
    # hiatus duration are updated later
    radio_proposed <- radio_likelihood(anchored_model = model_proposed,
                                       position_grid = position_grid,
                                       age = geochron_data$age,
                                       age_sd = geochron_data$age_sd,
                                       position = geochron_data$position)

    radio_current <- radio_likelihood(anchored_model = model_current,
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
    # update hiatus duration --------------------------------------------------
    if(n_hiatus > 0) {
      duration <- vector() # temporary vector of duration

      for(k in 1:n_hiatus) {
        duration[k] <- adaptive_update(chain = hiatus_storage[, k],
                                       i = j,
                                       start_index = burn / 2,
                                       initial_Cd = 0.00001,
                                       lower = 0,
                                       upper = 1e10)
      }

      # anchor the model in absolute time ---------------------------------------
      model_proposed <- anchor_sed_model(layer_boundaries = layer_boundaries,
                                         sed_rate = sed_rate[j, ],
                                         anchor_point = anchor_point[j],
                                         hiatus_duration = duration,
                                         position_grid = position_grid)

      model_current <- anchor_sed_model(layer_boundaries = layer_boundaries,
                                        sed_rate = sed_rate[j, ],
                                        anchor_point = anchor_point[j],
                                        hiatus_duration =
                                          hiatus_storage[j - 1, ],
                                        position_grid = position_grid)

      # calculate radiometric probability ---------------------------------------
      # only the duration(s) change at this point
      # anchor point is constant
      duration_proposed <- radio_likelihood(anchored_model = model_proposed,
                                            position_grid = position_grid,
                                            age = geochron_data$age,
                                            age_sd = geochron_data$age_sd,
                                            position = geochron_data$position)

      duration_current  <- radio_likelihood(anchored_model = model_current,
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

    ####FIX THIS FOR LISTS####

    for(i in seq_along(cyclostrat_data)) {
      tuned_cyclostrat[[i]][, j] <- f(cyclostrat_data[[i]]$position)
    }
    # tuned_cyclostrat[, j] <- f(cyclostrat_data$position)

  } # end of the main loop ####################################################

  # clean up the results and organize for output ------------------------------
  # calculate credible interval minus burn-in
  CI <- apply(X = model_storage[, burn:iterations],
              MARGIN = 1,
              FUN = quantile,
              prob = c(0.025, 0.5, 0.975),
              na.rm = TRUE) %>%
    t() %>%
    data.frame()

  # make it ggplot friendly
  names(CI) <- c('CI_2.5', 'median', 'CI_97.5')
  CI <- CI %>% add_column(position = position_grid)

  # calculate the credible interval for the cyclostratigraphy
  cyclostrat_CI <- list()
  for(i in seq_along(cyclostrat_data)) {
    CI_temp <- apply(X = tuned_cyclostrat[[i]][,burn:iterations],
                     MARGIN = 1,
                     FUN = quantile,
                     prob = c(0.025, 0.5, 0.975),
                     na.rm = TRUE)  %>%
      t() %>%
      as.data.frame()
    names(CI_temp) <- c('CI_2.5', 'median', 'CI_97.5')

    cyclostrat_CI[[i]] <- CI_temp %>%
      add_column(value = cyclostrat_data[[i]]$value)

  }

  # gather the inputs and outputs into a list
  output = list(CI = CI,                             # credible interval
                cyclostrat_CI = cyclostrat_CI,       # credible interval for cyclostrat data
                anchor_point = anchor_point,         # anchor point age
                iterations = iterations,             # number of iterations
                burn = burn,                         # burn-in
                model_iterations = model_storage,    # individual age models
                sed_rate = sed_rate,                 # sedimentation rate chains
                layer_boundaries = master_edges,     # layer_boundaries input
                layer_storage = layer_storage,     # layer boundary positions
                geochron_data = master_geochron,     # geochron input
                cyclostrat_data = cyclostrat_data,   # cyclostrat_data input
                target_frequency = target_frequency, # target frequencies input
                hiatus_durations = hiatus_storage)   # hiatus duration

  # assign a class and return
  class(output) <- "astroBayesModel"
  return(output)
}
