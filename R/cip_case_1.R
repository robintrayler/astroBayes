# process CIP Case 1 ----------------------------------------------------------
library(tidyverse)
library(astrochron)
library(astroBayes)

true_age <- read.csv(file = './data/CIP_data/case_1/true_age.csv')
tuning_frequency <- read.csv(file = './data/data_A/tuning_frequency.csv') %>%
  filter(orbital_cycle == 'precession')
rast <- raster::raster('./data/CIP_data/case_1/Signal_1_grey.tiff')

cyclostrat_data <- cbind(seq(0, 26.23, length = 1146), rast[,1]) %>%
  as.data.frame() %>%
  set_names(nm = c('position','value'))

cyclostrat_data <- bandpass(cyclostrat_data,
                            flow = 1,
                            fhigh = 7.5,
                            demean = TRUE,
                            detrend = TRUE,
                            genplot = FALSE)



n = 3 # number of points to generate
# pick some true ages
date_positions <- true_age[sample(seq_along(true_age[, 1]), n), ]
geochron_data  <- # assemble into synthetic geochronology
  data.frame(
    # set the mean age to the true age +-
    # case where accuracy reflects precision
    age = rnorm(n, mean = date_positions$age, sd = date_positions$age * 0.001),
    # set standard deviation to 0.1% of the age plus some noise
    age_sd = rnorm(n, date_positions$age * 0.001, sd = 0.0001),
    position = date_positions$position,
    thickness = 0,
    id = letters[1:n]) %>%
  arrange(position)


segment_edges <-
  data.frame(position = c(0, 5, 26.23),
             thickness = c(0,5,  0),
             hiatus_boundary = FALSE)

# run the model ---------------------------------------------------------------
model_output <- astro_bayes_model(geochron_data = geochron_data,
                                  cyclostrat_data = cyclostrat_data,
                                  tuning_frequency = tuning_frequency,
                                  segment_edges = segment_edges,
                                  iterations = 10000,
                                  burn = 2000,
                                  sed_prior_range = c(0, 20))

plot(model_output, type = 'age_depth') +
  geom_line(data = true_age,
            mapping = aes(y = age, x = position),
            inherit.aes = FALSE,
            color = 2,
            linetype = 'dashed')
plot(model_output, type = 'sed_rate')
plot(model_output, type = 'periodogram')
plot(model_output, type = 'trace')


new_positions <- data.frame(position = c(0, 4, 7.5, 15.0, 22.5, 25, 26.23),
                            thickness = c(0, 0,0,0,0,0,0),
                            id = letters[1:7])



case_1 <- tribble(~position, ~relative_age,
                  0,           6.000,
                  4,           6.296,
                  7.5,         6.562,
                  15.0,        7.124,
                  22.5,        7.733,
                  25.0,        7.904,
                  26.23,       8.000)

predictions <- predict(model_output, new_positions = new_positions)




(predictions$CI$median - case_1$relative_age)




