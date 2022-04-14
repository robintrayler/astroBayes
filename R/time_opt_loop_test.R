library(astroBayes)
library(astrochron)
library(tidyverse)

theme_set(theme_minimal())
source('./R/time_opt_likelihood.R')
source('./R/adaptive_update.R')
source('./R/pgram_likelihood.R')
# load the example data sets --------------------------------------------------
cyclostrat_data  <- read.csv(file = './data/data_A/cyclostratigraphic_record.csv')
tuning_frequency <- read.csv(file = './data/data_A/tuning_frequency.csv') %>%
  filter(orbital_cycle != 'obliquity')
segment_edges    <- read.csv(file = './data/data_A/segment_edges.csv')

# set up model parameters -----------------------------------------------------
iterations = 10000
burn  = 1000
# convert to work with default timeOpt frequencies
tuning_target <- 1/tuning_frequency$frequency * 1000

# set up storage and initial values
sed_rate <- matrix(nrow = iterations,
                   ncol = nrow(segment_edges) - 1)

sed_rate[1, ] <- rnorm(nrow(segment_edges) - 1, mean = 1, sd = 0.1)

# set up progress bar
pb <- progress::progress_bar$new(total = iterations,
                                 format = '[:bar] :percent eta: :eta')

for(i in 2:iterations) {
  pb$tick()
  sed_rate[i, ] <- sed_rate[i - 1, ]

  for(j in seq_along(sed_rate[i - 1, ])) {
    # propose new rate
    new_rate <- adaptive_update(chain = sed_rate[, j],
                                i = i,
                                lower = 0,
                                upper = 100,
                                initial_Cd = 0.001,
                                start_index = 100)

    # calculate probability of new rate
    current_cyclostrat <- cyclostrat_data %>%
      filter(position > segment_edges$position[j] &
               position < segment_edges$position[j + 1])

    LL_proposed <- time_opt_likelihood(current_cyclostrat,
                                       tuning_frequency = tuning_target,
                                       sed_rate = new_rate)$likelihood_combined

    LL_current <- time_opt_likelihood(current_cyclostrat,
                                      tuning_frequency = tuning_target,
                                      sed_rate = sed_rate[i - 1, j])$likelihood_combined

    p <- (LL_proposed) - (LL_current)

    if(!is.na(p)) {
      if(!is.infinite(p)) {
        if(exp(p) > runif(1)) {
          sed_rate[i, j] <- new_rate
        }
      }
    }
  }
}

sed_rate <- sed_rate %>%
  as.data.frame() %>%
  add_column(iteration = 1:iterations)
sed_rate <- sed_rate %>%
  pivot_longer(cols = colnames(.)[-8])

sed_rate %>%
  filter(iteration > burn) %>%
  ggplot(mapping = aes(x = iteration,
                       y = value * 10,
                       color = name)) +
  geom_line() +
  facet_wrap(name~.,
             scales = 'free')

sed_rate %>%
  filter(iteration > burn) %>%
  ggplot(mapping = aes(x = value * 10,
                       fill = name)) +
  geom_density() +
  facet_wrap(name~.,
             scales = 'free')


sed_rate %>%
  filter(iteration > burn) %>%
  group_by(name) %>%
  summarise(mean_rate = mean(value)*10,
            sd_rate = sd(value)*10)
