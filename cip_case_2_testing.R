library(astroBayes)
library(tidyverse)
theme_set(theme_minimal())
library(astrochron)
library(viridis)

# load CIP test case 2 --------------------------------------------------------
cyclostrat_data <- read.csv(file = './data/CIP_data/case_2/cyclostrat_data.csv')
true_age <- read.csv(file = './data/CIP_data/case_2/true_age.csv')

# tune on precession only
tuning_frequency <- read.csv(file = './data/data_A/tuning_frequency.csv') %>%
  filter(orbital_cycle == 'precession')

# periodogram -----------------------------------------------------------------
cyclostrat_data %>%
  periodogram(output = 1,
              genplot = FALSE,
              background = 1) %>%
  ggplot(mapping = aes(x = Frequency,
                       y = Power)) +
  geom_line() +
  xlim(0, 7.5) +
  geom_line(mapping = aes(y = AR1))

# filter out high frequency noise ---------------------------------------------
cyclostrat_data_filtered <- cyclostrat_data %>%
  bandpass(flow = 2.5,
           fhigh = 7.5,
           genplot = FALSE)

# eha analysis ----------------------------------------------------------------
eha_results <- astrochron::eha(cyclostrat_data_filtered,
                               win = 2,
                               step = 0.01,
                               genplot = FALSE,
                               output = 3,
                               verbose = FALSE,
                               ydir = -1,
                               fmax = 10)

eha_results <- eha_results %>% pivot_longer(cols = 2:(ncol(eha_results)),
                                            names_to = 'depth',
                                            values_to = 'amplitude')

eha_results$depth <- sub("X", "", eha_results$depth) %>% as.numeric()

# plot the results
eha_results %>%
  ggplot(mapping = aes(x = freq,
                       y = depth,
                       fill = amplitude)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis(option = 'mako') +
  xlab('Frequency (cycles/m)') +
  ylab('Depth (m)') +
  scale_y_reverse() +
  ggtitle('Evolutive Harmonic Analysis') +
  theme(legend.position = 'none')  +
  xlim(0, 7.5) +
  geom_hline(yintercept = 5 + 0.75,
             color = 'white')

# create synthetic geochronology ----------------------------------------------
n = 3 # number of points to generate
# pick some true ages
date_positions <- true_age[sample(seq_along(true_age[, 1]), n), ]
geochron_data  <- # assemble into synthetic geochronology
  data.frame(age = rnorm(n, mean = date_positions$age, sd = 0.01),
             age_sd = rnorm(n, date_positions$age * 0.01, sd = 0.001),
             position = date_positions$position,
             thickness = 0,
             id = letters[1:n]) %>%
  arrange(position)

# plot the geochronology ------------------------------------------------------
geochron_data %>%
  ggplot(mapping = aes(x = age,
                       y = position,
                       color = id)) +
  geom_point() +
  geom_linerange(mapping = aes(xmin = age - age_sd*2,
                               xmax = age + age_sd*2)) +
  scale_y_reverse()

# develop segment edges -------------------------------------------------------
segment_edges <-
  data.frame(position = c(0, 3, 5.75, 10),
             thickness = c(0, 0, 0.1, 0),
             hiatus_boundary = c(FALSE, FALSE, TRUE, FALSE))

# run the model ---------------------------------------------------------------
model_output <- astro_bayes_model(geochron_data = geochron_data,
                                  cyclostrat_data = cyclostrat_data_filtered,
                                  tuning_frequency = tuning_frequency,
                                  segment_edges = segment_edges,
                                  iterations = 50000,
                                  burn = 2000)

# plot the results ------------------------------------------------------------
pdf(file = 'cip_case_2_test_3_25_2022.pdf', width = 12, height = 12)
plot(model_output, type = 'age_depth') +
  geom_line(data = true_age,
            mapping = aes(x = position,
                          y = age),
            inherit.aes = FALSE,
            color = 2,
            linetype = 'dashed')
plot(model_output, type = 'trace')
plot(model_output, type = 'sed_rate')
plot(model_output, type = 'periodogram')
plot(model_output, type = 'cyclostrat')
dev.off()
