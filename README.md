# astroBayes

`astroBayes` is a Bayesian framework for combining radioisotopic dates and cyclostratigraphic data into age-depth models. The background of the modeling framework and details of model construction are covered in: 

**Trayler, R. B., Meyers, S. R., Sageman, B. B., Schmitz, M. D., (in prep) *Bayesian Integration of Astrochronology and Radioisotope Geochronology***. 

A draft version of this manuscript is available at [github.com/robintrayler/astroBayes_manuscript](https://github.com/robintrayler/astroBayes_manuscript). 

## Installation

`astroBayes` is as an `R` package. It can be installed using the [`devtools`](https://github.com/r-lib/devtools) R package. 

```r
devtools::install_github('robintrayler/astroBayes')
```

Once `astroBayes` is installed it can be loaded in the usual way for `R` packaged by adding `library(astroBayes)` to the beginning of an `R` script.

## Using `astroBayes`

`astroBayes` is designed to fit age-depth models to radioisotope geochronology and cyclostratigraphic data. `astroBayes` is not designed to test for the presence of astronomical signals in data, instead it is intended to be used in conjunction with the testing methods available (such as those available in [`{astrochron}`](https://geoscience.wisc.edu/research/x-ray-fluorescence-xrf-scanner-lab/astrochron-a-computational-tool-for-astrochronology/)). 

First load astroBayes and any other packages you want to use. 

```
library(astroBayes)        # for modeling
```

`astroBayes` includes a simple data set consisting of a set of radioisotopic dates, a cyclostratigraphic record, a set of target frequencies, and a set of layer boundaries. The example data can be loaded using the `data()` function. 

```
data("target_frequencies")
data("radioisotopic_dates")
data("cyclostratigraphic_data")
data("layer_boundaries")
```

The primary function in `astroBayes` is `astro_bayes_model()`. This function function uses a Metropolis-Hastings Markov Chain Monte Carlo algorithm to fit an age-depth model to the data by finding the most probable sedimentation rates for each model layer. Running the model will take a few minutes to several hours depending on the number of `iterations` and the complexity of the data (number of `dates`,  and number of layers). A progress bar should pop up in the `R` terminal with a rough estimate of time remaining.

```
age_model <- astro_bayes_model(geochron_data = dates,
                               cyclostrat_data = cyclostrat,
                               tuning_frequency = target_frequencies,
                               segment_edges = layer_boundaries,
                               iterations = 10000,
                               burn = 1000)
```

After the model has finished running you can visualize the results. 

The `age_depth` plot shows the age-depth model as a median (black line) and 95% credible interval (shaded grey region). The dates are shown as colored normal distribution. 

```
plot(age_model, type = 'age_depth')
```
![](./figures/age_depth.jpeg)