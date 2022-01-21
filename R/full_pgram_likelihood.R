# INPUTS
# rate = vector of sedimentation rates
# segment_edges = sedimentation rate change points
# cyclostrat = cyclostratigraphic record spanning the range of segment_edges
# tuning frequencies = vector of tuning frequencies to use
full_pgram_likelihood <- function(sed_rate,
                                  segment_edges,
                                  cyclostrat,
                                  tuning_frequency) {
  # prepare the data ----------------------------------------------------------
  age_model <- c(0, cumsum(diff(segment_edges) / sed_rate))

  f_1 <- approxfun(x = segment_edges,
                   y = age_model)

  f_2 <- cyclostrat |>
    mutate(position = f_1(position)) |>
    astrochron::linterp(genplot = FALSE,
            verbose = FALSE) |>
    astrochron::periodogram(output = 1,
                verbose = FALSE,
                genplot = FALSE,
                background = 1) |>
    mutate(probability = Power / AR1_Fit,
           probability = probability / sum(probability),
           time_freq   = Frequency) |>
    with(approxfun(x = time_freq,
                   y = probability))

  LL <- f_2(tuning_frequency) |> log() |> sum()
  return(LL)
}
