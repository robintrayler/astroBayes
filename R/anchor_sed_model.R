anchor_sed_model <- function(layer_boundaries,
                                    sed_rate,
                                    anchor_point,
                                    hiatus_duration = NA,
                                    position_grid) {

  # INPUTS:
  # layer_boundaries: a data frame of layer boundaries
  # sed_rate: a vector of n sedimentation rates that correspond to each segment
  # anchor_point: a vector of length 1 that anchors the top of the model
  # hiatus_duration: a vector of length n, that contains the duration of n
  # hiatuses
  # position_grid: evenly spaced grid of stratigraphic positions that span the
  # range of layer_boundaries$position

  # OUTPUTS:
  # anchored_model: a vector of evenly spaced ages that correspond to the
  # stratigraphic positions in position_grid

  # anchor the model ----------------------------------------------------------
  anchored_model <- c(0, cumsum(diff(layer_boundaries$position) / sed_rate)) +
    as.vector(anchor_point)

  # interpolate onto an even grid ---------------------------------------------
  f <- approxfun(x = layer_boundaries$position,
                 y = anchored_model)

  # store the results in a data frame
  grid <- data.frame(position = position_grid) %>%
    mutate(age = f(position))


  # find the position of any hiatuses -----------------------------------------
  if(any(layer_boundaries$hiatus_boundary)) {
    hiatus_position <- which(layer_boundaries$hiatus_boundary)

    # adjust the age at the hiatus points -------------------------------------
    for(i in seq_along(hiatus_position)) {
      grid$age[grid$position > layer_boundaries$position[hiatus_position[i]]] <-
        grid$age[grid$position > layer_boundaries$position[hiatus_position[i]]] +
        hiatus_duration[i]
    }
  }
  # return an anchored model --------------------------------------------------
  return(grid$age)

}
