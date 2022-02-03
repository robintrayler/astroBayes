anchor_hiatus_sed_model <- function(segment_edges,
                                    sed_rate,
                                    anchor_point,
                                    hiatus_duration = NA,
                                    position_grid) {

  # INPUTS:
  # segment_edges: a data frame of segment boundaries
  # sed_rate: a vector of n sedimentation rates that correspond to each segment
  # anchor_point: a vector of length 1 that anchors the top of the model
  # hiatus_duration: a vector of length n, that contains the duration of n
  # hiatuses
  # position_grid: evenly spaced grid of stratigraphic positions that span the
  # range of segment_edges$position

  # OUTPUTS:
  # anchored_model: a vector of evenly spaced ages that correspond to the
  # stratigraphic positions in position_grid

  # anchor the model ----------------------------------------------------------
  anchored_model <- c(0, cumsum(diff(segment_edges$position) / sed_rate)) +
    anchor_point

  # interpolate onto an even grid ---------------------------------------------
  f <- approxfun(x = segment_edges$position,
                 y = anchored_model)

  # store the results in a data frame
  grid <- data.frame(position = position_grid) %>%
    mutate(age = f(position))


  # find the position of any hiatuses -----------------------------------------
  if(any(segment_edges$hiatus_boundary)) {
    hiatus_position <- which(segment_edges$hiatus_boundary)

    # adjust the age at the hiatus points -------------------------------------
    for(i in seq_along(hiatus_position)) {
      grid$age[grid$position > segment_edges$position[hiatus_position[i]]] <-
        grid$age[grid$position > segment_edges$position[hiatus_position[i]]] +
        hiatus_duration[i]
    }
  }
  # return an anchored model --------------------------------------------------
  return(grid$age)

}
