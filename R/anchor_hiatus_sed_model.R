anchor_hiatus_sed_model <- function(segment_edges,
                                    sed_rate,
                                    anchor_point,
                                    hiatus_duration,
                                    position_grid) {

  # find the position of any hiatuses
  hiatus_position <- which(segment_edges$hiatus_boundary)

  # anchor the model
  anchored_model <- c(0, cumsum(diff(segment_edges$position) / sed_rate)) +
    anchor_point

  # interpolate onto an even grid
  f <- approxfun(x = segment_edges$position,
                 y = anchored_model)

  grid <- data.frame(position = position_grid) %>%
    mutate(age = f(position))

  # adjust the age at the hiatus points


  for(i in seq_along(hiatus_position)) {
    grid$age[grid$position > segment_edges$position[hiatus_position[i]]] <-
      grid$age[grid$position > segment_edges$position[hiatus_position[i]]] +
      hiatus_duration[i]
  }

  return(grid$age)

}
