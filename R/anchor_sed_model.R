anchor_sed_model <- function(segment_edges,
                             sed_rate,
                             anchor_point){
  anchored_model <- c(0, cumsum(diff(segment_edges) / sed_rate)) + anchor_point
  return(anchored_model)
}
