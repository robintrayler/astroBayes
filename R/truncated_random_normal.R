#' Truncated Random Normal
#' @name truncated_random_normal
#' @param mean mean value
#' @param sd standard deviation
#' @param low lower truncation value
#' @param high upper truncation value
#'
#' @return a random number
#'
#' @md
#' @export
##------------------------------------------------------------------------------------------
## Function for truncated random normal with a mean and sd
truncated_random_normal <- function(mean,
                                    sd,
                                    low  = 1e-10,
                                    high = 1e10){
  lowlimold <- (low - mean) / sd
  upplimold <- (high - mean) / sd
  y <- truncated_standard_normal(lowlimold, upplimold)
  newvalue <- mean + sd * y
  return(newvalue)
}
