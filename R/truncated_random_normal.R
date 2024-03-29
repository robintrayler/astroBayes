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


## Function for truncated random normal
truncated_standard_normal <- function(a, b){
  # this function draws a random number from a normal distribution centered at
  # 0 with a mean of 1, with truncations at a and b.
  # INPUTS: a = lower truncation
  #         b = upper truncation
  # OUTPUT: x = random number
  accept = FALSE
  A <- atan(a)
  B <- atan(b)
  maxA <- exp((-(a^2) / 4)) / cos(A)
  maxB <- exp((-(b^2) / 4)) / cos(B)
  maxR <- max(maxA, maxB)
  if((a < 1)& (b > 1)){
    maxR <- exp(-.25) * sqrt(2)
  }

  while (accept == FALSE) {
    r2 <- runif(1, 0, 1)
    r <- sqrt(r2) * maxR
    th <- runif(1, A, B)
    u <- r * cos(th)
    x <- tan(th)
    accept <- (x^2) < (log(u) * -4)
  }
  return(x)
}
