gen_cycles <- function(sed_rate, target_in, n_pts, dx)
{
  result <- matrix(0, n_pts, 2 * length(target_in))
  for (i in 1:length(target_in))
  {
    result[,2 * i - 1] <- cos( (2 * pi) / (target_in[i]) * (dx / sed_rate) * (1:n_pts))
    result[,2 * i]     <- sin( (2 * pi) / (target_in[i]) * (dx / sed_rate) * (1:n_pts))
  }
  return(result)
}

# function to perform fitting
#  dx, n_pts passed into function transparently
fit_it <- function(sed_rate,time_series,target_in, n_pts, dx)
{
  xm <- gen_cycles(sed_rate, target_in, n_pts, dx)

    fit_0 <- lm(time_series[, 2] ~ xm)

  # calculate rho and sigma based on residuals
  rho <- cor(fit_0$residuals[1:(n_pts-1)],fit_0$residuals[2:n_pts])
  sigma <- sd(fit_0$residuals)

  logLL <- calc_LL(fit_0$residuals, rho, sigma, n_pts)

  return(data.frame(sed_rate,
               logLL,
               rho,
               sigma))
}

# function to perform log-likelihood calculation, including
#  assessment of correlated residuals, assuming AR1 model
calc_LL <- function(residuals, rho, sigma, n_pts)
{
  # calculate Re^-1 (ReInv), as in EQ. A-7 of Malinverno & Briggs (2004)
  # set up array, with zeros
  ReInv <- double(n_pts * n_pts)
  dim(ReInv) <-c(n_pts, n_pts)
  # put 1+rho^2 on diagonal
  ReInv[row(ReInv)==col(ReInv)] = 1 + rho^2
  # except at (1,1) and (n_pts,n_pts), which have a value of 1
  ReInv[1, 1] = 1
  ReInv[n_pts, n_pts] = 1
  # put -rho on subdiagonal
  ReInv[(row(ReInv) - 1) == col(ReInv)] = -1 * rho
  # put -rho on superdiagonal
  ReInv[(row(ReInv) + 1) == col(ReInv)] = -1 * rho
  # now multiple matrix by 1/(1-rho^2)
  ReInv = ReInv/(1-rho^2)
  # calculate log-likelihood
  logLL = ( n_pts * log(2*pi) ) + ( 2*n_pts*log(sigma) ) + ( (n_pts-1)*log(1-rho^2) )
  logLL = logLL + (1/(sigma^2)) * t(residuals) %*% ReInv %*% residuals
  logLL = -0.5 * logLL
  return(logLL)
}
