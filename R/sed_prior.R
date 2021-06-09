sed_prior <- function(rate, min, max){
  LL <- dunif(rate, min, max) %>% log() %>% sum()
  return(LL)
}