#' A faster way to baseline hazard function for Cox model
#' 
#' @export
#' @param a coef for log_lambda 
#' @param b coef for Lambda
#' @param t observation time
#' @param cum output cumulative/instantaneous baseline hazard function
#' @example 
#' library(survival)
#' setwd(123)
#' N <- 5000
#' x <- rnorm(N)
#' error <- log( - log( runif(N) ) )
#' t <- exp( - x * 2 + error)  # Cox model with exponential baseline
#' cen <- rexp( N )
#' time <- pmin(t, cen)
#' status <- t < cen
#' fit <- coxph(Surv(time, status) ~ x)
#' fit
#'
#' ## Derived Baseline Hazard
#' baz0 = basehaz_ab( a = status, b = exp(fit$coef * x), t = time)
#' ## Buildin Baseline Hazard
#' baz1 = basehaz(fit, centered = F)
#' 
#' ## Check Equivalence
#' plot( time, baz0)
#' lines(baz1$time, baz1$hazard, col = 2)
#' 
#' max( abs(sort(baz0) - baz1$hazard) )

basehaz_ab <- function(a, b, t, cum = TRUE){
  row = 1:length(a)
  t_order <- order(t)
  a_sort = a[t_order]
  b_sort = b[t_order]
  
  denom = sum(b_sort) -  cumsum(b_sort) + b_sort # at risk cumulative sum
  lambda = a_sort / ( denom )
  if(cum == TRUE) lambda =  cumsum(lambda)
  lambda = lambda[ order(t_order) ]
  lambda
}






