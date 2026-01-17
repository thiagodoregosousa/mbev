#' Generate random iid bivariate vectors following MBEV distribution
#'
#' @param n sample size
#' @param (mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep) mu,delta,sigma,xi are interpreted as the parameters 
#' of the marginal BGEV and dep is the dependency parameter passed to evd::rbvevd.
#' @param D vector of means mu1 and mu2 as in pg 44 of YASMIN_TESE
#' @param modelo (only "log" implemented and tested) passed to evd::rbvevd. Only tested with `log` (symmetric logistic)
#'
#' @return xy (nx2)-matrix simulated iid sample of size n following MBEV distribution
#' as a matrix 
#'
#' @details
#' gumbel bivariada from evd::rbvevd
#' 
#' @references 
#' YASMIN_TESE
#' 
#' MBEV_ARTIGO
#' 
#'
#' @export
rmbev <- function (n = 50,
                   mu1 = 0,
                   mu2 = 0,
                   delta1 = 0,
                   delta2 = 0,
                   sigma1 = 1,
                   sigma2 = 1,
                   xi1 = 0.5,
                   xi2 = .5,
                   modelo = "log",
                   dep = 1) {
  
  # check input parameters 
  condition = (n > 1) &&
              delta1 > -1 &&  delta2 > -1 && 
              sigma1 > 0  &&  sigma2 > 0  &&
              modelo == "log" &&
              dep > 0 && dep <= 1  
  if(!condition)
    stop("Invalid input parameters")
  
  
  x <- evd::rbvevd(
    n = n,
    dep = dep,
    model = modelo,
    mar1 = c(0, sigma1, xi1),
    mar2 = c(0, sigma2, xi2)
  )
  y1 <- sign(x[, 1]) * abs(x[, 1])^(1 / (delta1 + 1)) + mu1
  y2 <- sign(x[, 2]) * abs(x[, 2])^(1 / (delta2 + 1)) + mu2
  
  
  y <- matrix(
    c(y1, y2),
    ncol = 2,
    nrow = n,
    dimnames = list(c(), c("y1", "y2"))
  )
  return(y)
}



