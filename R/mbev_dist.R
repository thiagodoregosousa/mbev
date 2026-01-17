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
#' Uses gumbel bivariada from evd::rbvevd
#' Parameter support:
#' \itemize{
#'   \item \code{n > 1}
#'   \item \code{delta1 > -1}, \code{delta2 > -1}
#'   \item \code{sigma1 > 0}, \code{sigma2 > 0}
#'   \item \code{modelo = "log"} (logistic symmetric EV dependence)
#'   \item \code{dep >= 1} (dependence parameter)
#' }
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
                   xi2 = 0.5,
                   dep = 1,
                   modelo = "log") {
  
  # check input parameters 
  condition = (n > 1) &&
              delta1 > -1 &&  delta2 > -1 && 
              sigma1 > 0  &&  sigma2 > 0  &&
              modelo == "log" &&
              dep > 0 && dep <= 1  
  if(!condition)
    stop("Invalid input parameters")
  
  # dep is in [1,Inf) but evd::rbvevd uses r in (0,1] so that is why we use 1/dep (see Details of ?dbvevd for the model = "log")
  x <- evd::rbvevd(
    n = n,
    dep = 1/dep,
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



