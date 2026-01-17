#' Compute log likelihood of data following MBEV distribution
#'
#' @param Y (nx2)-matrix sample
#' @param pars parameter vector in the following order (mu1, mu2, delta1, delta2, sigma1, sigma2, x1, x2, dep)
#'
#' @return blogl numeric loglikelihood of MBEV model evaluated at pars
#' 
#' @seealso \link{rmbev}
#' 
#' @references
#' YASMIN_TESE
#' 
#' MBEV_ARTIGO
#' 
#' @export
mbev_log_likelihood <- function(Y, pars = c(0,0,0,0,1,1,0.5,0.5,1)) {

  # get parameters
  mu1    <- pars[1]
  mu2    <- pars[2]
  delta1 <- pars[3]
  delta2 <- pars[4]
  sigma1 <- pars[5]
  sigma2 <- pars[6]
  xi1    <- pars[7]
  xi2    <- pars[8]
  dep    <- pars[9]
  n <- dim(Y)[1]
  
  # check input parameters 
  condition = all.equal(dim(Y),c(n,2)) &&
    delta1 > -1 &&  delta2 > -1 && 
    sigma1 > 0  &&  sigma2 > 0  &&
    dep >= 1 
  if(!condition)
    return(1e99)
  
  # set minimum value for alpha
  #alfa <- ifelse(dep > 0, dep, 0.00001)
  
  # Compute auxiliary variables:
  X1 = Y[, 1]
  X2 = Y[, 2]
  
  t1 <- (X1 - mu1) * ((abs(X1 - mu1)) ** delta1)
  t2 <- (X2 - mu2) * ((abs(X2 - mu2)) ** delta2)
  
  # Compute density points
  tm = matrix(c(t1, t2), ncol = 2)
  
  bdgdelta <- dbvevd(
    x = tm,
    dep = 1/alpha,
    model = "log",
    mar1 = c(0, sigma1, xi1),
    mar2 = c(0, sigma2, xi2)
  ) * 
  ((delta1 + 1) * (delta2 + 1)) * (abs(X1 - mu1) ** delta1) * (abs(X2 - mu2) ** delta2)
  # Log:
  blogl <- sum(log(bdgdelta), na.rm = TRUE)
  
  # Return negative Value for maximization
  return(-blogl)
}


#' Estimate parameters of MBEV distribution from data
#'
#' @param Y (nx2)-matrix sample
#' @param par_start starting value for optimization as a vector in the following order:
#' \code{mu1, mu2, delta1, delta2, sigma1, sigma2, x1, x2, dep}
#' @param par_lower lower bounds in the same format as \code{par_start}
#' @param par_upper upper bounds in the same format as \code{par_start}
#' @return optim_result result of the optim function (see \link{optim}).
#' Use \code{optim_result$par} to get par the estimated parameters
#' 
#' @seealso \link{rmbev}
#' 
#' @references
#' YASMIN_TESE
#' 
#' MBEV_ARTIGO
#' 
#' @export
mbev_estimation = function(Y, 
                           par_start,
                           par_lower,
                           par_upper) {
  
  # get marginal vectors
  Y1=Y[,1]
  Y2=Y[,2]
  
  # estimate marginal distribution
  c(mu1, sigma1, xi1, delta1) <- bgev.mle(Y1, itermax = 20)$par # returns (mu, sigma, xi, delta)
  c(mu2, sigma2, xi2, delta2) <- bgev.mle(Y2, itermax = 20)$par
  
  
  # get initial, lower and upper bound for optimization
  # format (mu1, mu2, delta1, delta2, sigma1, sigma2, x1, x2, dep)
  par_start = c(mu1, mu2, delta1, delta2, sigma1, sigma2, x1, x2, 2)
  par_lower = c(mu1 - 5*sigma1, mu2 - 5*sigma2, -1 + 0.001, -1 + 0.001, 0.001, 0.001, x1-5*abs(x1), x1-5*abs(x2), 1)
  par_upper = c(mu1 + 5*sigma1, mu2 + 5*sigma2, max(delta1,10), max(delta2,10), 5*sigma1, 5*sigma2, x1+5*abs(x1), x1+5*abs(x2), 10)
  
  optim_result <- optim(par = thetaalfa, 
           fn=function(par) mbev_log_likelihood(Y, par), 
           Y=Y,
           lower=c(0.001),
           upper=c(1),
           method="L-BFGS-B")$par
  
  # return 
  return(optim_result)
}



