#' Compute log likelihood of data following MBEV distribution
#'
#' @param Y (nx2)-matrix sample
#' @param pars parameter vector in the following order (mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep)
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
  condition = is.matrix(Y) && ncol(Y) == 2 && n > 1 &&
    delta1 > -1 &&  delta2 > -1 && 
    sigma1 > 0  &&  sigma2 > 0  &&
    dep > 0 && dep <= 1  
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
  
  bdgdelta <- evd::dbvevd(
    x = tm,
    dep = dep,
    model = "log",
    mar1 = c(0, sigma1, xi1),
    mar2 = c(0, sigma2, xi2)
  ) * 
  ((delta1 + 1) * (delta2 + 1)) * (abs(X1 - mu1) ** delta1) * (abs(X2 - mu2) ** delta2)
  # Log:
  log_likelihood <- sum(log(bdgdelta), na.rm = TRUE)
  
  # Return
  return(log_likelihood)
}


#' Estimate parameters of MBEV distribution from data
#'
#' @param Y (nx2)-matrix sample
#' @param par_start starting value for optimization as a vector in the following order:
#' \code{mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep}
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
  
  # estimate marginal distribution to get good starting values
  Y1_estimators <- bgev::bgev.mle(x = as.vector(Y1), deoptim.itermax = 20)$par # returns (mu, sigma, xi, delta)
  Y2_estimators <- bgev::bgev.mle(x = as.vector(Y2), deoptim.itermax = 20)$par
  mu1 = Y1_estimators[1]; sigma1 = Y1_estimators[2]; xi1 = Y1_estimators[3]; delta1 = Y1_estimators[4]
  mu2 = Y2_estimators[1]; sigma2 = Y2_estimators[2]; xi2 = Y2_estimators[3]; delta2 = Y2_estimators[4]
  
  
  # use DEOPTIM to get good starting values for dep parameter in (0,1]
  dep_start = DEoptim::DEoptim(fn = function(dep) -mbev_log_likelihood(Y, c(mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep)), 
                   lower = 0.001, upper = 1, control = DEoptim::DEoptim.control(itermax = 5, 
                                                                      trace = FALSE))$optim$bestmem
  
  # get initial, lower and upper bound for optimization
  # format (mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep)
  par_start = c(mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep_start)
  par_lower = c(mu1 - 5*sigma1, mu2 - 5*sigma2, -0.999, -0.999, 0.001, 0.001, xi1-5*abs(xi1), xi2-5*abs(xi2), 0.001)
  par_upper = c(mu1 + 5*sigma1, mu2 + 5*sigma2, delta1*5, delta2*5, 5*sigma1, 5*sigma2, xi1+5*abs(xi1), xi2+5*abs(xi2), 0.999)
  
  
  fn <- function(par) {
    val <- -mbev_log_likelihood(Y, par)
    if (!is.finite(val)) return(1e99)
    val
  }
  
  optim_result <- optim(par = par_start, 
           fn = fn,
           lower=par_lower,
           upper=par_upper,
           method="L-BFGS-B")
  
  # return 
  return(optim_result)
}



