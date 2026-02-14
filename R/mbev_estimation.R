# XXXX TMP 
# TEMPORARY SOURCE FILES FOR MODIFIED bgev_mle
source("../bgev_github/R/bgev_domain.R")
source("../bgev_github/R/bgev_distribution.R")
source("../bgev_github/R/bgev_estimation.R")



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





#' Check if parameters of mbev distribution are valid
#' 
#' @param pars \code{mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep}
#' @return boolean TRUE if pars are valid or FALSE otherwise
mbev_valid_pars = function(pars){
  delta1 = pars[3]
  delta2 = pars[4]
  sigma1 = pars[5]
  sigma2 = pars[6]
  dep = pars[9]
  
  if(delta1 > -1 &&  delta2 > -1 && 
     sigma1 > 0  &&  sigma2 > 0  &&
     dep > 0 && dep <= 1  )
    return(TRUE)
  else
    return(FALSE)
}




#' Maximum Likelihood Estimation for the BGEV distribution
#' 
#' Estimate parameters of MBEV distribution from data
#'
#' @param Y (nx2)-matrix sample
#' @param par_start starting value for optimization as a vector in the following order:
#' \code{mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep}
#' @param control_DEoptim_step1 List of type DEoptim::DEoptim.control (PUT LINK HERE) for step1 estimation
#' @param DEoptim_replicates_step1 Number of DEoptim runs to get best likelihood
#' @param lower Optional vector of lower bounds for the parameters \code{mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep}
#' @param upper Optional vector of upper bounds as in \code{lower}
#' @return best An object of class DEOptim with the final estimation
#' 
#' @author Thiago do Rego Sousa and Yasmin Lirio
#' 
#' @seealso \link{rmbev}
#' 
#' @references
#' YASMIN_TESE
#' 
#' MBEV_ARTIGO
#' 
#' @export
mbev_estimation <- function (Y, lower = c(-5,-5, -0.99, -0.99, 0.01, 0.01, -5, -5, 0.01), 
                                 upper = c( 5, 5, 5,    5,    5,   5,   5,   5, 0.99),
                                 control_DEoptim_step1 = DEoptim::DEoptim.control(itermax = 100, NP = 100, trace = FALSE),
                                 control_DEoptim_step2 = DEoptim::DEoptim.control(itermax = 50, NP = 100, trace = FALSE),
                                 control_DEoptim_step3 = DEoptim::DEoptim.control(itermax = 500, NP = 100, trace = FALSE),
                                 DEoptim_replicates_step1 = 5,
                                 DEoptim_replicates_step3 = 10,
                                 trace = TRUE)
{
  # validate input parameters
  if (is.null(as.vector(Y)) || anyNA(Y)) {
    stop("Y must be a non-null matrix vector with no missing values.", call. = FALSE)
  }

  # step1: estimate each one-dimensional vector as bgev using specified lower and upper bounds
  # parameter order: (mu, sigma, xi, delta)
  Y1_est <- bgev_mle(x = as.vector(Y[,1]), lower = lower[c(1,5,7,3)], upper = upper[c(1,5,7,3)], 
                                         control = control_DEoptim_step1, 
                                         DEoptim_replicates = DEoptim_replicates_step1)$optim$bestmem
  Y2_est <- bgev_mle(x = as.vector(Y[,2]), lower = lower[c(2,6,8,4)], upper = upper[c(2,6,8,4)], 
                            control = control_DEoptim_step1, 
                            DEoptim_replicates = DEoptim_replicates_step1)$optim$bestmem 
  if(trace){
    print(c("Y1_est", as.character(Y1_est)))
    print(c("Y2_est", as.character(Y2_est)))
  }
  

  # step2: fix estimated marginal parameters and estimate dependency parameter
  # parameter order   {mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep}
  dep_start = DEoptim::DEoptim(fn = function(dep) -mbev_log_likelihood(Y, c(Y1_est[1], Y2_est[1],
                                                                            Y1_est[4], Y2_est[4],
                                                                            Y1_est[2], Y2_est[2],
                                                                            Y1_est[3], Y2_est[3],
                                                                            dep)), 
                               lower = 0.001, upper = 1, control = control_DEoptim_step2)$optim$bestmem
  if(trace)
    print(paste("dep_start", dep_start))
  
  
  # step3: use mbev full log-likelihood to estimate all parameters jointly on a shrinked grid
  mbev_log_likelihood_negative = function(pars) {
    if(!mbev_valid_pars(pars))
      return(1e+100)
    val = -mbev_log_likelihood(Y, pars)
    if (!is.finite(val)) 
      return(1e+100)
    else return(val)
  }
  
  sz = 2 # size of shrinked grid for final tunning of estimated parameters
  lower_step3 =   c(Y1_est[1] - sz,             Y2_est[1] - sz,
                    max(Y1_est[4] - sz, -0.99), max(Y2_est[4] - sz, -0.99),
                    max(Y1_est[2] - sz, 0.01),  max(Y2_est[2] - sz, 0.01),
                    Y1_est[3] - sz,             Y2_est[3] - sz,
                    0.01) 
  
  upper_step3 =   c(Y1_est[1] + sz,             Y2_est[1] + sz,
                    Y1_est[4] + sz,             Y2_est[4] + sz,
                    Y1_est[2] + sz,             Y2_est[2] + sz, 
                    Y1_est[3] + sz,             Y2_est[3] + sz,
                    0.99) 

  # construct joint log likelihood mbev
  fits <- replicate(DEoptim_replicates_step3, {
    DEoptim::DEoptim(fn = mbev_log_likelihood_negative, control = control_DEoptim_step3, 
                     lower = lower_step3, upper = upper_step3)
  }, simplify = FALSE)
  best <- fits[[ which.min(sapply(fits, function(f) f$optim$bestval)) ]]
  
  if(trace)
    print(paste("best", dep_start))
  
  best
}










