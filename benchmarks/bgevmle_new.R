bgev.mle.new <- function (x, method_envstats = "mle", deoptim.itermax = 200, 
          optim.method = "L-BFGS-B", start = NULL, lower = NULL, upper = NULL, verbose = FALSE) 
{
  likbgev2 = function(theta) {
    val = -likbgev(x, theta)
    if (!is.finite(val)) 
      return(1e+100)
    else return(val)
  }
  if (is.null(start)) {
    fit_gev <- try(EnvStats::egevd(x = x, method = method_envstats))
    if (inherits(fit_gev, "try-error")) {
      if(verbose)
        message("I tried to get good starting values using egevd and it failed. Try using a different method for method_envstats, e.g, pwme. See the help of egevd for the available options.")
      fit_gev <- try(EnvStats::egevd(x = x, method = "pwme"))
      #return(NULL)
    }
    mu_start = fit_gev$parameters[1]
    sd_start = fit_gev$parameters[2]
    xi_start = -fit_gev$parameters[3]
    starts = c(mu_start, sd_start, xi_start, 0.1)
    xi_min = xi_start/5
    xi_max = xi_start * 5
    if (xi_start < 0) {
      xi_max = xi_start/5
      xi_min = xi_start * 5
    }
  }
  if (is.null(lower) | is.null(upper)) {
    lower = c(mu_start - 2 * sd_start, sd_start/5, xi_min, 
              -0.99)
    upper = c(mu_start + 2 * sd_start, 5 * sd_start, xi_max, 
              20)
  }
  starts.DEoptim = DEoptim::DEoptim(fn = likbgev2, lower, upper, 
                                    control = DEoptim::DEoptim.control(itermax = deoptim.itermax, 
                                                                       trace = FALSE))
  esti <- stats::optim(par = starts.DEoptim$optim$bestmem, 
                       fn = likbgev2, method = optim.method, lower = lower, 
                       upper = upper)
  esti
}