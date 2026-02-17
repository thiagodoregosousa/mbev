source("R/mbev_dist.R")
source("R/mbev_estimation.R")


# create design
library(tibble)
library(dplyr)
theta_tbl <- tibble(
  theta_id = paste0("Theta", 1:6),
  mu1    = c(0,   0,   0,   0,   0,   0),
  mu2    = c(0,   0,   0,   0,   0,   0),
  delta1 = c(0,   0,   0,   1,   1,   0.8),
  delta2 = c(0,   1,   1,   2,   1,   2),
  sigma1 = c(1,   1,   1,   1,   1,   1),
  sigma2 = c(1,   1,   1,   1,   1,   1),
  xi1    = c(0.5, 0.5, 1,  -1,  0.5, 0),
  xi2    = c(0.5, 0.5, 0,  -1,  0.5, 0),
  dep  = c(1,   1,   1/2,   1,   1,   1/2.5)
)
n_vals <- c(50, 100, 500)
Design <- tidyr::crossing(
  n = n_vals,
  theta_tbl
) %>%
  select(
    n, mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep
  )




Generate <- function(condition, fixed_objects) {
  dat <- with(condition, rmbev(n = n, mu1 = mu1, mu2 = mu2, delta1 = delta1, delta2 = delta2, sigma1 = sigma1, sigma2 = sigma2, xi1 = xi1, xi2 = xi2, dep = dep)  )
  dat
}


Analyse <- function(condition, dat, fixed_objects) {
  
  par_names = c("mu1", "mu2", "delta1", "delta2", "sigma1", "sigma2", "xi1", "xi2", "dep") 
  ret_error <- rep(NA_real_, 9)
  names(ret_error) <- par_names
  est <- tryCatch(mbev_estimation(dat, 
                                  control_DEoptim_step1 = DEoptim::DEoptim.control(itermax = 100, NP = 100, trace = FALSE),
                                  control_DEoptim_step2 = DEoptim::DEoptim.control(itermax = 50, NP = 100, trace = FALSE),
                                  control_DEoptim_step3 = DEoptim::DEoptim.control(itermax = 500, NP = 100, trace = FALSE)), error = function(e) NULL)
  if(is.null(est))
    return(ret_error)
  if (is.null(est$optim$bestmem) ||
      length(est$optim$bestmem) != 9 ||
      any(!is.finite(est$optim$bestmem))) {
    return(ret_error)
  }
  ret_success <- as.numeric(est$optim$bestmem)
  names(ret_success) <- par_names
  return(ret_success)
}


Summarise <- function(condition, results, fixed_objects) {
  
  true_mu1 <- condition$mu1
  true_mu2 <- condition$mu2
  true_delta1 <- condition$delta1
  true_delta2 <- condition$delta2
  true_sigma1 <- condition$sigma1
  true_sigma2 <- condition$sigma2
  true_xi1 <- condition$xi1
  true_xi2 <- condition$xi2
  true_dep <- condition$dep
  
  ret <- c(
    
    # Bias
    bias_mu1 = mean(results[, "mu1"]) - true_mu1,
    bias_mu2 = mean(results[, "mu2"]) - true_mu2,
    bias_delta1 = mean(results[, "delta1"]) - true_delta1,
    bias_delta2 = mean(results[, "delta2"]) - true_delta2,
    bias_sigma1 = mean(results[, "sigma1"]) - true_sigma1,
    bias_sigma2 = mean(results[, "sigma2"]) - true_sigma2,
    bias_xi1 = mean(results[, "xi1"]) - true_xi1,
    bias_xi2 = mean(results[, "xi2"]) - true_xi2,
    bias_dep = mean(results[, "dep"]) - true_dep,
    
    # Empirical SE (what referee wants)
    SE_mu1 = sd(results[, "mu1"]),
    SE_mu2 = sd(results[, "mu2"]),
    SE_delta1 = sd(results[, "delta1"]),
    SE_delta2 = sd(results[, "delta2"]),
    SE_sigma1 = sd(results[, "sigma1"]),
    SE_sigma2 = sd(results[, "sigma2"]),
    SE_xi1 = sd(results[, "xi1"]),
    SE_xi2 = sd(results[, "xi2"]),
    SE_dep = sd(results[, "dep"]),
    
    # RMSE
    RMSE_mu1 = sqrt(mean((results[, "mu1"] - true_mu1)^2)),
    RMSE_mu2 = sqrt(mean((results[, "mu2"] - true_mu2)^2)),
    RMSE_delta1 = sqrt(mean((results[, "delta1"] - true_delta1)^2)),
    RMSE_delta2 = sqrt(mean((results[, "delta2"] - true_delta2)^2)),
    RMSE_sigma1 = sqrt(mean((results[, "sigma1"] - true_sigma1)^2)),
    RMSE_sigma2 = sqrt(mean((results[, "sigma2"] - true_sigma2)^2)),
    RMSE_xi1 = sqrt(mean((results[, "xi1"] - true_xi1)^2)),
    RMSE_xi2 = sqrt(mean((results[, "xi2"] - true_xi2)^2)),
    RMSE_dep = sqrt(mean((results[, "dep"] - true_dep)^2))
  )
  
  return(ret)
}



run_monte_carlo_mbev <- function(lines_to_use,
                                 person,
                                 ncores = 7) {
  
  # ---- Worker map ----
  worker_map <- c(
    thiago = 1,
    yasmin = 2,
    cira   = 3,
    raul   = 4
  )
  
  if (!person %in% names(worker_map))
    stop("Unknown person name.")
  
  worker_id <- worker_map[[person]]
  
  # ---- Seed structure ----
  seed <- lines_to_use +
    1000  * worker_id
  
  # ---- Run simulation ----
  Final <- SimDesign::runSimulation(
    design        = Design[lines_to_use, ],
    replications  = 500,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    progress      = TRUE,
    verbose       = TRUE,
    store_results = TRUE,
    parallel      = TRUE,
    ncores        = ncores,
    save_results  = TRUE,
    seed          = seed,
    save_details = list(save_results_dirname = paste0(
      "mc_lines_",
      paste(lines_to_use, collapse = "_"),
      "_",
      person))
  )
  
  # ---- Save file ----
  filename <- paste0(
    "benchmarks/mc_lines_",
    paste(lines_to_use, collapse = "_"),
    "_",
    person,
    ".rds"
  )
  
  saveRDS(Final, file = filename)
  
  print(paste("Finished:", filename))
  
  invisible(Final)
}














