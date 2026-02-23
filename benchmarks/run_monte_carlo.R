rm(list = ls(all.names = TRUE))
source("benchmarks/monte_carlo_mbev_estimation.R")
# dones: lines 1-14 # 1 at the end 
run_monte_carlo_mbev(lines_to_use = 1, person = "thiago", ncores = 8)
Sys.sleep(10)
run_monte_carlo_mbev(lines_to_use = 17, person = "thiago", ncores = 7)
Sys.sleep(10)
run_monte_carlo_mbev(lines_to_use = 18, person = "thiago", ncores = 7)
Sys.sleep(10)

