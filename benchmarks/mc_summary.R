library(dplyr)

# read mc
mc = readRDS("benchmarks/mc_summary_all.rds")
mc = mc[-c(6,12,18),] # remove scenario of theta_6
# create column called theta_1,2,3,4,5
mc$theta = factor(rep(1:5,3), levels = 1:5,
                             labels = paste0("Theta[", 1:5, "]")  )


mc[c(1:5),c(10,2,6,8,4,3,7,9,5)]

# -----------------------------
# Convert wide RMSE table to long format
# -----------------------------
df_long <- mc %>%
  pivot_longer(
    cols = starts_with("RMSE_"),
    names_to = "parameter",
    values_to = "RMSE"
  ) %>%
  mutate(
    parameter = gsub("RMSE_", "", parameter),
    
    # 🔥 FORCE correct order here
    parameter = factor(
      parameter,
      levels = c(
        "dep",        # rho first,
        "mu1", "sigma1", "xi1", "delta1",
        "mu2", "sigma2", "xi2", "delta2"
      )
    ),
    
    n = factor(n, levels = c(50, 100, 500)),
    
    theta = theta)

# -----------------------------
# Create parsed math labels for parameters
# -----------------------------
param_labels <- c(
  mu1     = "mu[1]",
  mu2     = "mu[2]",
  delta1  = "delta[1]",
  delta2  = "delta[2]",
  sigma1  = "sigma[1]",
  sigma2  = "sigma[2]",
  xi1     = "xi[1]",
  xi2     = "xi[2]",
  dep     = "r"   # change if needed
)

# -----------------------------
# Plot RMSE heatmap
# -----------------------------
ggplot(df_long, aes(x = n, y = parameter, fill = RMSE)) +
  geom_tile(color = "white") +
  facet_wrap(~ theta, ncol = 3, labeller = label_parsed) +
  scale_y_discrete(
    labels = function(x) parse(text = param_labels[x]),
    limits = rev(levels(df_long$parameter))
  ) +
  scale_fill_viridis_c(option = "C") +
  labs(
    x = "Sample Size (n)",
    y = NULL,
    fill = "RMSE"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )





# BIAS and SE share 
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# -----------------------------
# Pivot RMSE
# -----------------------------
rmse_long <- mc %>%
  pivot_longer(
    cols = starts_with("RMSE_"),
    names_to = "parameter",
    values_to = "RMSE"
  ) %>%
  mutate(parameter = gsub("RMSE_", "", parameter))

# -----------------------------
# Pivot Bias  (lowercase bias_)
# -----------------------------
bias_long <- mc %>%
  pivot_longer(
    cols = starts_with("bias_"),
    names_to = "parameter",
    values_to = "Bias"
  ) %>%
  mutate(parameter = gsub("bias_", "", parameter))

# -----------------------------
# Pivot SE  (uppercase SE_)
# -----------------------------
se_long <- mc %>%
  pivot_longer(
    cols = starts_with("SE_"),
    names_to = "parameter",
    values_to = "SE"
  ) %>%
  mutate(parameter = gsub("SE_", "", parameter))

# -----------------------------
# Merge and compute shares
# -----------------------------
df_all <- rmse_long %>%
  left_join(bias_long, by = c("theta", "n", "parameter")) %>%
  left_join(se_long,   by = c("theta", "n", "parameter")) %>%
  mutate(
    BiasShare = (Bias^2) / (RMSE^2),
    SEShare   = (SE^2)   / (RMSE^2)
  )

df_all <- df_all %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c(
        "dep",
        "mu1", "mu2",
        "delta1", "delta2",
        "sigma1", "sigma2",
        "xi1", "xi2"
      )
    ), 
    # 🔥 FORCE correct order here
    parameter = factor(
      parameter,
      levels = c(
        "dep",        # rho first,
        "mu1", "sigma1", "xi1", "delta1",
        "mu2", "sigma2", "xi2", "delta2"
      )
    ),
    n = factor(n, levels = c(50, 100, 500))
  )


ggplot(df_all, aes(x = n, y = parameter, fill = BiasShare)) +
  geom_tile(color = "white") +
  facet_wrap(~ theta, ncol = 3, labeller = label_parsed) +
  scale_y_discrete(
    labels = function(x) parse(text = param_labels[x]),
    limits = rev(levels(df_all$parameter))
  ) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1)
  ) +
  labs(
    x = "Sample Size (n)",
    y = NULL,
    fill = "Bias Share"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )








ggplot(df_all, aes(x = n, y = parameter, fill = SEShare)) +
  geom_tile(color = "white") +
  facet_wrap(~ theta, ncol = 3, labeller = label_parsed) +
  scale_y_discrete(
    labels = function(x) parse(text = param_labels[x]),
    limits = rev(levels(df_all$parameter))
  ) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1)
  ) +
  labs(
    x = "Sample Size (n)",
    y = NULL,
    fill = "SE Share"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )































# RMSE HEAT MAP
# -----------------------------
# Load packages
# -----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# -----------------------------
# Convert wide RMSE table to long format
# -----------------------------
df_long <- mc %>%
  pivot_longer(
    cols = starts_with("RMSE_"),
    names_to = "parameter",
    values_to = "RMSE"
  ) %>%
  mutate(
    parameter = gsub("RMSE_", "", parameter),
    
    # 🔥 FORCE correct order here
    parameter = factor(
      parameter,
      levels = c(
        "dep",        # rho first,
        "mu1", "sigma1", "xi1", "delta1",
        "mu2", "sigma2", "xi2", "delta2"
      )
    ),
    
    n = factor(n, levels = c(50, 100, 500))
  )

# -----------------------------
# Create parsed math labels for parameters
# -----------------------------
param_labels <- c(
  mu1     = "mu[1]",
  mu2     = "mu[2]",
  delta1  = "delta[1]",
  delta2  = "delta[2]",
  sigma1  = "sigma[1]",
  sigma2  = "sigma[2]",
  xi1     = "xi[1]",
  xi2     = "xi[2]",
  dep     = "r"   # change if needed
)

# -----------------------------
# Plot RMSE heatmap
# -----------------------------
ggplot(df_long, aes(x = n, y = parameter, fill = RMSE)) +
  geom_tile(color = "white") +
  facet_wrap(~ theta, ncol = 3, labeller = label_parsed) +
  scale_y_discrete(
    labels = function(x) parse(text = param_labels[x]),
    limits = rev(levels(df_long$parameter))
  ) +
  scale_fill_viridis_c(option = "C") +
  labs(
    x = "Sample Size (n)",
    y = NULL,
    fill = "RMSE"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )





# BIAS and SE share 
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# -----------------------------
# Pivot RMSE
# -----------------------------
rmse_long <- mc %>%
  pivot_longer(
    cols = starts_with("RMSE_"),
    names_to = "parameter",
    values_to = "RMSE"
  ) %>%
  mutate(parameter = gsub("RMSE_", "", parameter))

# -----------------------------
# Pivot Bias  (lowercase bias_)
# -----------------------------
bias_long <- mc %>%
  pivot_longer(
    cols = starts_with("bias_"),
    names_to = "parameter",
    values_to = "Bias"
  ) %>%
  mutate(parameter = gsub("bias_", "", parameter))

# -----------------------------
# Pivot SE  (uppercase SE_)
# -----------------------------
se_long <- mc %>%
  pivot_longer(
    cols = starts_with("SE_"),
    names_to = "parameter",
    values_to = "SE"
  ) %>%
  mutate(parameter = gsub("SE_", "", parameter))

# -----------------------------
# Merge and compute shares
# -----------------------------
df_all <- rmse_long %>%
  left_join(bias_long, by = c("theta", "n", "parameter")) %>%
  left_join(se_long,   by = c("theta", "n", "parameter")) %>%
  mutate(
    BiasShare = (Bias^2) / (RMSE^2),
    SEShare   = (SE^2)   / (RMSE^2)
  )

df_all <- df_all %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c(
        "mu1", "mu2",
        "delta1", "delta2",
        "sigma1", "sigma2",
        "xi1", "xi2",
        "dep"
      )
    ),
    n = factor(n, levels = c(50, 100, 500)),
    theta = factor(theta,
                   levels = 1:5,
                   labels = paste0("Theta[", 1:5, "]"))
  )


ggplot(df_all, aes(x = n, y = parameter, fill = BiasShare)) +
  geom_tile(color = "white") +
  facet_wrap(~ theta, ncol = 3, labeller = label_parsed) +
  scale_y_discrete(
    labels = function(x) parse(text = param_labels[x]),
    limits = rev(levels(df_all$parameter))
  ) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1)
  ) +
  labs(
    x = "Sample Size (n)",
    y = NULL,
    fill = "Bias Share"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )








ggplot(df_all, aes(x = n, y = parameter, fill = SEShare)) +
  geom_tile(color = "white") +
  facet_wrap(~ theta, ncol = 3, labeller = label_parsed) +
  scale_y_discrete(
    labels = function(x) parse(text = param_labels[x]),
    limits = rev(levels(df_all$parameter))
  ) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1)
  ) +
  labs(
    x = "Sample Size (n)",
    y = NULL,
    fill = "SE Share"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )






# find out if errors ocurred and NA was returned 
for (f in list.files("benchmarks/individual_results", pattern = ".rds", recursive = TRUE)){
  print(f)
  
  res = readRDS(paste0("benchmarks/individual_results/",f))
  print(paste("how many NAs:", sum(is.na(res$results))))
  
}
