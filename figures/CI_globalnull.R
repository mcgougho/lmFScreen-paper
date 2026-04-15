
# make sure to run renv::activate()!!!!
# (Version guard) Ensure we have the correct lmFScreen
expected <- "0.2.0" 
actual   <- as.character(utils::packageVersion("lmFScreen"))
if (actual != expected) {
  stop(sprintf(
    "lmFScreen version mismatch: expected %s but loaded %s.\nPlease run renv::restore() and re-open the project.",
    expected, actual
  ), call. = FALSE)
}


library(lmFScreen)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

set.seed(1)

# Set parameters
n <- 100
p <- 5
beta <- rep(0,p)
sigma <- 1  
test_col <- 1
B <- 10000
alpha_ov <- 0.05
colors <- c(
  "Selective CI" = "#E69F00",
  "Standard CI" = "#CC79A7"
)
labels <- c(
  "Selective CI " = "Selective~CI",
  "Standard CI" = "Standard~CI"
)


############################# coverage ###################################

n_iter <- 1000

# Sequence of alpha values to test
alpha_seq <- seq(0.01, 0.95, length.out = 10)

# Prepare lists to store confidence intervals for each alpha value
CIs_list <- vector("list", length(alpha_seq))
CIs_naive_list <- vector("list", length(alpha_seq))  # Naïve CI storage

# Loop over each alpha value
for (a in seq_along(alpha_seq)) {
  current_alpha <- alpha_seq[a]
  cat("Running for alpha =", current_alpha, "\n")
  
  # Create matrices to store CIs for the current alpha value
  CIs <- matrix(NA, nrow = n_iter, ncol = 2)
  CIs_naive <- matrix(NA, nrow = n_iter, ncol = 2)  # Naïve CI storage
  
  # Run iterations for the current alpha value
  for (iter in 1:n_iter) {
    if(iter %% 100 == 0) cat("Iteration:", iter, "\n")
    
    repeat {
      X <- matrix(rnorm(n * p), ncol = p)
      y <- X %*% beta + rnorm(n) * sigma
      
      # Project out the intercept
      Xy_info <- lmFScreen:::get_Xy_centered(X,y)
      X <- Xy_info$X
      y <- Xy_info$y
      
      # Compute test statistic
      U <- svd(X)$u
      yPy <- sum((t(U) %*% y)^2)
      rss <- sum(y^2) - yPy
      F_statistic <- (yPy / rss)
      
      # Check if we pass the F-test threshold
      cc <- p/(n-p-1) * qf(1 - alpha_ov, p, (n - p - 1))
      if (F_statistic >= cc) {
        break
      }
    }
    
    # Get confidence intervals using the unknown variance method (with known sigma^2)
    out <- lmFScreen.fit(X, y, test_cols = test_col, alpha = current_alpha, alpha_ov = alpha_ov, B = B)[["selective CIs"]]
    CIs[iter, 1] <- out[test_col, 1]
    CIs[iter, 2] <- out[test_col, 2]
    
    # Compute Naïve 95% Confidence Interval
    lm_naive <- lm(y ~ X + 0)  # Fit linear model without intercept
    beta_hat <- coef(lm_naive)[test_col]
    se_beta <- summary(lm_naive)$coefficients[test_col, 2]  # Standard error
    z_crit <- qnorm(1 - current_alpha / 2)
    
    # Store Naïve CI
    CIs_naive[iter, 1] <- beta_hat - z_crit * se_beta
    CIs_naive[iter, 2] <- beta_hat + z_crit * se_beta
  }
  
  # Save the results for this alpha
  CIs_list[[a]] <- CIs
  CIs_naive_list[[a]] <- CIs_naive  # Store naïve CIs
}

# Function to compute coverage probability
compute_coverage <- function(CI_list, true_value = 0) {
  sapply(CI_list, function(CI_matrix) {
    mean(CI_matrix[, 1] <= true_value & CI_matrix[, 2] >= true_value, na.rm = TRUE)
  })
}

# Compute empirical coverage for all methods
coverage <- compute_coverage(CIs_list)
coverage_naive <- compute_coverage(CIs_naive_list)  # Compute naïve coverage

# Prepare data for plotting
coverage_results <- data.frame(
  nominal_coverage = 1 - alpha_seq,  # Nominal coverage (1 - alpha)
  coverage = coverage, # Empirical coverage for unknown variance method
  coverage_naive = coverage_naive    # Empirical coverage for Naïve method
)

# Reshape data using tidyr::pivot_longer
coverage_results_renamed <- coverage_results %>%
  rename(
    "Selective CI" = coverage,
    "Standard CI" = coverage_naive
  ) %>%
  pivot_longer(-nominal_coverage, names_to = "Method", values_to = "Coverage")


p1 <- ggplot(coverage_results_renamed, aes(x = nominal_coverage, y = Coverage, color = Method)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Reference line y = x
  labs(
    x = expression("Nominal Coverage (1 - " * alpha * ")"),
    y = "Empirical Coverage"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    aspect.ratio = 1,
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Box around plot
  ) +
  scale_color_manual(values = colors, labels = parse(text = labels))

print(p1)

########################## widths, changing beta1 #########################

n_iter <- 1000

beta1_values <- seq(-1, 1, length.out = 40)  # True beta1 values
alpha_ci <- 0.05  # Only compute 95% CIs

# Prepare lists to store CI widths for each beta1 value
widths_lmfscreen_globalnull <- numeric(length(beta1_values))

# Function to compute CI width
compute_CI_width <- function(CI_matrix) {
  mean(CI_matrix[, 2] - CI_matrix[, 1], na.rm = TRUE)
}

# Loop over each true beta1 value
for (b in seq_along(beta1_values)) {
  beta1 <- beta1_values[b]
  
  cat("Running for beta1 =", beta1, "\n")
  
  # Set up beta vector for Global Null case
  beta_globalnull <- rep(0, p)
  beta_globalnull[1] <- beta1
  
  # Create matrices to store CIs
  CIs_lmfscreen_globalnull <- matrix(NA, nrow = n_iter, ncol = 2)
  
  # Run iterations
  for (iter in 1:n_iter) {
    if (iter %% 100 == 0) cat("Iteration:", iter, "\n")
    
    repeat {
      X <- matrix(rnorm(n * p), ncol = p)
      y_globalnull <- X %*% beta_globalnull + rnorm(n) * sigma
      
      # Project out the intercept
      Xy_info <- lmFScreen:::get_Xy_centered(X,y_globalnull)
      X <- Xy_info$X
      y_globalnull <- Xy_info$y
      
      # Compute test statistic
      U <- svd(X)$u
      yPy <- sum((t(U) %*% y_globalnull)^2)
      rss <- sum(y_globalnull^2) - yPy
      F_statistic <- (yPy / rss)
      
      # Check if we pass the F-test threshold
      cc <- p / (n - p - 1) * qf(1 - alpha_ov, p, (n - p - 1))
      if (F_statistic >= cc) {
        break
      }
    }
    
    # lmFScreen CI (default)
    CIs_lmfscreen_globalnull[iter, ] <- lmFScreen.fit(X, y_globalnull, test_cols = test_col,
                                                      alpha = alpha_ci, alpha_ov = alpha_ov, B = B)[["selective CIs"]][test_col, ]
    
  }
  
  # Compute average CI widths
  widths_lmfscreen_globalnull[b] <- compute_CI_width(CIs_lmfscreen_globalnull)
}


# Prepare data
widths_globalnull <- data.frame(
  beta1 = beta1_values,
  lmfscreen = widths_lmfscreen_globalnull,
  naive = NA
)

widths_globalnull_renamed <- widths_globalnull %>%
  rename(
    "Selective CI" = lmfscreen,
    "Standard CI" = naive
  ) %>%
  pivot_longer(-beta1, names_to = "Method", values_to = "Width")

p2 <- ggplot(widths_globalnull_renamed, aes(x = beta1, y = Width, color = Method)) +
  geom_line(linewidth = 1.2) +
  labs(
    x = expression("True " * beta[1]),
    y = "Average 95% CI Width"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    aspect.ratio = 1,
    legend.position = "none",  # Keep legend only in this plot
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_color_manual(values = colors, labels = parse(text = labels))

print(p2)

######################### widths, changing n ################################

n_values <- seq(50, 500, length.out = 10)

# Prepare lists to store CI widths for each sample size
widths_lmfscreen_n <- numeric(length(n_values))

# Loop over each sample size
for (i in seq_along(n_values)) {
  n <- n_values[i]
  
  cat("Running for sample size n =", n, "\n")
  
  # Set up beta vector for Global Null case
  beta_null <- rep(0, p)
  
  # Create matrices to store CIs
  CIs_lmfscreen_n <- matrix(NA, nrow = n_iter, ncol = 2)
  
  # Run iterations
  for (iter in 1:n_iter) {
    if (iter %% 100 == 0) cat("Iteration:", iter, "\n")
    
    repeat {
      X <- matrix(rnorm(n * p), ncol = p)
      y <- X %*% beta_null + rnorm(n) * sigma  # Under the global null
      
      # Project out the intercept
      Xy_info <- lmFScreen:::get_Xy_centered(X,y)
      X <- Xy_info$X
      y <- Xy_info$y
      
      # Compute test statistic
      U <- svd(X)$u
      yPy <- sum((t(U) %*% y)^2)
      rss <- sum(y^2) - yPy
      F_statistic <- (yPy / rss)
      
      # Check if we pass the F-test threshold
      cc <- p / (n - p - 1) * qf(1 - alpha_ov, p, (n - p - 1))
      if (F_statistic >= cc) {
        break
      }
    }
    
    # lmFScreen CI (default)
    CIs_lmfscreen_n[iter, ] <- lmFScreen.fit(X, y, test_cols = test_col,
                                             alpha = alpha_ci, alpha_ov = alpha_ov, B = B)[["selective CIs"]][test_col, ]
    
  }
  
  # Compute average CI widths
  widths_lmfscreen_n[i] <- compute_CI_width(CIs_lmfscreen_n)
}

widths_n <- data.frame(
  n = n_values,
  lmfscreen = widths_lmfscreen_n,
  naive = NA
)

widths_n_renamed <- widths_n %>%
  rename(
    "Selective CI" = lmfscreen,
    "Standard CI" = naive
  ) %>%
  pivot_longer(-n, names_to = "Method", values_to = "Width")

p3 <- ggplot(widths_n_renamed, aes(x = n, y = Width, color = Method)) +
  geom_line(linewidth = 1.2) +
  labs(
    x = expression("Sample Size " * n),
    y = "Average 95% CI Width"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Box around plot
    aspect.ratio = 1  # Make plot square
  ) +
  scale_color_manual(values = colors, labels = parse(text = labels))

print(p3)

################################ final plot ######################################

# Combine plots into one figure with a shared legend
final_plot <- p1 + p2 + p3 +
  plot_layout(guides = "collect") &
  theme(
    text = element_text(family = "Helvetica"),
    labels = parse(text = labels),
    legend.position = "bottom",
    legend.background = element_rect(color = "black", fill = "white", size = 0.8),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(5, "mm"),
    legend.text = element_text(size = 10),
    legend.title = element_blank()
  )

# Display final plot
print(final_plot)

pdf("CI_globalnull.pdf", width = 8, height = 3)
print(final_plot)
dev.off()






