# make sure to run renv::activate()!!!!
# (Version guard) Ensure we have the correct lmFScreen
expected <- "0.1.0"   
actual   <- as.character(utils::packageVersion("lmFScreen"))
if (actual != expected) {
  stop(sprintf(
    "lmFScreen version mismatch: expected %s but loaded %s.\nPlease run renv::restore() and re-open the project.",
    expected, actual
  ), call. = FALSE)
}


library(lmFScreen)
library(ggplot2)
library(ggExtra)

# set parameters
seed <- 112233
n <- 10
p <- 2
b_true <- 0
sigma <- 2
beta <- rep(0, p)
beta[1] <- b_true
alpha <- 0.05
alpha_ov <- 0.05
B <- 10000
F.quantile <- qf(1 - alpha, p, (n - p - 1))
cc <- p/(n-p-1) * F.quantile
test_cols <- 1

################################# oracle variance ############################

# We will store the conditional and unconditional MLEs
beta_mle_vals_cond <- beta_mle_vals_uncond <- numeric(B)

for (rep_idx in 1:B) {
  # generate data that passes F-screening
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
    if (F_statistic >= cc) {
      break
    }
    # Otherwise, generate a new (X, y)
  }

  # get unconditional MLE (closed-form)
  naive_lm_coefs <- summary(lm(y ~ X + 0))$coefficients[test_cols, , drop = FALSE]
  naive_point_est <- naive_lm_coefs[1, 1]
  beta_mle_vals_uncond[rep_idx] <- naive_point_est

  # get interval to search for conditional MLE in
  naive_se_est <- naive_lm_coefs[1, 2]
  interval <- c(naive_point_est - 10 * naive_se_est, naive_point_est + 10 * naive_se_est)

  # get conditional MLE
  point_est <- as.numeric(lmFScreen:::compute_MLE(X, y, sigma_sq = sigma^2, interval = interval, seed = seed*rep_idx, alpha_ov = alpha_ov, B = B)[[1]])
  beta_mle_vals_cond[rep_idx] <- point_est
}

# Find limits to ensure equal scaling
xlim_range <- range(beta_mle_vals_uncond)
ylim_range <- range(beta_mle_vals_cond)
max_range <- range(c(xlim_range, ylim_range))  # Get common range

# Create scatter plot
scatter_plot <- ggplot(data = data.frame(beta_mle_vals_cond, beta_mle_vals_uncond),
                       aes(x = beta_mle_vals_uncond, y = beta_mle_vals_cond)) +
  geom_point(color = "blue", alpha = 0.1, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Unconditional MLE", y = "Conditional MLE", title = " ") +
  coord_fixed(ratio = 1, xlim = max_range, ylim = max_range)  # Set equal limits

# Add marginal histograms
ggExtra::ggMarginal(scatter_plot, type = "histogram", fill = "gray", bins = 35)


######################### with debiased variance estimate ############################

# We will store the conditional and unconditional MLEs
beta_mle_vals_cond2 <- beta_mle_vals_uncond2 <- numeric(B)

for (rep_idx in 1:B) {
  # generate data that passes F-screening
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
    if (F_statistic >= cc) {
      break
    }
    # Otherwise, generate a new (X, y)
  }

  # get unconditional MLE (closed-form)
  naive_lm_coefs <- summary(lm(y ~ X + 0))$coefficients[test_cols, , drop = FALSE]
  naive_point_est <- naive_lm_coefs[1, 1]
  beta_mle_vals_uncond2[rep_idx] <- naive_point_est

  # get variance estimate
  sigma_sq <- lmFScreen:::get_variance_estimate(X,y,alpha_ov)

  # get interval to search for conditional MLE in
  naive_se_est <- naive_lm_coefs[1, 2]
  interval <- c(naive_point_est - 10 * naive_se_est, naive_point_est + 10 * naive_se_est)

  # get conditional MLE
  point_est <- as.numeric(lmFScreen:::compute_MLE(X, y, sigma_sq = sigma_sq, interval = interval, seed = seed*rep_idx, alpha_ov = alpha_ov, B = B)[[1]])
  beta_mle_vals_cond2[rep_idx] <- point_est
}

# Find limits to ensure equal scaling
xlim_range <- range(beta_mle_vals_uncond2)
ylim_range <- range(beta_mle_vals_cond2)
max_range <- range(c(xlim_range, ylim_range))  # Get common range

# Create scatter plot
scatter_plot <- ggplot(data = data.frame(beta_mle_vals_cond2, beta_mle_vals_uncond2),
                       aes(x = beta_mle_vals_uncond2, y = beta_mle_vals_cond2)) +
  geom_point(color = "blue", alpha = 0.1, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Unconditional MLE", y = "Conditional MLE", title = " ") +
  coord_fixed(ratio = 1, xlim = max_range, ylim = max_range)  # Set equal limits

# Add marginal histograms
ggExtra::ggMarginal(scatter_plot, type = "histogram", fill = "gray", bins = 35)


################################# make final plot #############################################

# Find global axis limits for both plots
xlim_range <- range(c(beta_mle_vals_uncond, beta_mle_vals_uncond2))
ylim_range <- range(c(beta_mle_vals_cond, beta_mle_vals_cond2))
max_range <- range(c(xlim_range, ylim_range))  # Ensure common scale for both axes

# Function to create scatter plot (without `ggMarginal`)
create_scatter <- function(mle_vals_uncond, mle_vals_cond) {
  ggplot(data = data.frame(mle_vals_uncond, mle_vals_cond),
         aes(x = mle_vals_uncond, y = mle_vals_cond)) +
    geom_point(color = "blue", alpha = 0.1, size = 1) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    labs(x = "Unconditional MLE", y = "Conditional MLE", title = " ") +
    coord_fixed(ratio = 1, xlim = max_range, ylim = max_range)  # Set same limits
}

# Create scatter plots
plot1 <- create_scatter(beta_mle_vals_uncond, beta_mle_vals_cond)
plot2 <- create_scatter(beta_mle_vals_uncond2, beta_mle_vals_cond2)

# Convert ggMarginal() outputs into ggplot-compatible objects
p1_marginal <- as.ggplot(ggMarginal(plot1, type = "histogram", fill = "gray", bins = 35))
p2_marginal <- as.ggplot(ggMarginal(plot2, type = "histogram", fill = "gray", bins = 35))

# Arrange the two plots side by side
final_plot <-  p2_marginal + p1_marginal + plot_layout(nrow = 1)

# Display the final combined plot
print(final_plot)

