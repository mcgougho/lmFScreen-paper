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
library(ggExtra)

seed <- 112233
n <- 50
p <- 2
b_true <- 0
sigma <- 2
beta <- c(b_true, rep(0, p-1))
alpha_ov <- 0.05
B_rep <- 10000
F.quantile <- qf(1 - alpha_ov, p, (n - p - 1))
cc <- p/(n-p-1) * F.quantile

beta_mle_vals_uncond <- beta_mle_vals_cond <- numeric(B_rep)
sigma_sq_mle_vals_uncond <- sigma_sq_mle_vals_cond <- numeric(B_rep)

for (rep_idx in 1:B_rep) {
  
  repeat {
    X <- matrix(rnorm(n * p), ncol = p)
    y <- X %*% beta + rnorm(n) * sigma
    Xy_info <- lmFScreen:::get_Xy_centered(X, y)
    X <- Xy_info$X
    y <- Xy_info$y
    U <- svd(X)$u
    yPy <- sum((t(U) %*% y)^2)
    rss <- sum(y^2) - yPy
    F_statistic <- yPy / rss
    if (F_statistic >= cc) break
  }
  
  fit_naive <- lm(y ~ X + 0)
  coefs <- summary(fit_naive)$coefficients[1, , drop = FALSE]
  beta_hat <- coefs[1, 1]
  se_hat   <- coefs[1, 2]
  sigma_sq_hat <- sum(resid(fit_naive)^2) / (nrow(X)-p)
  interval_beta <- c(beta_hat - 10 * se_hat, beta_hat + 10 * se_hat)
  
  beta_mle_vals_uncond[rep_idx]  <- beta_hat
  sigma_sq_mle_vals_uncond[rep_idx] <- sigma_sq_hat
  
  interval_beta <- c(beta_hat - 10 * se_hat, beta_hat + 10 * se_hat)
  interval_sigma_sq <- c(sigma_sq_hat / 10, sigma_sq_hat * 10)
  
  est <- compute_MLE(
    X, y, 
    alpha_ov = alpha_ov,
    beta_init = beta_hat,
    sigma_sq_init = sigma_sq_hat,
    interval_beta = interval_beta,
    interval_sigma_sq = interval_sigma_sq,
    B = 10000
  )
  
  beta_mle_vals_cond[rep_idx]   <- est$beta1
  sigma_sq_mle_vals_cond[rep_idx] <- est$sigma_sq
}

# beta plot
xlim_range <- range(beta_mle_vals_uncond)
ylim_range <- range(beta_mle_vals_cond)
max_range <- range(c(xlim_range, ylim_range))  # Get common range

# Create scatter plot
beta_scatter_plot <- ggplot(data = data.frame(beta_mle_vals_cond, beta_mle_vals_uncond),
                            aes(x = beta_mle_vals_uncond, y = beta_mle_vals_cond)) +
  geom_point(color = "blue", alpha = 0.1, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Unconditional MLE", y = "Conditional MLE", title = " ") +
  coord_fixed(ratio = 1, xlim = max_range, ylim = max_range)  # Set equal limits

# Add marginal histograms
ggExtra::ggMarginal(beta_scatter_plot, type = "histogram", fill = "gray", bins = 35)


# sigma^2 plot
xlim_range <- range(sigma_sq_mle_vals_uncond)
ylim_range <- range(sigma_sq_mle_vals_cond)
max_range <- range(c(xlim_range, ylim_range))  # Get common range

# Create scatter plot
sigma_sq_scatter_plot <- ggplot(data = data.frame(sigma_sq_mle_vals_cond, sigma_sq_mle_vals_uncond),
                                aes(x = sigma_sq_mle_vals_uncond, y = sigma_sq_mle_vals_cond)) +
  geom_point(color = "blue", alpha = 0.1, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Unconditional MLE", y = "Conditional MLE", title = " ") +
  coord_fixed(ratio = 1, xlim = max_range, ylim = max_range)  # Set equal limits

# Add marginal histograms
ggExtra::ggMarginal(sigma_sq_scatter_plot, type = "histogram", fill = "gray", bins = 35)



################################# non-zero beta #############################################


b_true <- 0.5
beta <- c(b_true, rep(0, p-1))

beta_mle_vals_uncond.5 <- beta_mle_vals_cond.5 <- numeric(B_rep)
sigma_sq_mle_vals_uncond.5 <- sigma_sq_mle_vals_cond.5 <- numeric(B_rep)

for (rep_idx in 1:B_rep) {
  
  repeat {
    X <- matrix(rnorm(n * p), ncol = p)
    y <- X %*% beta + rnorm(n) * sigma
    Xy_info <- lmFScreen:::get_Xy_centered(X, y)
    X <- Xy_info$X
    y <- Xy_info$y
    U <- svd(X)$u
    yPy <- sum((t(U) %*% y)^2)
    rss <- sum(y^2) - yPy
    F_statistic <- yPy / rss
    if (F_statistic >= cc) break
  }
  
  fit_naive <- lm(y ~ X + 0)
  coefs <- summary(fit_naive)$coefficients[1, , drop = FALSE]
  beta_hat <- coefs[1, 1]
  se_hat   <- coefs[1, 2]
  sigma_sq_hat <- sum(resid(fit_naive)^2) / (nrow(X)-p)
  interval_beta <- c(beta_hat - 10 * se_hat, beta_hat + 10 * se_hat)
  
  beta_mle_vals_uncond.5[rep_idx]  <- beta_hat
  sigma_sq_mle_vals_uncond.5[rep_idx] <- sigma_sq_hat
  
  interval_beta <- c(beta_hat - 10 * se_hat, beta_hat + 10 * se_hat)
  interval_sigma_sq <- c(sigma_sq_hat / 10, sigma_sq_hat * 10)
  
  est <- compute_MLE(
    X, y, 
    alpha_ov = alpha_ov,
    beta_init = beta_hat,
    sigma_sq_init = sigma_sq_hat,
    interval_beta = interval_beta,
    interval_sigma_sq = interval_sigma_sq,
    B = 10000
  )
  
  beta_mle_vals_cond.5[rep_idx]   <- est$beta1
  sigma_sq_mle_vals_cond.5[rep_idx] <- est$sigma_sq
}

# beta plot
xlim_range <- range(beta_mle_vals_uncond.5)
ylim_range <- range(beta_mle_vals_cond.5)
max_range <- range(c(xlim_range, ylim_range))  # Get common range

# Create scatter plot
beta_scatter_plot.5 <- ggplot(data = data.frame(beta_mle_vals_cond.5, beta_mle_vals_uncond.5),
                              aes(x = beta_mle_vals_uncond.5, y = beta_mle_vals_cond.5)) +
  geom_point(color = "blue", alpha = 0.1, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Unconditional MLE", y = "Conditional MLE", title = " ") +
  coord_fixed(ratio = 1, xlim = max_range, ylim = max_range)  # Set equal limits

# Add marginal histograms
ggExtra::ggMarginal(beta_scatter_plot, type = "histogram", fill = "gray", bins = 35)


# sigma^2 plot
xlim_range <- range(sigma_sq_mle_vals_uncond.5)
ylim_range <- range(sigma_sq_mle_vals_cond.5)
max_range <- range(c(xlim_range, ylim_range))  # Get common range

# Create scatter plot
sigma_sq_scatter_plot.5 <- ggplot(data = data.frame(sigma_sq_mle_vals_cond.5, sigma_sq_mle_vals_uncond.5),
                                  aes(x = sigma_sq_mle_vals_uncond.5, y = sigma_sq_mle_vals_cond.5)) +
  geom_point(color = "blue", alpha = 0.1, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Unconditional MLE", y = "Conditional MLE", title = " ") +
  coord_fixed(ratio = 1, xlim = max_range, ylim = max_range)  # Set equal limits

# Add marginal histograms
ggExtra::ggMarginal(sigma_sq_scatter_plot, type = "histogram", fill = "gray", bins = 35)





################################# make final plot #############################################

library(patchwork) 

# Convert ggMarginal() outputs into ggplot-compatible objects
p1_marginal <- ggExtra::ggMarginal(beta_scatter_plot, type = "histogram", fill = "gray", bins = 35)
p2_marginal <- ggExtra::ggMarginal(beta_scatter_plot.5, type = "histogram", fill = "gray", bins = 35)

# Arrange the two plots side by side
final_plot <- wrap_elements(p1_marginal) + 
  wrap_elements(p2_marginal) +
  plot_layout(nrow = 1)

# Display the final combined plot
print(final_plot)

pdf("point_est.pdf", width = 6, height = 3)
print(final_plot)
dev.off()




