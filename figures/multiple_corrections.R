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
library(dplyr)
library(tidyr)
library(patchwork)

# Simulation parameters
nreps <- 1e5
n_per_group <- 30
n_group <- 4
sigma <- 1

# Pre-allocate vectors for the p-values:
p.overall <- rep(NA, nreps)  # overall F-test p-value (group effect)
p.1n <- rep(NA, nreps)       # unadjusted p-value for group B (vs. baseline A)
p.1b <- rep(NA, nreps)       # Bonferroni-adjusted (×2)
p.1b3 <- rep(NA, nreps)      # Bonferroni-adjusted (×3)
p.1s <- rep(NA, nreps)       # Scheffé-adjusted
p.1closed <- rep(NA, nreps)  # Closed testing
psel.1 <- rep(NA, nreps)     # Selection-adjusted

# make design matrix
group <- factor(rep(letters[1:n_group], each = n_per_group))
X <- model.matrix(~ group)
X.1 <- X[,-1, drop = FALSE]
beta <- rep(0, n_group)

for (i in 1:nreps) {
  # Generate the group factor and design matrix X (baseline is group A)
  y <- X %*% beta + rnorm(n_group * n_per_group) * sigma
  
  # Fit the linear model using lm:
  lm_model <- lm(y ~ X + 0)
  
  # Overall test: compare the full model to the null model (only intercept)
  lm_null <- lm(y ~ 1)
  anova_out <- anova(lm_null, lm_model)
  # Extract the p-value for the additional predictors (i.e. the group effect)
  p.overall[i] <- anova_out$`Pr(>F)`[2]
  
  # Extract the summary to get the coefficient for groupB
  lm_summary <- summary(lm_model)
  # The coefficient named "XgroupB" compares group B vs. baseline A.
  p_raw <- lm_summary$coefficients["Xgroupb", "Pr(>|t|)"]
  p.1n[i] <- p_raw
  
  # Bonferroni correction for 2 and 3 tests (overall and pairwise)
  p.1b[i] <- ifelse(p_raw <= 1/2, 2 * p_raw, 1)
  p.1b3[i] <- ifelse(p_raw <= 1/3, 3 * p_raw, 1)
  
  # Scheffé adjustment:
  # 1. Extract the t-statistic for the groupB coefficient.
  t_stat <- lm_summary$coefficients["Xgroupb", "t value"]
  # 2. Convert the t-statistic to an F-statistic:
  F_stat <- t_stat^2 / (n_group - 1)  # here k = n_group
  # 3. Calculate the residual degrees of freedom:
  df_error <- n_per_group * n_group - ncol(X)
  # 4. Compute the Scheffé p-value from the F distribution:
  p.1s[i] <- 1 - pf(F_stat, df1 = n_group - 1, df2 = df_error)
  
  # Closed testing p-value
  p.1closed[i] <- max(p.overall[i], p_raw)
  
  # Compute the Selective p-value:
  Xy_centered <- lmFScreen:::get_Xy_centered(X.1,y)
  X_centered <- Xy_centered$X
  y_centered <- Xy_centered$y
  pselb <- lmFScreen:::get_pselb(X=X_centered,y=y_centered,sigma_sq=sigma^2,alpha_ov=0.05,B=100000,min_select=10000, seed = 123*i )
  psel.1[i] <- pselb(0)
}



# Define quantiles for empirical CDF plots:
qs <- seq(0, 1, length.out = 1000)
theorqs <- qunif(qs)

# Compute empirical CDF quantiles for the three sets of p-values:
emp.cdf.n <- ecdf(p.1n)
empqs.n <- quantile(emp.cdf.n, qs)

emp.cdf.b <- ecdf(p.1b)
empqs.b <- quantile(emp.cdf.b, qs)

emp.cdf.b3 <- ecdf(p.1b3)
empqs.b3 <- quantile(emp.cdf.b3, qs)

emp.cdf.s <- ecdf(p.1s)
empqs.s <- quantile(emp.cdf.s, qs)

emp.cdf.closed <- ecdf(p.1closed)
empqs.closed <- quantile(emp.cdf.closed, qs)

emp.cdf.sel <- ecdf(psel.1)
empqs.sel <- quantile(emp.cdf.sel, qs)

mask <- p.overall < (0.05)

emp.cdf.n.mask <- ecdf(p.1n[mask])
empqs.n.mask <- quantile(emp.cdf.n.mask, qs)

emp.cdf.b.mask <- ecdf(p.1b[mask])
empqs.b.mask <- quantile(emp.cdf.b.mask, qs)

emp.cdf.b3.mask <- ecdf(p.1b3[mask])
empqs.b3.mask <- quantile(emp.cdf.b3.mask, qs)

emp.cdf.s.mask <- ecdf(p.1s[mask])
empqs.s.mask <- quantile(emp.cdf.s.mask, qs)

emp.cdf.closed.mask <- ecdf(p.1closed[mask])
empqs.closed.mask <- quantile(emp.cdf.closed.mask, qs)

emp.cdf.sel.mask <- ecdf(psel.1[mask])
empqs.sel.mask <- quantile(emp.cdf.sel.mask, qs)


# Convert empirical distributions to data frames for ggplot
df_unconditional <- data.frame(
  Theoretical = theorqs,
  No_Correction = empqs.n,
  Bonferroni_2 = empqs.b,
  Bonferroni_3 = empqs.b3,
  Scheffe = empqs.s,
  Closed = empqs.closed,
  Selective = empqs.sel
) %>%
  pivot_longer(cols = -Theoretical, names_to = "Method", values_to = "Empirical") %>%
  mutate(Method = factor(Method, levels = c("No_Correction", "Bonferroni_2", "Bonferroni_3", "Scheffe", "Closed", "Selective")))

df_conditional <- data.frame(
  Theoretical = theorqs,
  No_Correction = empqs.n.mask,
  Bonferroni_2 = empqs.b.mask,
  Bonferroni_3 = empqs.b3.mask,
  Scheffe = empqs.s.mask,
  Closed = empqs.closed.mask,
  Selective = empqs.sel.mask
) %>%
  pivot_longer(cols = -Theoretical, names_to = "Method", values_to = "Empirical") %>%
  mutate(Method = factor(Method, levels = c("No_Correction", "Bonferroni_2", "Bonferroni_3", "Scheffe", "Closed", "Selective")))

colors <- c(
  "No_Correction" = "#000000",    # Black
  "Bonferroni_2" = "#E69F00",     # Orange
  "Bonferroni_3" = "#009E73",     # Teal/Green
  "Scheffe" = "#0072B2",          # Blue
  "Closed" = "#F0E442",           # Yellow
  "Selective" = "#CC79A7"         # Purple/Magenta
)

labels <- c(
  "No_Correction" = expression(p[H[0]^M] ~ "in (5)"),
  "Bonferroni_2" = expression(2 * p[H[0]^M]),
  "Bonferroni_3" = expression(3 * p[H[0]^M]),
  "Scheffe" = "Scheffé",
  "Closed" = "Closed",
  "Selective" = "Selective"
)


# Base theme modifications
base_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Box around plot
    legend.position = "none",  # Remove legend from individual plots
    aspect.ratio = 1,  # Ensure square plots
    plot.title = element_blank(),  # Remove plot title
    axis.title = element_text(size = 10)   # Keep axis labels but make them smaller
  )

# Create the first plot (Unconditional p-values)
p1 <- ggplot(df_unconditional, aes(x = Theoretical, y = Empirical, color = Method)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = colors, labels = labels) +
  labs(x = "Theoretical Quantiles", y = "Empirical Quantiles") +
  base_theme

# Create the second plot (Conditional p-values)
p2 <- ggplot(df_conditional, aes(x = Theoretical, y = Empirical, color = Method)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = colors, labels = labels) +
  labs(x = "Theoretical Quantiles", y = "Empirical Quantiles") +
  base_theme

# Arrange side-by-side with a single common legend
final_plot <- p1 + p2 +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.background = element_rect(color = "black", fill = "white", size = 0.8),  # Box around legend
    legend.key = element_rect(fill = "black"),  # Ensure color boxes are filled
    legend.key.size = unit(5, "mm"),  # Reduce overall size but keep visibility
    legend.text = element_text(size = 9),  # Reduce legend text size
    legend.spacing.x = unit(2, "mm"),  # Reduce spacing between legend items
    legend.margin = margin(3, 3, 3, 3),  # Compact legend margin
    legend.title = element_blank()  # Remove "Method" label from legend
  ) &
  guides(color = guide_legend(
    override.aes = list(fill = colors, shape = 22, size = 5)  # Force legend to use filled squares
  ))

# Display the final plot
print(final_plot)

pdf("multiple_corrections.pdf", width = 8, height = 4)
print(final_plot)
dev.off()

