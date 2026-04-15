library(lmFScreen)
library(ggplot2)
library(patchwork)


n <- 30
p <- 10
sigma <- 1
nreps <- 10000
beta <- rep(0, p)
alpha_levels <- c(0.01, 0.05, 0.10)


##################### overall test ############################

pvals_F <- numeric(nreps)
pvals_chisq_plugin <- numeric(nreps)
pvals_chisq_conserv <- numeric(nreps)
pvals_chisq_oracle <- numeric(nreps)

for(iter in 1:nreps){
  
  X <- matrix(rnorm(n * p), ncol = p)
  y <- X %*% beta + rnorm(n)*sigma
  
  # center X and y
  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- as.numeric(scale(y, center = TRUE, scale = FALSE))
  
  U <- svd(Xc)$u
  yPy <- sum((crossprod(U, yc))^2)
  rss <- sum(yc^2) - yPy
  
  # estimated sigma^2
  sigmahat_sq <- rss / (n - p - 1)
  sigmahat_sq_conserv <- sum((yc - mean(yc))^2)/(n-1)
  
  # exact F test
  Fstat <- (yPy / p) / sigmahat_sq
  pvals_F[iter] <- 1 - pf(Fstat, df1 = p, df2 = n - p - 1)
  
  # chi-square test
  chisq_plugin_stat <- yPy / sigmahat_sq
  pvals_chisq_plugin[iter] <- 1 - pchisq(chisq_plugin_stat, df = p)
  chisq_conserv_stat <- yPy / sigmahat_sq_conserv
  pvals_chisq_conserv[iter] <- 1 - pchisq(chisq_conserv_stat, df = p)
  chisq_oracle_stat <- yPy / sigma^2
  pvals_chisq_oracle[iter] <- 1 - pchisq(chisq_oracle_stat, df = p)
  
}

par(pty = "s")

plot(pvals_F, pvals_chisq_oracle,
     xlab = "F-statistic p-values",
     ylab = "Chi-squared p-values",
     xlim = c(0,1), ylim = c(0,1),
     pch = 16, cex = 0.5, col = "lightpink")
points(pvals_F, pvals_chisq_conserv,
       pch = 16, cex = 0.5, col = "#0072B2")
points(pvals_F, pvals_chisq_plugin,
       pch = 16, cex = 0.5, col = "#009E73")
abline(0, 1, col = "red", lwd =2)
abline(h = 0.05, col = "red", lty = 2, lwd =2)
abline(v = 0.05, col = "red", lty = 2, lwd =2)


############################ post hoc test ###############################

simulate_data <- function(n, p, beta, sigma){
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  y <- X %*% beta + rnorm(n, sd=sigma)
  list(X=X, y=y)
}

# truncated normal CDF on (-Inf,B] U [A,Inf)
ptrunc_union <- function(x, mu, sd, B, A){
  
  if(!is.finite(A) && !is.finite(B))
    return(pnorm((x-mu)/sd))
  
  zB <- (B-mu)/sd
  zA <- (A-mu)/sd
  zx <- (x-mu)/sd
  
  # denominator: P(T in (-Inf,B] U [A,Inf))
  denom <- pnorm(zB) + (1 - pnorm(zA))
  
  if(x <= B){
    numer <- pnorm(zx)
  } else if(x < A){
    numer <- pnorm(zB) # flat over the gap (B, A)
  } else{
    numer <- pnorm(zB) + (pnorm(zx) - pnorm(zA))
  }
  
  numer / denom
}

compute_trunc_pval <- function(beta_hat, XtX_inv, sigma2_hat, alpha, p){
  
  Sigma_hat <- sigma2_hat * XtX_inv
  K_hat <- solve(Sigma_hat)
  
  s0 <- qchisq(1-alpha, df=p)
  
  Tobs <- as.numeric(beta_hat[1])
  Sigma11 <- Sigma_hat[1,1]
  
  cvec <- Sigma_hat[,1,drop=FALSE] / Sigma11
  W <- beta_hat - cvec * Tobs
  
  a <- as.numeric(t(cvec)%*%K_hat%*%cvec)
  b <- as.numeric(t(W)%*%K_hat%*%cvec)
  d <- as.numeric(t(W)%*%K_hat%*%W - s0)
  
  Delta <- 4*b^2 - 4*a*d
  
  if(Delta < 0){
    B <- Inf
    A <- - Inf
  } else{
    B <- (-2*b - sqrt(Delta))/(2*a)
    A <- (-2*b + sqrt(Delta))/(2*a)
  }
  
  sdT <- sqrt(Sigma11)
  cdf_val <- ptrunc_union(Tobs, 0, sdT, B, A)
  p_right <- 1 - cdf_val
  pval <- 2 * min(cdf_val, 1 - cdf_val)
  min(pval, 1)
}

compute_trunc_pval_onesided <- function(beta_hat, XtX_inv, sigma2_hat, alpha, p){
  
  Sigma_hat <- sigma2_hat * XtX_inv
  K_hat <- solve(Sigma_hat)
  
  s0 <- qchisq(1-alpha, df=p)
  
  Tobs <- as.numeric(beta_hat[1])
  Sigma11 <- Sigma_hat[1,1]
  
  cvec <- Sigma_hat[,1,drop=FALSE] / Sigma11
  W <- beta_hat - cvec * Tobs
  
  a <- as.numeric(t(cvec)%*%K_hat%*%cvec)
  b <- as.numeric(t(W)%*%K_hat%*%cvec)
  d <- as.numeric(t(W)%*%K_hat%*%W - s0)
  
  Delta <- 4*b^2 - 4*a*d
  
  if(Delta < 0){
    B <- -Inf
    A <- Inf
  } else{
    B <- (-2*b - sqrt(Delta))/(2*a)
    A <- (-2*b + sqrt(Delta))/(2*a)
  }
  
  sdT <- sqrt(Sigma11)
  cdf_val <- ptrunc_union(Tobs, 0, sdT, B, A)
  1 - cdf_val
}


simulate_trunc <- function(nsim=5000, n=30, p=10,
                           beta=rep(0,p), sigma=1, alpha=0.05){
  
  pvals_heller <- numeric(nsim)
  pvals_heller_onesided <- numeric(nsim)
  pvals_heller_oracle <- numeric(nsim)
  pvals_heller_conserv <- numeric(nsim)
  pvals_heller_switch <- numeric(nsim)
  pvals_fscreen <- numeric(nsim)
  i <- 1
  F_cut <- qf(1 - alpha, df1 = p, df2 = n - p - 1)
  
  while(i <= nsim){
    
    dat <- simulate_data(n,p,beta,sigma)
    centered <- lmFScreen:::get_Xy_centered(dat$X, dat$y)
    
    X <- centered$X
    y <- centered$y
    
    XtX <- crossprod(X)
    XtX_inv <- solve(XtX)
    
    beta_hat <- XtX_inv %*% crossprod(X,y)
    
    rss <- sum((y - X %*% beta_hat)^2)
    sigma2_hat <- rss/(n-p-1)
    
    Sigma_hat <- sigma2_hat * XtX_inv
    K_hat <- solve(Sigma_hat)
    
    S <- as.numeric(t(beta_hat) %*% K_hat %*% beta_hat)
    F_stat <- S / p
    
    sigma2_hat_v2 <- sum((y - mean(y))^2)/(n-1)
    
    if(F_stat >= F_cut){
      alpha2 <- 0.05^3
      if(F_stat >=  qf(1 - alpha2, df1 = p, df2 = n - p - 1)){
        pvals_heller_switch[i] <- compute_trunc_pval(beta_hat, XtX_inv, sigma2_hat, alpha, p)
      } else{
        pvals_heller_switch[i] <- compute_trunc_pval(beta_hat, XtX_inv, sigma2_hat_v2, alpha, p)
      }
      pvals_heller[i] <- compute_trunc_pval(beta_hat, XtX_inv, sigma2_hat, alpha, p)
      pvals_heller_onesided[i] <- compute_trunc_pval_onesided(beta_hat, XtX_inv, sigma2_hat, alpha, p)
      pvals_heller_conserv[i] <- compute_trunc_pval(beta_hat, XtX_inv, sigma2_hat_v2, alpha, p)
      pvals_heller_oracle[i] <- compute_trunc_pval(beta_hat, XtX_inv, sigma^2, alpha, p)
      pselb <- lmFScreen:::get_pselb(X=X,y=y,alpha_ov=0.05)
      pvals_fscreen[i] <- pselb(0)
      i <- i + 1
    }
  }
  
  return(list(pvals_heller = pvals_heller, pvals_heller_onesided = pvals_heller_onesided, pvals_fscreen = pvals_fscreen, pvals_heller_conserv = pvals_heller_conserv, pvals_heller_oracle = pvals_heller_oracle, pvals_heller_switch = pvals_heller_switch))
}

# run simulation
pvals <- simulate_trunc()
pvals_heller_onesided <- pvals$pvals_heller_onesided
pvals_heller <- pvals$pvals_heller
pvals_heller_conserv <- pvals$pvals_heller_conserv
pvals_heller_oracle <- pvals$pvals_heller_oracle
pvals_heller_switch <- pvals$pvals_heller_switch
pvals_fscreen <- pvals$pvals_fscreen

qs <- seq(from = 0, to = 1, length.out = 1000)
theorqs <- qunif(qs)
emp.cdf.tn <- ecdf(pvals_heller)
empqs.tn <- quantile(emp.cdf.tn, qs)
emp.cdf.tn.onesided <- ecdf(pvals_heller_onesided)
empqs.tn.onesided <- quantile(emp.cdf.tn.onesided, qs)
emp.cdf.fs <- ecdf(pvals_fscreen)
empqs.fs <- quantile(emp.cdf.fs, qs)
emp.cdf.tn_conserv <- ecdf(pvals_heller_conserv)
empqs.tn_conserv <- quantile(emp.cdf.tn_conserv, qs)
emp.cdf.tn_oracle <- ecdf(pvals_heller_oracle)
empqs.tn_oracle <- quantile(emp.cdf.tn_oracle, qs)
emp.cdf.tn_switch <- ecdf(pvals_heller_switch)
empqs.tn_switch <- quantile(emp.cdf.tn_switch, qs)

par(pty = "s")
plot(theorqs, empqs.fs, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
points(theorqs, empqs.tn_oracle, type = "l", lwd=3, col = "lightpink")
points(theorqs, empqs.tn, type = "l", lwd=3, col = "#009E73")
points(theorqs, empqs.tn_conserv, type = "l", lwd=3, col = "#0072B2")
points(theorqs, empqs.tn_switch, type = "l", lwd=3, col = "orange")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "lightpink","#009E73", "#0072B2", "orange"), c("Our proposal",expression("Heller's proposal (" * sigma^2 * ")"), expression("Heller's proposal (" * hat(sigma)^2 * ")"), expression("Heller's proposal (" * tilde(sigma)^2 * ")"), expression("Heller's proposal (" * breve(sigma)^2 * ")")))


######################### final plot ########################

pdf("heller_figs.pdf", width = 8, height = 3)

layout(matrix(1:3, nrow = 1),
       widths = c(1, 1, 0.5))

# --- plot 1 ---
par(mar = c(4, 4, 1.5, 1), pty = "s")
plot(pvals_F, pvals_chisq_oracle,
     xlab = "F-statistic p-values",
     ylab = "Chi-squared p-values",
     xlim = c(0,1), ylim = c(0,1),
     pch = 16, cex = 0.5, col = "lightpink")
points(pvals_F, pvals_chisq_conserv,
       pch = 16, cex = 0.5, col = "#0072B2")
points(pvals_F, pvals_chisq_plugin,
       pch = 16, cex = 0.5, col = "#009E73")
abline(0, 1, col = "red", lwd =2)
abline(h = 0.05, col = "red", lty = 2, lwd =2)
abline(v = 0.05, col = "red", lty = 2, lwd =2)

# --- plot 2 ---
par(mar = c(4, 4, 1.5, 1), pty = "s")
plot(theorqs, empqs.fs, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
points(theorqs, empqs.tn_oracle, type = "l", lwd=3, col = "lightpink")
points(theorqs, empqs.tn, type = "l", lwd=3, col = "#009E73")
points(theorqs, empqs.tn_conserv, type = "l", lwd=3, col = "#0072B2")
segments(0, 0, 1, 1, col = "red")

# --- legend panel ---
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", pch=15, col=c("black", "lightpink","#009E73", "#0072B2"), c("(i) Our proposal",expression("(ii) Heller's proposal (" * sigma^2 * ")"), expression("(iii) Heller's proposal (" * hat(sigma)^2 * ")"), expression("(iv) Heller's proposal (" * tilde(sigma)^2 * ")")))

dev.off()

