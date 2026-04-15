library(lmFScreen)
library(KScorrect)
library(VGAM)


##################### heteroskedasticity #####################


simulate_anova_het <- function(nA, nB, nC,
                               mu = c(0, 0, 0),
                               sigma = c(1, 2, 1.5)) {
  stopifnot(length(mu) == 3, length(sigma) == 3)
  
  group <- factor(rep(c("A", "B", "C"), c(nA, nB, nC)),
                  levels = c("A", "B", "C"))
  
  y <- c(
    mu[1] + rnorm(nA, sd = sigma[1]),
    mu[2] + rnorm(nB, sd = sigma[2]),
    mu[3] + rnorm(nC, sd = sigma[3])
  )
  
  X_raw <- cbind(
    I_B = as.numeric(group == "B"),
    I_C = as.numeric(group == "C")
  )
  
  centered <- lmFScreen:::get_Xy_centered(X_raw, y) # this makes it centered wrt group 1
  
  list(
    X = centered$X,
    y = centered$y,
    group = group,
    y_raw = y,
    X_raw = X_raw
  )
}


nA <- 20
nB <- 20
nC <- 20
n_total <- nA + nB + nC
p <- 2
alpha <- 0.05
niter <- 5000

# need mu_A = mu_B, but mu_C can be whatever you want
mu <- c(0, 0, 0)
sigma_grp <- c(1, 2, 1.5)   # heteroskedasticity by group

F.quantile <- qf(1 - alpha, df1 = p, df2 = n_total - p - 1)

p_overall <- c()
p_naive <- numeric(niter)
psel <- numeric(niter)
stats <- numeric(niter)

for (iter in 1:niter) {
  repeat {
    dat <- simulate_anova_het(
      nA = nA, nB = nB, nC = nC,
      mu = mu,
      sigma = sigma_grp
    )
    X <- dat$X
    y <- dat$y
    U <- svd(X)$u
    yPy <- sum((crossprod(U, y))^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n_total - p - 1) / p * (yPy / rss)
    p_overall <- c(p_overall, 1 - pf(F_statistic, df1 = p, df2 = n_total - p - 1))
    if (F_statistic >= F.quantile) {
      break
    }
  }
  
  # Naive test of coef 1 = 0, i.e. mu_B = mu_A
  mod <- lm(y ~ X + 0)
  p_naive[iter] <- summary(mod)$coef[1, 4]
  test_stat <- summary(mod)$coef[1, 3]
  stats[iter] <- test_stat^2
  
  # Selective p-value for coef 1 = 0
  pselb_fun <- lmFScreen:::get_pselb(X = X, y = y, alpha_ov = alpha, min_select = 10)
  psel[iter] <- pselb_fun(0)
}


qs <- seq(from = 0, to = 1, length.out = 5000)
theorqs <- qunif(qs)
emp.cdf.overall <- ecdf(p_overall)
empqs.overall <- quantile(emp.cdf.overall, qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf <- ecdf(psel)
empqs <- quantile(emp.cdf, qs)

pdf("heterosked.pdf", width = 8, height = 4)

par(pty = "s", mfrow = c(1, 2))

plot(theorqs, empqs.overall, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
abline(a = 0.05, b = 0, col = "red", lwd = 1, lty = 2)
segments(0, 0, 1, 1, col = "red")


plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))

dev.off()




########################### non normality ###########################


simulate_data_nonnormal <- function(n, p, beta, sigma, error_dist = c("t", "chisq", "exp", "laplace", "beta", "unif", "pois", "twopt")) {
  error_dist <- match.arg(error_dist)
  
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  
  if (error_dist == "t") {
    df <- 2.5
    # Var(t_df) = df / (df - 2), so rescale to variance sigma^2
    eps <- rt(n, df = df) * sigma / sqrt(df / (df - 2))
  } else if (error_dist == "chisq") {
    # centered + scaled chi-square with variance sigma^2
    eps <- (rchisq(n, df = 1) - 1) * sigma / sqrt(2)
  } else if (error_dist == "exp") {
    # centered + scaled exponential with variance sigma^2
    eps <- (rexp(n) - 1) * sigma
  } else if (error_dist == "laplace") {
    # Laplace distribution with mean 0 and variance sigma^2
    b <- sigma / sqrt(2)  # scale parameter for Laplace
    eps <- rlaplace(n, location = 0, scale = b)
  } else if(error_dist == "beta") {
    # centered + scaled beta distribution with variance sigma^2
    a <- 2
    b <- 1
    eps <- (rbeta(n, shape1 = a, shape2 = b) - a / (a + b)) * sigma / sqrt(a * b / ((a + b)^2 * (a + b + 1)))
  } else if (error_dist == "unif") {
    # centered + scaled uniform distribution with variance sigma^2
    eps <- (runif(n) - 0.5) * sigma * sqrt(12)
  } else if (error_dist == "pois"){
    lambda <- 0.2
    eps <- (rpois(n, lambda) - lambda) * sigma / sqrt(lambda)
  } else if (error_dist == "twopt") {
    eps <- sample(c(-sigma, 20*sigma), size = n, replace = TRUE)
  } 
  y <- X %*% beta + eps
  list(X = X, y = as.numeric(y))
}

# laplace 

niter <- 5000
n <- 30
p <- 2
sigma <- 1
beta <- rep(0,p)
alpha <- 0.05
F.quantile <- qf(1 - alpha, df1 = p, df2 = n - p - 1)

p_naive <- numeric(niter)
psel <- numeric(niter)
samples <- c()

for (iter in 1:niter) {
  repeat {
    dat <- simulate_data_nonnormal(n, p, beta, sigma, error_dist = "laplace")
    
    Xy_centered <- lmFScreen:::get_Xy_centered(dat$X, dat$y)
    X <- Xy_centered$X
    y <- Xy_centered$y
    
    U <- svd(X)$u
    yPy <- sum((crossprod(U, y))^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n - p - 1) / p * (yPy / rss)
    
    if (F_statistic >= F.quantile) break
  }
  
  samples <- c(samples, y)
  mod <- lm(y ~ X + 0)
  p_naive[iter] <- summary(mod)$coef[1, 4]
  
  pselb_fun <- lmFScreen:::get_pselb(X = X, y = y, alpha_ov = alpha, min_select = 100)
  psel[iter] <- pselb_fun(0)
}

qs <- seq(from = 0, to = 1, length.out = 5000)
theorqs <- qunif(qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf <- ecdf(psel)
empqs <- quantile(emp.cdf, qs)

par(pty = "s")

plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))



pdf("non_normal_laplace.pdf", width = 5, height = 4)
par(pty = "s")

plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))

dev.off()




# exponential 

niter <- 5000
n <- 30
p <- 2
sigma <- 1
beta <- rep(0,p)
alpha <- 0.05
F.quantile <- qf(1 - alpha, df1 = p, df2 = n - p - 1)

p_naive <- numeric(niter)
psel <- numeric(niter)
samples <- c()

for (iter in 1:niter) {
  repeat {
    dat <- simulate_data_nonnormal(n, p, beta, sigma, error_dist = "exp")
    
    Xy_centered <- lmFScreen:::get_Xy_centered(dat$X, dat$y)
    X <- Xy_centered$X
    y <- Xy_centered$y
    
    U <- svd(X)$u
    yPy <- sum((crossprod(U, y))^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n - p - 1) / p * (yPy / rss)
    
    if (F_statistic >= F.quantile) break
  }
  
  samples <- c(samples, y)
  mod <- lm(y ~ X + 0)
  p_naive[iter] <- summary(mod)$coef[1, 4]
  
  pselb_fun <- lmFScreen:::get_pselb(X = X, y = y, alpha_ov = alpha, min_select = 100)
  psel[iter] <- pselb_fun(0)
}

qs <- seq(from = 0, to = 1, length.out = 5000)
theorqs <- qunif(qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf <- ecdf(psel)
empqs <- quantile(emp.cdf, qs)

par(pty = "s")

plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))



pdf("non_normal_exp.pdf", width = 5, height = 4)
par(pty = "s")

plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))

dev.off()



# beta 


niter <- 5000
n <- 30
p <- 2
sigma <- 1
beta <- rep(0,p)

p_naive <- numeric(niter)
psel <- numeric(niter)
samples <- c()

for (iter in 1:niter) {
  repeat {
    dat <- simulate_data_nonnormal(n, p, beta, sigma, error_dist = "beta")
    
    Xy_centered <- lmFScreen:::get_Xy_centered(dat$X, dat$y)
    X <- Xy_centered$X
    y <- Xy_centered$y
    
    U <- svd(X)$u
    yPy <- sum((crossprod(U, y))^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n - p - 1) / p * (yPy / rss)
    
    if (F_statistic >= F.quantile) break
  }
  
  samples <- c(samples, y)
  mod <- lm(y ~ X + 0)
  p_naive[iter] <- summary(mod)$coef[1, 4]
  
  pselb_fun <- lmFScreen:::get_pselb(X = X, y = y, alpha_ov = alpha, min_select = 100)
  psel[iter] <- pselb_fun(0)
}

qs <- seq(from = 0, to = 1, length.out = 5000)
theorqs <- qunif(qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf <- ecdf(psel)
empqs <- quantile(emp.cdf, qs)

par(pty = "s")


plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))



pdf("non_normal_beta.pdf", width = 5, height = 4)
par(pty = "s")


plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))

dev.off()





# rademacher

niter <- 5000
n <- 30
p <- 2
sigma <- 1
beta <- rep(0,p)

p_naive <- numeric(niter)
psel <- numeric(niter)
samples <- c()

for (iter in 1:niter) {
  repeat {
    dat <- simulate_data_nonnormal(n, p, beta, sigma, error_dist = "twopt")
    
    Xy_centered <- lmFScreen:::get_Xy_centered(dat$X, dat$y)
    X <- Xy_centered$X
    y <- Xy_centered$y
    
    U <- svd(X)$u
    yPy <- sum((crossprod(U, y))^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n - p - 1) / p * (yPy / rss)
    
    if (F_statistic >= F.quantile) break
  }
  
  samples <- c(samples, y)
  mod <- lm(y ~ X + 0)
  p_naive[iter] <- summary(mod)$coef[1, 4]
  
  pselb_fun <- lmFScreen:::get_pselb(X = X, y = y, alpha_ov = alpha, min_select = 100)
  psel[iter] <- pselb_fun(0)
}

qs <- seq(from = 0, to = 1, length.out = 5000)
theorqs <- qunif(qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf <- ecdf(psel)
empqs <- quantile(emp.cdf, qs)

par(pty = "s")


plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))


pdf("non_normal_twopt.pdf", width = 5, height = 4)
par(pty = "s")


plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))

dev.off()



# use anova design matrix 



simulate_data_nonnormal_anova <- function(n, p, beta, sigma, error_dist = c("t", "chisq", "exp", "laplace", "beta", "unif", "pois", "twopt")) {
  error_dist <- match.arg(error_dist)
  
  nA <- nB <- nC <- n / 3
  group <- factor(rep(c("A", "B", "C"), c(nA, nB, nC)),
                  levels = c("A", "B", "C"))
  X <- cbind(
    I_B = as.numeric(group == "B"),
    I_C = as.numeric(group == "C")
  )
  
  
  if (error_dist == "t") {
    df <- 2.5
    # Var(t_df) = df / (df - 2), so rescale to variance sigma^2
    eps <- rt(n, df = df) * sigma / sqrt(df / (df - 2))
  } else if (error_dist == "chisq") {
    # centered + scaled chi-square with variance sigma^2
    eps <- (rchisq(n, df = 1) - 1) * sigma / sqrt(2)
  } else if (error_dist == "exp") {
    # centered + scaled exponential with variance sigma^2
    eps <- (rexp(n) - 1) * sigma
  } else if (error_dist == "laplace") {
    # Laplace distribution with mean 0 and variance sigma^2
    b <- sigma / sqrt(2)  # scale parameter for Laplace
    eps <- rlaplace(n, location = 0, scale = b)
  } else if(error_dist == "beta") {
    # centered + scaled beta distribution with variance sigma^2
    a <- 2
    b <- 1
    eps <- (rbeta(n, shape1 = a, shape2 = b) - a / (a + b)) * sigma / sqrt(a * b / ((a + b)^2 * (a + b + 1)))
  } else if (error_dist == "unif") {
    # centered + scaled uniform distribution with variance sigma^2
    eps <- (runif(n) - 0.5) * sigma * sqrt(12)
  } else if (error_dist == "pois"){
    lambda <- 0.2
    eps <- (rpois(n, lambda) - lambda) * sigma / sqrt(lambda)
  } else if (error_dist == "twopt") {
    eps <- sample(c(-sigma, 20*sigma), size = n, replace = TRUE)
  } 
  y <- X %*% beta + eps
  list(X = X, y = as.numeric(y))
}



# laplace with anova


niter <- 5000
n <- 30
p <- 2
sigma <- 1
beta <- rep(0,p)
alpha <- 0.05
F.quantile <- qf(1 - alpha, df1 = p, df2 = n - p - 1)

p_naive <- numeric(niter)
psel <- numeric(niter)
samples <- c()

for (iter in 1:niter) {
  repeat {
    dat <- simulate_data_nonnormal_anova(n, p, beta, sigma, error_dist = "laplace")
    
    Xy_centered <- lmFScreen:::get_Xy_centered(dat$X, dat$y)
    X <- Xy_centered$X
    y <- Xy_centered$y
    
    U <- svd(X)$u
    yPy <- sum((crossprod(U, y))^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n - p - 1) / p * (yPy / rss)
    
    if (F_statistic >= F.quantile) break
  }
  
  samples <- c(samples, y)
  mod <- lm(y ~ X + 0)
  p_naive[iter] <- summary(mod)$coef[1, 4]
  
  pselb_fun <- lmFScreen:::get_pselb(X = X, y = y, alpha_ov = alpha, min_select = 100)
  psel[iter] <- pselb_fun(0)
}

qs <- seq(from = 0, to = 1, length.out = 5000)
theorqs <- qunif(qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf <- ecdf(psel)
empqs <- quantile(emp.cdf, qs)

par(pty = "s")

plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))



pdf("non_normal_laplace_anova.pdf", width = 5, height = 4)
par(pty = "s")

plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))

dev.off()




# exponential with anova

niter <- 5000
n <- 30
p <- 2
sigma <- 1
beta <- rep(0,p)
alpha <- 0.05
F.quantile <- qf(1 - alpha, df1 = p, df2 = n - p - 1)

p_naive <- numeric(niter)
psel <- numeric(niter)
samples <- c()

for (iter in 1:niter) {
  repeat {
    dat <- simulate_data_nonnormal_anova(n, p, beta, sigma, error_dist = "exp")
    
    Xy_centered <- lmFScreen:::get_Xy_centered(dat$X, dat$y)
    X <- Xy_centered$X
    y <- Xy_centered$y
    
    U <- svd(X)$u
    yPy <- sum((crossprod(U, y))^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n - p - 1) / p * (yPy / rss)
    
    if (F_statistic >= F.quantile) break
  }
  
  samples <- c(samples, y)
  mod <- lm(y ~ X + 0)
  p_naive[iter] <- summary(mod)$coef[1, 4]
  
  pselb_fun <- lmFScreen:::get_pselb(X = X, y = y, alpha_ov = alpha, min_select = 100)
  psel[iter] <- pselb_fun(0)
}

qs <- seq(from = 0, to = 1, length.out = 5000)
theorqs <- qunif(qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf <- ecdf(psel)
empqs <- quantile(emp.cdf, qs)

par(pty = "s")

plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))



pdf("non_normal_exp_anova.pdf", width = 5, height = 4)
par(pty = "s")

plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))

dev.off()





# beta with anova

niter <- 5000
n <- 30
p <- 2
sigma <- 1
beta <- rep(0,p)

p_naive <- numeric(niter)
psel <- numeric(niter)
samples <- c()

for (iter in 1:niter) {
  repeat {
    dat <- simulate_data_nonnormal_anova(n, p, beta, sigma, error_dist = "beta")
    
    Xy_centered <- lmFScreen:::get_Xy_centered(dat$X, dat$y)
    X <- Xy_centered$X
    y <- Xy_centered$y
    
    U <- svd(X)$u
    yPy <- sum((crossprod(U, y))^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n - p - 1) / p * (yPy / rss)
    
    if (F_statistic >= F.quantile) break
  }
  
  samples <- c(samples, y)
  mod <- lm(y ~ X + 0)
  p_naive[iter] <- summary(mod)$coef[1, 4]
  
  pselb_fun <- lmFScreen:::get_pselb(X = X, y = y, alpha_ov = alpha, min_select = 100)
  psel[iter] <- pselb_fun(0)
}

qs <- seq(from = 0, to = 1, length.out = 5000)
theorqs <- qunif(qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf <- ecdf(psel)
empqs <- quantile(emp.cdf, qs)

par(pty = "s")


plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))



pdf("non_normal_beta_ANOVA.pdf", width = 5, height = 4)
par(pty = "s")


plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))

dev.off()



# rademacher with anova 

niter <- 5000
n <- 30
p <- 2
sigma <- 1
beta <- rep(0,p)

p_naive <- numeric(niter)
psel <- numeric(niter)
samples <- c()

for (iter in 1:niter) {
  repeat {
    dat <- simulate_data_nonnormal_anova(n, p, beta, sigma, error_dist = "twopt")
    
    Xy_centered <- lmFScreen:::get_Xy_centered(dat$X, dat$y)
    X <- Xy_centered$X
    y <- Xy_centered$y
    
    U <- svd(X)$u
    yPy <- sum((crossprod(U, y))^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n - p - 1) / p * (yPy / rss)
    
    if (F_statistic >= F.quantile) break
  }
  
  samples <- c(samples, y)
  mod <- lm(y ~ X + 0)
  p_naive[iter] <- summary(mod)$coef[1, 4]
  
  pselb_fun <- lmFScreen:::get_pselb(X = X, y = y, alpha_ov = alpha, min_select = 100)
  psel[iter] <- pselb_fun(0)
}

qs <- seq(from = 0, to = 1, length.out = 5000)
theorqs <- qunif(qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf <- ecdf(psel)
empqs <- quantile(emp.cdf, qs)

par(pty = "s")


plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))


pdf("non_normal_twopt_ANOVA.pdf", width = 5, height = 4)
par(pty = "s")


plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs, lwd=3, col = "purple2")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E])))

dev.off()



