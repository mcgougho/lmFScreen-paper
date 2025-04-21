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

############################### Global variables ###################################

set.seed(1)
n <- 100
p <- 10
sigma <- 1
beta <- rep(0,p)
Zs <- rchisq(10000, n-p-1)
As <- rchisq(10000, p)


############################## Type 1 error figure ################################

niter <- 10000
alpha <- 0.05
F.quantile <- qf(1 - alpha, p, (n - p - 1))
cc <- p / (n - p - 1) * F.quantile
scaling <- mean(Zs[As >= cc*Zs])
p_naive <- psel_oracle <- psel_DB <- psel_plugin <- rep(NA, niter)
for(iter in 1:niter){
  repeat {
    X <- matrix(rnorm(n * p), ncol = p)
    y <- X %*% beta + rnorm(n, mean=0, sd=sigma)
    Xy_centered <- lmFScreen:::get_Xy_centered(X,y)
    X <- Xy_centered$X
    y <- Xy_centered$y
    U <- svd(X)$u
    yPy <- sum((t(U)%*%y)^2)
    rss <- sum(y^2) - yPy
    F_statistic <- (n - p - 1) / p * (yPy / rss)
    if (F_statistic >= F.quantile) {
      break
    }
  }
  mod <- lm(y~X+0)
  p_naive[iter] <- summary(mod)$coef[1,4]
  sigma_sq <- sigma^2
  psel <- lmFScreen:::get_pselb(X=X,y=y,sigma_sq=sigma_sq,yPy=yPy,rss=rss,seed=iter,alpha_ov=0.05,min_select=100)
  psel_oracle[iter] <- psel(0)
  sigma_sq <- rss/scaling
  pselb <- lmFScreen:::get_pselb(X=X,y=y,sigma_sq=sigma_sq,yPy=yPy,rss=rss,seed=iter,alpha_ov=0.05,min_select=100)
  psel_DB[iter] <- pselb(0)
  sigma_sq <- rss/(n-p-1)
  pselb <- lmFScreen:::get_pselb(X=X,y=y,sigma_sq=sigma_sq,yPy=yPy,rss=rss,seed=iter,alpha_ov=0.05,min_select=100)
  psel_plugin[iter] <- pselb(0)
}

qs <- seq(from = 0, to = 1, length.out = 1000)
theorqs <- qunif(qs)
emp.cdf.naive <- ecdf(p_naive)
empqs.naive <- quantile(emp.cdf.naive, qs)
emp.cdf.oracle <- ecdf(psel_oracle)
empqs.oracle <- quantile(emp.cdf.oracle, qs)
emp.cdf.DB <- ecdf(psel_DB)
empqs.DB <- quantile(emp.cdf.DB, qs)
emp.cdf.plugin <- ecdf(psel_plugin)
empqs.plugin <- quantile(emp.cdf.plugin, qs)

par(pty = "s")
plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs.oracle, lwd=3, col = "purple2")
lines(theorqs, empqs.DB, lwd=3, col = "blue")
lines(theorqs, empqs.plugin,lwd=3, col = "lightblue")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2", "blue" , "lightblue"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E]), expression(p[H[0]^M*"|"*E]^{tilde(sigma)^2}), expression(p[H[0]^M*"|"*E]^{hat(sigma)^2})))



############################## Power figures ################################

pr_reject_split <- function(beta1,n,p,alpha_ov,prop_train=0.5,sigma=1){
  if (beta1 == 0){nreps = 50000}
  else{nreps <- round(5000/beta1)}
  computeFstat <- rep(FALSE, nreps)
  pvals.split  <- rep(NA, nreps)
  beta[1] <- beta1
  for(iter in 1:nreps){
    X <- matrix(rnorm(n*p),ncol=p) # make design matrix
    y <- X%*%beta+rnorm(n)*sigma # construct y
    smp_size <- floor(prop_train * n) # make sample size
    train_ind <- sample(seq_len(n), size = smp_size) # select train indices
    X_train <- X[train_ind,] # make X_train
    X_train <- scale(X_train, T, F)
    X_test <- X[-train_ind,] # make X_test
    X_test <- scale(X_test, T, F)
    y_train <- y[train_ind] # make y_train
    y_train <- scale(y_train, T, F)
    y_test <- y[-train_ind] # make y_test
    y_test <- scale(y_test, T, F)
    mod_train <- lm(y_train~X_train+0)
    overall.F <- summary(mod_train)$fstat[1]
    if(overall.F >= qf(1-alpha_ov, p, smp_size-p-1)){ # if F-stat is significant
      computeFstat[iter] <- TRUE
      mod_test <- lm(y_test~X_test+0)
      psplit <- summary(mod_test)$coef[1,4]
      pvals.split[iter] <- psplit
    }
  }
  pvals.split.noNA <- pvals.split[computeFstat]
  pr_reject_test <- mean(pvals.split.noNA <= 0.05)
  pr_reject_train <- mean(computeFstat)
  return(list(pr_reject_test = pr_reject_test, pr_reject_train = pr_reject_train))
}


pr_reject_sel <- function(beta1,n,p,alpha_ov,sigma=1,oracle=TRUE){
  if (beta1 == 0){nreps = 50000}
  else{nreps <- round(5000/beta1)}
  beta[1] <- beta1
  F.quantile <- qf(1 - alpha_ov, p, (n - p - 1))
  cc <- p / (n - p - 1) * F.quantile
  scaling <- mean(Zs[As >= cc*Zs])
  pvals.sel <- rep(NA, nreps)
  compute_F <- rep(FALSE, nreps)
  sigma_est <- rep(NA, nreps)
  for(i in 1:nreps){
    X <- matrix(rnorm(n * p), ncol = p)
    y  <- X %*%beta+rnorm(n)*sigma
    Xy_info <- lmFScreen:::get_Xy_centered(X,y) # from lmFScreen package
    X <- Xy_info$X
    y <- Xy_info$y
    U <- svd(X)$u
    yPy <- sum((t(U)%*%y)^2)
    rss <- sum(y^2) - yPy
    if (oracle == TRUE){
      sigma_sq <- sigma^2
    }
    else{
      sigma_sq <- rss/scaling
    }
    F_overall <- (n-p-1)/p * yPy/rss
    if(F_overall >= F.quantile){
      compute_F[i] = TRUE
      sigma_est[i] <- sigma_sq
      pselb <- lmFScreen:::get_pselb(X=X,y=y,sigma_sq=sigma_sq,yPy=yPy,rss=rss,
                                     alpha_ov=alpha_ov,min_select=50,seed=i) # from lmFScreen package
      pvals.sel[i] <- pselb(0)
    }
  }
  pvals.sel.noNA <- pvals.sel[compute_F]
  pr_reject_sel <- mean(pvals.sel.noNA <= 0.05)
  pr_reject_overall <- mean(compute_F)
  return(list(pr_reject_sel = pr_reject_sel, pr_reject_overall = pr_reject_overall))
}


compute_psel <- function(alpha_ov, beta1_val, oracle = TRUE) {
  pr_reject <- pr_reject_sel(beta1 = beta1_val, n = 100, p = 10, alpha_ov = alpha_ov, oracle = oracle)
  return(c(pr_reject[1], pr_reject[2]))
}

compute_psplit <- function(alpha_ov, beta1_val) {
  pr_reject <- pr_reject_split(beta1 = beta1_val, n = 100, p = 10, alpha_ov = alpha_ov)
  return(c(pr_reject[1], pr_reject[2]))
}

run_simulation <- function(beta1_val) {
  print(beta1_val)
  alpha_levels <- c(0.05, 0.1, 0.5)
  results_psel <- sapply(alpha_levels, function(a) compute_psel(a, beta1_val, oracle = TRUE))
  results_psel_hat <- sapply(alpha_levels, function(a) compute_psel(a, beta1_val, oracle = FALSE))
  results_psplit <- sapply(alpha_levels, function(a) compute_psplit(a, beta1_val))
  return(list(
    pr_reject_overall = results_psel[2, ],
    pr_reject_sel = results_psel[1, ],
    pr_reject_sel_hat = results_psel_hat[1, ],
    pr_reject_split = results_psplit[1, ],
    pr_reject_train = results_psplit[2, ]
  ))
}

# Run simulation for each beta1 value
beta1 <- seq(0, 0.8, 0.1)
results <- lapply(beta1, run_simulation)


# Extract vectors from results into individual named vectors for plotting
pr_reject_overall_05 <- sapply(results, function(res) res$pr_reject_overall[[1]])
pr_reject_overall_1  <- sapply(results, function(res) res$pr_reject_overall[[2]])
pr_reject_overall_5  <- sapply(results, function(res) res$pr_reject_overall[[3]])

pr_reject_sel_05     <- sapply(results, function(res) res$pr_reject_sel[[1]])
pr_reject_sel_1      <- sapply(results, function(res) res$pr_reject_sel[[2]])
pr_reject_sel_5      <- sapply(results, function(res) res$pr_reject_sel[[3]])

pr_reject_sel_hat_05 <- sapply(results, function(res) res$pr_reject_sel_hat[[1]])
pr_reject_sel_hat_1  <- sapply(results, function(res) res$pr_reject_sel_hat[[2]])
pr_reject_sel_hat_5  <- sapply(results, function(res) res$pr_reject_sel_hat[[3]])

pr_reject_split_05   <- sapply(results, function(res) res$pr_reject_split[[1]])
pr_reject_split_1    <- sapply(results, function(res) res$pr_reject_split[[2]])
pr_reject_split_5    <- sapply(results, function(res) res$pr_reject_split[[3]])

pr_reject_train_05   <- sapply(results, function(res) res$pr_reject_train[[1]])
pr_reject_train_1    <- sapply(results, function(res) res$pr_reject_train[[2]])
pr_reject_train_5    <- sapply(results, function(res) res$pr_reject_train[[3]])




par(mfrow = c(1, 2))

plot(beta1, pr_reject_train_05, ylim=c(0, 1), pch=2,type = "b",xlab=expression("Values of " ~ beta[1]), ylab="Probability of rejecting overall F-test", lty = 3,col="#E69F00")
points(beta1, pr_reject_train_1, ylim=c(0, 1), pch=2,type = "b", lty = 2,col="#E69F00")
points(beta1, pr_reject_train_5, ylim=c(0, 1), pch=2,type = "b", col="#E69F00")
points(beta1, pr_reject_overall_05, ylim=c(0, 1), pch = 4, type = "b",lty = 3, col="#009E73")
points(beta1, pr_reject_overall_1, ylim=c(0, 1), pch = 4, type = "b",lty = 2, col="#009E73")
points(beta1, pr_reject_overall_5, ylim=c(0, 1), pch = 4, type = "b", col="#009E73")
legend("bottomright", pch = c(4, 2, NA, NA,NA), col=c("#009E73","#E69F00","black", "black", "black"), c(expression(F[H[0]^{1:p}]), bquote(F[H[0]^{1:p}]^{scriptscriptstyle("train")}),expression(alpha[0] == 0.05), expression(alpha[0] == 0.1), expression(alpha[0] == 0.5)), lty = c(NA, NA, 3,2,1))

plot(beta1, pr_reject_split_05, ylim=c(0, 1),pch=2, type = "b",xlab=expression("Values of " ~ beta[1]), ylab=expression("Conditional probability of rejecting " ~ H[0] * ": " ~ beta[1] == 0), lty = 3,col="#E69F00")
abline(h=0.05, col="darkgrey", lty=2)
points(beta1, pr_reject_split_1, ylim=c(0, 1), pch=2,type = "b", lty = 2, col="#E69F00")
points(beta1, pr_reject_split_5, ylim=c(0, 1), pch=2,type = "b", col="#E69F00")
points(beta1, pr_reject_sel_hat_05, ylim=c(0, 1), type = "b", lty = 3,col="#0072B2")
points(beta1, pr_reject_sel_hat_1, ylim=c(0, 1), type = "b", lty = 2, col="#0072B2")
points(beta1, pr_reject_sel_hat_5, ylim=c(0, 1), type = "b", col="#0072B2")
points(beta1, pr_reject_sel_05, ylim=c(0, 1),pch = 4, type = "b", lty = 3, col= "#009E73")
points(beta1, pr_reject_sel_1, ylim=c(0, 1), pch = 4, type = "b", lty = 2, col= "#009E73")
points(beta1, pr_reject_sel_5, ylim=c(0, 1), pch = 4, type = "b", col= "#009E73")
legend("bottomright", pch = c(4, 1, 2, NA, NA,NA), col=c( "#009E73",  "#0072B2","#E69F00","black", "black", "black"), c(expression(p[H[0]^M*"|"*E]), expression(p[H[0]^M*"|"*E]^{tilde(sigma)^2}), bquote(p[H[0]^M]^{scriptscriptstyle("test")}), expression(alpha[0] == 0.05), expression(alpha[0] == 0.1), expression(alpha[0] == 0.5)), lty = c(NA, NA, NA, 3,2,1))


################################ final plot #########################################

pdf("power_t1error_global_null.pdf", width = 8, height = 3)

par(mfrow = c(1, 3), pty = "s")

plot(theorqs, empqs.naive, type = "l", lwd=3, xlab = "Uniform Theoretical Quantiles", ylab = "Empirical Quantiles", xlim = c(0, 1), ylim = c(0, 1), col = "black")
lines(theorqs, empqs.oracle, lwd=3, col = "purple2")
lines(theorqs, empqs.DB, lwd=3, col = "blue")
lines(theorqs, empqs.plugin,lwd=3, col = "lightblue")
segments(0, 0, 1, 1, col = "red")
legend("topleft", pch=15, col=c("black", "purple2", "blue" , "lightblue"), c(expression(p[H[0]^M]), expression(p[H[0]^M*"|"*E]), expression(p[H[0]^M*"|"*E]^{tilde(sigma)^2}), expression(p[H[0]^M*"|"*E]^{hat(sigma)^2})))


plot(beta1, pr_reject_train_05, ylim=c(0, 1), pch=2,type = "b",xlab=expression("Values of " ~ beta[1]), ylab="Probability of rejecting overall F-test", lty = 3,col="#E69F00")
points(beta1, pr_reject_train_1, ylim=c(0, 1), pch=2,type = "b", lty = 2,col="#E69F00")
points(beta1, pr_reject_train_5, ylim=c(0, 1), pch=2,type = "b", col="#E69F00")
points(beta1, pr_reject_overall_05, ylim=c(0, 1), pch = 4, type = "b",lty = 3, col="#009E73")
points(beta1, pr_reject_overall_1, ylim=c(0, 1), pch = 4, type = "b",lty = 2, col="#009E73")
points(beta1, pr_reject_overall_5, ylim=c(0, 1), pch = 4, type = "b", col="#009E73")
legend("bottomright", pch = c(4, 2, NA, NA,NA), col=c("#009E73","#E69F00","black", "black", "black"), c(expression(F[H[0]^{1:p}]), bquote(F[H[0]^{1:p}]^{scriptscriptstyle("train")}),expression(alpha[0] == 0.05), expression(alpha[0] == 0.1), expression(alpha[0] == 0.5)), lty = c(NA, NA, 3,2,1))

plot(beta1, pr_reject_split_05, ylim=c(0, 1),pch=2, type = "b",xlab=expression("Values of " ~ beta[1]), ylab=expression("Conditional probability of rejecting " ~ H[0] * ": " ~ beta[1] == 0), lty = 3,col="#E69F00")
abline(h=0.05, col="darkgrey", lty=2)
points(beta1, pr_reject_split_1, ylim=c(0, 1), pch=2,type = "b", lty = 2, col="#E69F00")
points(beta1, pr_reject_split_5, ylim=c(0, 1), pch=2,type = "b", col="#E69F00")
points(beta1, pr_reject_sel_hat_05, ylim=c(0, 1), type = "b", lty = 3,col="#0072B2")
points(beta1, pr_reject_sel_hat_1, ylim=c(0, 1), type = "b", lty = 2, col="#0072B2")
points(beta1, pr_reject_sel_hat_5, ylim=c(0, 1), type = "b", col="#0072B2")
points(beta1, pr_reject_sel_05, ylim=c(0, 1),pch = 4, type = "b", lty = 3, col= "#009E73")
points(beta1, pr_reject_sel_1, ylim=c(0, 1), pch = 4, type = "b", lty = 2, col= "#009E73")
points(beta1, pr_reject_sel_5, ylim=c(0, 1), pch = 4, type = "b", col= "#009E73")
legend("bottomright", pch = c(4, 1, 2, NA, NA,NA), col=c( "#009E73",  "#0072B2","#E69F00","black", "black", "black"), c(expression(p[H[0]^M*"|"*E]), expression(p[H[0]^M*"|"*E]^{tilde(sigma)^2}), bquote(p[H[0]^M]^{scriptscriptstyle("test")}), expression(alpha[0] == 0.05), expression(alpha[0] == 0.1), expression(alpha[0] == 0.5)), lty = c(NA, NA, NA, 3,2,1),bg = "white")


dev.off()





