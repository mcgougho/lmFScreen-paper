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

fisher_info_sel <- function(X, Y, beta1, beta_minus_1, n, p, sigma){
  beta <- as.matrix(c(beta1,beta_minus_1))
  U <- svd(X[,-1])$u
  yPXminus1y <- sum((t(U) %*% Y)^2)
  dtilde2 <- yPXminus1y
  e <- sigma^2 * qchisq(0.95, df=p) - dtilde2
  Xperm <- cbind(X[,-1], X[,1])
  Q <- qr.Q(qr(Xperm))     
  U2 <- Q[, p]  # orthonormal basis for (P_X - P_{X_{-1}})
  U2X1 <- t(U2)%*%X[,1]
  front_coef <- sum(U2X1^2)/sigma^2
  if (e <= 0) {
    return(as.numeric(front_coef))
  }
  a <- (-sqrt(e) - U2X1 * beta1)/sigma
  b <- (sqrt(e) - U2X1 * beta1)/sigma
  num1 <- (dnorm(a) - dnorm(b))^2
  denom1 <- (1 - pnorm(b) + pnorm(a))^2
  num2 <- (a * dnorm(a) - b * dnorm(b))
  denom2 <- (1 - pnorm(b) + pnorm(a))
  back_coef <- 1 - num1/denom1 - num2/denom2
  info=front_coef * back_coef
  return(as.numeric(info))
}


fisher_info_split <- function(X_test, beta1, beta_minus_1, n, p, sigma) {
  Xperm <- cbind(X_test[,-1], X_test[,1])
  Q <- qr.Q(qr(Xperm))     
  U2 <- Q[, p]  # orthonormal basis for (P_X - P_{X_{-1}})
  U2X1_test <- t(U2)%*%X_test[,1]
  info <- sum(U2X1_test^2) / sigma^2
  return(as.numeric(info))
}


power_overall_F <- function(X, n, p, beta1, beta_minus_1, sigma, alpha = 0.05) {
  beta <- as.matrix(c(beta1, beta_minus_1))
  df1 <- p
  df2 <- n-p-1
  mu <- as.vector(X %*% beta)
  lambda <- sum(mu^2) / sigma^2
  fcrit <- qf(1 - alpha, df1, df2)
  power <- 1 - pf(fcrit, df1, df2, ncp = lambda)
  return(as.numeric(power))
}


compute_fisher_info <- function(X, p, sigma, nrep, beta1, s_vals, rhos, alpha = 0.05, seed_splits = 1) {
  n <- nrow(X); stopifnot(ncol(X) == p)
  
  # fixed splits once
  set.seed(seed_splits)
  X_tr_list <- list(); X_te_list <- list()
  for (k in seq_along(rhos)) {
    n_train <- max(p + 1, floor(rhos[k] * n))
    idx     <- sample.int(n, n)
    idx_tr  <- idx[seq_len(n_train)]
    idx_te  <- idx[(n_train + 1):n]
    key <- as.character(rhos[k])
    X_tr_list[[key]] <- X[idx_tr, , drop = FALSE]
    X_te_list[[key]] <- X[idx_te, , drop = FALSE]
  }
  
  mean_fi_sel   <- numeric(length(s_vals))
  rej_prob      <- matrix(NA_real_, nrow = length(s_vals),
                          ncol = 1 + length(rhos),
                          dimnames = list(NULL, c("Selective", paste0("Split, rho=", rhos))))
  fi_split_mat  <- matrix(NA_real_, nrow = length(s_vals), ncol = length(rhos))
  
  for (j in seq_along(s_vals)) {
    s <- s_vals[j]
    beta_minus_1 <- rep(s, p - 1)
    
    # MC over Y for selective FI
    fi_vec <- numeric(nrep)
    for (i in seq_len(nrep)) {
      Y <- as.vector(X %*% c(beta1, beta_minus_1) + sigma * rnorm(n))
      U <- svd(X)$u
      yPy <- sum((t(U) %*% Y)^2)
      if (yPy > sigma^2 * qchisq(0.95, df = p)) {
        fi_vec[i] <- fisher_info_sel(X, Y, beta1, beta_minus_1, n, p, sigma)
      } else fi_vec[i] <- NA_real_
    }
    mean_fi_sel[j] <- mean(fi_vec, na.rm = TRUE)
    
    # deterministic overall power (full data)
    rej_prob[j, "Selective"] <- power_overall_F(
      X, n = n, p = p, beta1 = beta1, beta_minus_1 = beta_minus_1, sigma = sigma, alpha = alpha
    )
    
    # deterministic split FI and power (fixed splits)
    for (k in seq_along(rhos)) {
      key    <- as.character(rhos[k])
      X_tr   <- X_tr_list[[key]]
      X_te   <- X_te_list[[key]]
      n_train <- nrow(X_tr)
      
      fi_split_mat[j, k] <- fisher_info_split(X_te, beta1, beta_minus_1, n, p, sigma)
      rej_prob[j, paste0("Split, rho=", rhos[k])] <- power_overall_F(
        X_tr, n = n_train, p = p, beta1 = beta1, beta_minus_1 = beta_minus_1, sigma = sigma, alpha = alpha
      )
    }
  }
  
  list(rej_prob = rej_prob,
       fi_sel = mean_fi_sel,
       fi_split = fi_split_mat,
       beta1 = beta1)
}


draw_panel_fi <- function(res, s_vals, rhos, cols, ltys, col_idx) {
  mar_left  <- if (col_idx == 1) 2.0 else 1.2
  mar_right <- if (col_idx == 3) 0.2 else 0.3
  par(pty="s", mar=c(1.8, mar_left, 0.6, mar_right), mgp=c(1.0, 0.35, 0), xaxs="i", yaxs="i")
  yr <- range(c(res$fi_sel, res$fi_split), na.rm=TRUE)
  pad <- 0.12 * diff(yr); if (!is.finite(pad)) pad <- 1
  yr <- c(yr[1]-pad, yr[2]+pad)
  plot(s_vals, res$fi_sel, type="l", ylim=yr, axes=TRUE, ann=FALSE,
       col=cols[1], lwd=2)  # force green (Set2[1]) not black
  grid(col=rgb(0,0,0,0.1), lty=3)
  for (k in seq_along(rhos)) {
    lines(s_vals, res$fi_split[, k], col=cols[k+1], lwd=2, lty=ltys[k+1])
  }
  title(xlab="s", ylab="Leftover Fisher info", line=1.0, cex.lab=0.95)
}

draw_panel_prob <- function(res, s_vals, rhos, alpha, cols, ltys, col_idx) {
  # narrow margins: left more generous only for first column
  mar_left  <- if (col_idx == 1) 2.0 else 1.2
  mar_right <- if (col_idx == 3) 0.2 else 0.3
  par(pty="s", mar=c(1.8, mar_left, 0.6, mar_right), mgp=c(1.0, 0.35, 0), xaxs="i", yaxs="i")
  plot(NA, xlim=range(s_vals), ylim=c(0,1), axes=TRUE, ann=FALSE)
  grid(col=rgb(0,0,0,0.1), lty=3)
  lines(s_vals, res$rej_prob[, "Selective"], col=cols[1], lwd=2, lty=ltys[1])
  for (k in seq_along(rhos)) {
    lines(s_vals, res$rej_prob[, paste0("Split, rho=", rhos[k])],
          col=cols[k+1], lwd=2, lty=ltys[k+1])
  }
  title(xlab="s", ylab="Pr(reject overall)", line=1.0, cex.lab=0.95)
}

draw_panel_weighted <- function(res, s_vals, rhos, cols, ltys, col_idx) {
  mar_left  <- if (col_idx == 1) 2.0 else 1.2
  mar_right <- if (col_idx == 3) 0.2 else 0.3
  par(pty="s", mar=c(1.8, mar_left, 0.6, mar_right), mgp=c(1.0, 0.35, 0), xaxs="i", yaxs="i")
  expected_sel   <- res$rej_prob[, "Selective"] * res$fi_sel
  expected_split <- sapply(seq_along(rhos), function(k)
    res$rej_prob[, paste0("Split, rho=", rhos[k])] * res$fi_split[, k])
  yr <- range(c(expected_sel, expected_split), na.rm=TRUE)
  pad <- 0.12 * diff(yr); if (!is.finite(pad)) pad <- 1
  yr <- c(yr[1]-pad, yr[2]+pad)
  plot(NA, xlim=range(s_vals), ylim=yr, axes=TRUE, ann=FALSE)
  grid(col=rgb(0,0,0,0.1), lty=3)
  lines(s_vals, expected_sel, col=cols[1], lwd=2, lty=ltys[1])
  for (k in seq_along(rhos)) {
    lines(s_vals, expected_split[, k], col=cols[k+1], lwd=2, lty=ltys[k+1])
  }
  title(xlab="s", ylab="Screening-weighted info", line=1.0, cex.lab=0.95)
}



final_plot <- function(res_list, s_vals, rhos, alpha = 0.05) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Install RColorBrewer")
  cols <- RColorBrewer::brewer.pal(max(3, length(rhos) + 1), "Set2")
  ltys <- c(1, 3, 2, 4, 5, 6)[seq_len(length(rhos) + 1)]
  
  R <- length(res_list); C <- 3
  
  # regular panel ids
  panel_ids <- matrix(seq_len(R * C), nrow = R, ncol = C, byrow = TRUE)
  # legend row: one id repeated across all C columns (spans full width)
  legend_id <- R * C + 1L
  lay <- rbind(panel_ids, rep(legend_id, C))
  
  layout(lay,
         widths  = rep(1, C),
         heights = c(rep(1, R), 0.35))
  
  op <- par(mgp = c(1.6, 0.5, 0), oma = c(0,0,0,0))
  on.exit(par(op), add = TRUE)
  
  for (r in seq_len(R)) {
    res <- res_list[[r]]
    
    
    ## (1) Leftover Fisher information
    par(mar = c(2.6, 2.6, 1, 0.6), mgp = c(1.6, 0.5, 0), pty = "s")
    yr <- range(c(res$fi_sel, res$fi_split), na.rm = TRUE)
    pad <- 0.15 * diff(yr); if (!is.finite(pad)) pad <- 1
    yr <- c(yr[1] - pad, yr[2] + pad)
    plot(s_vals, res$fi_sel, type = "l", ylim = yr, col = cols[1], lwd = 2,
         xlab = "", ylab = "")
    grid(col = rgb(0,0,0,0.1), lty = 3)
    for (k in seq_along(rhos)) {
      lines(s_vals, res$fi_split[, k], col = cols[k+1], lwd = 2, lty = ltys[k+1])
    }
    title(xlab = "s", line = 1.6)

    
    ## (2) Overall rejection probability
    par(mar = c(2.6, 2.2, 1, 0.2), pty = "s")
    plot(NA, xlim = range(s_vals), ylim = c(0, 1), xlab = "", ylab = "")
    grid(col = rgb(0,0,0,0.1), lty = 3)
    lines(s_vals, res$rej_prob[, "Selective"], col = cols[1], lwd = 2, lty = ltys[1])
    for (k in seq_along(rhos)) {
      lines(s_vals, res$rej_prob[, paste0("Split, rho=", rhos[k])],
            col = cols[k+1], lwd = 2, lty = ltys[k+1])
    }
    title(xlab = "s", line = 1.6)
    
    
    ## (3) Screening-weighted information
    par(mar = c(2.6, 2.6, 1, 0.6), mgp = c(1.6, 0.5, 0), pty = "s")
    expected_sel   <- res$rej_prob[, "Selective"] * res$fi_sel
    expected_split <- sapply(seq_along(rhos), function(k)
      res$rej_prob[, paste0("Split, rho=", rhos[k])] * res$fi_split[, k])
    yr3 <- range(c(expected_sel, expected_split), na.rm = TRUE)
    pad3 <- 0.15 * diff(yr3); if (!is.finite(pad3)) pad3 <- 1
    yr3 <- c(yr3[1] - pad3, yr3[2] + pad3)
    plot(NA, xlim = range(s_vals), ylim = yr3, xlab = "", ylab = "")
    grid(col = rgb(0,0,0,0.1), lty = 3)
    lines(s_vals, expected_sel, col = cols[1], lwd = 2, lty = ltys[1])
    for (k in seq_along(rhos)) {
      lines(s_vals, expected_split[, k], col = cols[k+1], lwd = 2, lty = ltys[k+1])
    }
    title(xlab = "s", line = 1.6)
  } # end loop over rows r
  
  ## legend (occupies the spanning bottom row)
  par(mar = rep(0, 4)); plot.new()
  legend_expr <- do.call(expression, c(list(quote("Selective")), lapply(rhos, function(r) bquote(rho == .(r)))))
  legend("center", legend = legend_expr,
         col = cols[1:(length(rhos) + 1)], lwd = 2, lty = ltys,
         horiz = TRUE, bty = "n", xpd = NA, cex = 1.2)
}

final_plot_wrap <- function(X, p, sigma, nrep, beta1_vals, s_vals, rhos,
                            alpha = 0.05, seed_splits = 1) {
  res_list <- lapply(beta1_vals, function(b1) {
    compute_fisher_info(
      X = X, p = p, sigma = sigma, nrep = nrep,
      beta1 = b1,               
      s_vals = s_vals, rhos = rhos,
      alpha = alpha, seed_splits = seed_splits
    )
  })
  final_plot(res_list, s_vals, rhos, alpha)
}

## ===== run it =====
set.seed(2025)
n <- 100; p <- 5
X <- matrix(rnorm(n*p), n, p)
X <- 10 * qr.Q(qr(X))   
sigma <- 1
nrep  <- 1000
s_vals <- seq(-1, 1, length.out = 50)
rhos <- c(0.1, 0.5, 0.9)

beta1_vals <- c(.5, .25, 0)
final_plot_wrap(X, p, sigma, nrep, beta1_vals, s_vals, rhos, alpha = 0.05, seed_splits = 1)


pdf("fisher_info.pdf", width = 6.5, height = 6.5)
final_plot_wrap(X, p, sigma, nrep, beta1_vals, s_vals, rhos, alpha = 0.05, seed_splits = 1)
dev.off()

