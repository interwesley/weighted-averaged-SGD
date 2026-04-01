library(mvtnorm)

## =========================
## Global Settings
## =========================
## This is the code for d = 500 and n = 100000. Set d = 100 and n_iter = 5000 to reproduce the other scenario.
d       <- 500
n_iter  <- 100000
Nrep    <- 500
s       <- 0.1
xt      <- seq(0, 1, length.out = d)
alpha_exp <- 0.8
t_tail <- n_iter / 2L
eta_pda <- 3      # for polynomial-decay averaging


## Linear gradient
grad_linear <- function(x, a, b) {
  diff <- sum(a * x) - b
  diff * a
}

## Sample (a, b)
sample_ab <- function() {
  a <- rnorm(d, mean = 0, sd = sqrt(s))
  b <- sum(a * xt) + rnorm(1)
  list(a = a, b = b)
}


## =================================================
## 1. OW Weighted SGD
## =================================================

run_ow_poly_one <- function() {
  x_old <- rep(0, d)
  
  ## t=2 initialize
  sb2 <- sample_ab()
  a2 <- sb2$a; b2 <- sb2$b
  step2 <- 2^(-alpha_exp)
  g2 <- grad_linear(x_old, a2, b2)
  x_2 <- x_old - step2 * g2
  
  sol <- 2^(alpha_exp - 1) * x_2
  x_old <- x_2
  
  ## t ≥ 3
  for (i in 3:n_iter) {
    sb <- sample_ab()
    a_i <- sb$a; b_i <- sb$b
    
    step_i <- i^(-alpha_exp)
    g_i <- grad_linear(x_old, a_i, b_i)
    x_new <- x_old - step_i * g_i
    
    sol <- ((i - 1)/i) * sol +
      (1 - i^alpha_exp) * (x_old / i) +
      (i^(alpha_exp - 1)) * x_new
    
    x_old <- x_new
  }
  
  mean((sol - xt)^2)
}


## =================================================
## 2. Simple Averaging
## =================================================

run_simple_avg_one <- function() {
  x_old <- rep(0, d)
  avg_x <- rep(0, d)
  
  for (i in 1:n_iter) {
    sb <- sample_ab()
    a_i <- sb$a; b_i <- sb$b
    
    step_i <- i^(-alpha_exp)
    g_i <- grad_linear(x_old, a_i, b_i)
    x_new <- x_old - step_i * g_i
    
    avg_x <- ((i - 1)/i) * avg_x + (1/i) * x_new
    x_old <- x_new
  }
  
  mean((avg_x - xt)^2)
}


## =================================================
## 3. Tail Averaging (decaying step-size)
## =================================================

run_tail_avg_one <- function() {
  x_old <- rep(0, d)
  sum_tail <- rep(0, d)
  
  for (i in 1:n_iter) {
    sb <- sample_ab()
    a_i <- sb$a; b_i <- sb$b
    
    step_i <- i^(-alpha_exp)
    g_i <- grad_linear(x_old, a_i, b_i)
    x_new <- x_old - step_i * g_i
    
    if (i > t_tail) {
      sum_tail <- sum_tail + x_new
    }
    
    x_old <- x_new
  }
  
  tail_avg <- sum_tail / (n_iter - t_tail)
  mean((tail_avg - xt)^2)
}


## =================================================
## 4. Polynomial-Decay Averaging (PDA, η=3)
## =================================================
## w̄_1 = w_1,  and for t > 1:
## w̄_t = (1 - (η+1)/(t+η)) w̄_{t-1} + (η+1)/(t+η) w_t

run_pda_one <- function() {
  x_old <- rep(0, d)
  sb <- sample_ab()
  a1 <- sb$a; b1 <- sb$b
  
  step1 <- 1^(-alpha_exp)
  g1 <- grad_linear(x_old, a1, b1)
  x1 <- x_old - step1 * g1
  
  wbar <- x1
  x_old <- x1
  
  for (t in 2:n_iter) {
    sb <- sample_ab()
    a_i <- sb$a; b_i <- sb$b
    
    step_i <- t^(-alpha_exp)
    g_i <- grad_linear(x_old, a_i, b_i)
    x_new <- x_old - step_i * g_i
    
    wbar <- (1 - (eta_pda + 1)/(t + eta_pda)) * wbar +
      (eta_pda + 1)/(t + eta_pda) * x_new
    
    x_old <- x_new
  }
  
  mean((wbar - xt)^2)
}


## =================================================
## 5. Run 500 trials
## =================================================

set.seed(2025)

mse_ow   <- numeric(Nrep)
mse_simp <- numeric(Nrep)
mse_tail <- numeric(Nrep)
mse_pda  <- numeric(Nrep)

for (j in 1:Nrep) {
  mse_ow[j]   <- run_ow_poly_one()
  mse_simp[j] <- run_simple_avg_one()
  mse_tail[j] <- run_tail_avg_one()
  mse_pda[j]  <- run_pda_one()
  
  if (j %% 50 == 0) cat("rep", j, "done\n")
}

cat("\n===== Linear regression, 500-run results =====\n")
cat("OW weighted:               mean =", mean(mse_ow),
    " sd =", sd(mse_ow), "\n")
cat("Simple averaged:           mean =", mean(mse_simp),
    " sd =", sd(mse_simp), "\n")
cat("Tail-averaged (decay lr):  mean =", mean(mse_tail),
    " sd =", sd(mse_tail), "\n")
cat("Polynomial-decay avg (η=3):mean =", mean(mse_pda),
    " sd =", sd(mse_pda), "\n")



