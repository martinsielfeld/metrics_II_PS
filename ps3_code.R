######################################
##
## Problem Set 3
## Question 3 (Diebold PS - Problem 3)
##
## Author: Martin Sielfeld
## Last edition: 02/19/2026
##
######################################

rm(list = ls())
options(scipen = 999)

## Packages
library(data.table)
library(numDeriv)
library(ggplot2)
library(KFAS)

## Parameters
set.seed(123)
START_DATE <- as.IDate("1970-01-01")
USE_LOG_DIFF <- TRUE # TRUE: growth in log differences
ANNUALIZE <- FALSE # TRUE: 400*diff(log), FALSE: 100*diff(log)
NSIM_STATES <- 300 # posterior draws for bands

############################
## 0) Load + merge levels ##
############################
gdp <- fread("GDP.csv")
gdi <- fread("GDI.csv")

## Robust date parsing:
if (
  !("observation_date" %in% names(gdp)) || !("observation_date" %in% names(gdi))
) {
  stop("Expected column 'observation_date' in both GDP.csv and GDI.csv")
}
gdp[, observation_date := as.IDate(observation_date)]
gdi[, observation_date := as.IDate(observation_date)]

setorder(gdp, observation_date)
setorder(gdi, observation_date)

## Robust rename value columns:
val_col <- function(dt) setdiff(names(dt), "observation_date")[1]
setnames(gdp, val_col(gdp), "GDP")
setnames(gdi, val_col(gdi), "GDI")

data <- merge(
  gdp[, .(observation_date, GDP)],
  gdi[, .(observation_date, GDI)],
  by = "observation_date",
  all = FALSE
)

## Keep 1970+:
data <- data[observation_date >= START_DATE]


####################################################################
## 0.1) Construct quarterly series (growth of GDP and GDI levels) ##
####################################################################
Y_dt <- copy(data)

if (USE_LOG_DIFF) {
  scale_fac <- if (ANNUALIZE) 400 else 100
  Y_dt[, `:=`(
    GDPE = scale_fac * (log(GDP) - shift(log(GDP), 1)),
    GDPI = scale_fac * (log(GDI) - shift(log(GDI), 1))
  )]
} else {
  Y_dt[, `:=`(
    GDPE = GDP - shift(GDP, 1),
    GDPI = GDI - shift(GDI, 1)
  )]
}

## Drop first observation created by diff
Y_dt <- Y_dt[!is.na(GDPE) & !is.na(GDPI)]

## Time fields for ts start
Y_dt[, `:=`(
  year = as.integer(format(observation_date, "%Y")),
  month = as.integer(format(observation_date, "%m"))
)]
Y_dt[, quarter := ((month - 1L) %/% 3L) + 1L]

y_mat <- as.matrix(Y_dt[, .(GDPE, GDPI)])
colnames(y_mat) <- c("GDPE", "GDPI")
y_ts <- ts(y_mat, start = c(Y_dt$year[1], Y_dt$quarter[1]), frequency = 4)

Tobs <- nrow(Y_dt)


###############################################################################
## 1) Model
##
## Measurement:
##   yE_t = GDP_t + eE_t
##   yI_t = GDP_t + eI_t
## State:
##   GDP_t = mu(1-rho) + rho GDP_{t-1} + eG_t
##
## Equivalent coding with x_t = GDP_t - mu:
##   x_t = rho x_{t-1} + eG_t
##   yE_t = mu + x_t + eE_t
##   yI_t = mu + x_t + eI_t
###############################################################################

build_model <- function(y, mu, rho, sigG2, sigE2, sigI2) {
  stopifnot(is.finite(mu), is.finite(rho), abs(rho) < 1)
  stopifnot(
    all(is.finite(c(sigG2, sigE2, sigI2))),
    all(c(sigG2, sigE2, sigI2) > 0)
  )

  n <- nrow(y)

  Z <- array(1, dim = c(2, 1, n))
  Tm <- matrix(rho, 1, 1)
  R <- matrix(1, 1, 1)
  Q <- matrix(sigG2, 1, 1)
  H <- diag(c(sigE2, sigI2), 2, 2)

  ## observation intercept u_t = (mu, mu)'
  u <- array(mu, dim = c(2, 1, n))

  ## stationary init for x_t
  a1 <- matrix(0, 1, 1)
  P1 <- matrix(sigG2 / (1 - rho^2), 1, 1)

  SSModel(
    y ~ -1 +
      SSMcustom(
        Z = Z,
        T = Tm,
        R = R,
        Q = Q,
        a1 = a1,
        P1 = P1,
        P1inf = matrix(0, 1, 1)
      ),
    H = H,
    u = u
  )
}

## par = (mu, p_rho, log_sigG2, log_sigE2, log_sigI2)
unpack_par <- function(par) {
  eps <- 1e-8
  mu <- par[1]
  rho <- max(min(tanh(par[2]), 1 - eps), -1 + eps)
  sigG2 <- exp(par[3])
  sigE2 <- exp(par[4])
  sigI2 <- exp(par[5])
  list(mu = mu, rho = rho, sigG2 = sigG2, sigE2 = sigE2, sigI2 = sigI2)
}

loglik_kf <- function(par, y) {
  th <- unpack_par(par)

  if (!is.finite(th$mu) || !is.finite(th$rho) || abs(th$rho) >= 0.999999) {
    return(-Inf)
  }
  if (
    any(!is.finite(c(th$sigG2, th$sigE2, th$sigI2))) ||
      any(c(th$sigG2, th$sigE2, th$sigI2) <= 0)
  ) {
    return(-Inf)
  }

  mod <- build_model(y, th$mu, th$rho, th$sigG2, th$sigE2, th$sigI2)
  ll <- KFAS::KFS(mod, filtering = "state", smoothing = "none")$logLik
  if (is.null(ll) || !is.finite(ll)) {
    return(-Inf)
  }
  as.numeric(ll)
}

negloglik <- function(par, y) {
  ll <- loglik_kf(par, y)
  if (!is.finite(ll)) {
    return(1e12)
  }
  -ll
}

grad_negloglik <- function(par, y) numDeriv::grad(negloglik, par, y = y)


###################################
## 1(a) MLE: Nelder-Mead -> BFGS ##
###################################
mu0 <- mean(rowMeans(y_mat), na.rm = TRUE)
rho0 <- 0.7

v1 <- var(y_mat[, 1], na.rm = TRUE)
v2 <- var(y_mat[, 2], na.rm = TRUE)

sigG2_0 <- 0.10 * mean(c(v1, v2))
sigE2_0 <- 0.45 * v1
sigI2_0 <- 0.45 * v2

par0 <- c(mu0, atanh(rho0), log(sigG2_0), log(sigE2_0), log(sigI2_0))

mle_nm <- optim(
  par0,
  negloglik,
  y = y_ts,
  method = "Nelder-Mead",
  control = list(maxit = 4000)
)

mle <- optim(
  mle_nm$par,
  negloglik,
  gr = grad_negloglik,
  y = y_ts,
  method = "BFGS",
  control = list(maxit = 4000, reltol = 1e-12)
)

th_hat <- unpack_par(mle$par)

cat("\n====================================\n")
cat("MLE results (Problem 3 growth model)\n")
cat("====================================\n")
cat("mu    :", th_hat$mu, "\n")
cat("rho   :", th_hat$rho, "\n")
cat("sigG2 :", th_hat$sigG2, "\n")
cat("sigE2 :", th_hat$sigE2, "\n")
cat("sigI2 :", th_hat$sigI2, "\n")
cat("convergence code:", mle$convergence, "\n")
cat("max logLik:", -mle$value, "\n\n")

## Hessian SEs (unconstrained)
Hess <- numDeriv::hessian(negloglik, mle$par, y = y_ts)
Vhat <- tryCatch(solve(Hess), error = function(e) NULL)
if (!is.null(Vhat)) {
  se_uncon <- sqrt(diag(Vhat))
  names(se_uncon) <- c("mu", "p_rho", "log_sigG2", "log_sigE2", "log_sigI2")
  cat("Approx SEs (unconstrained params):\n")
  print(se_uncon)
}


##########################################
## Extracted GDP: filtered and smoothed ##
##########################################
model_hat <- build_model(
  y_ts,
  th_hat$mu,
  th_hat$rho,
  th_hat$sigG2,
  th_hat$sigE2,
  th_hat$sigI2
)
kfs_hat <- KFAS::KFS(model_hat, filtering = "state", smoothing = "state")

## Smoothed state: alphahat is x_t (deviation), so GDP_t = mu + x_t
x_smooth_dev <- as.numeric(drop(kfs_hat$alphahat))
Y_dt[, GDP_smooth := th_hat$mu + x_smooth_dev]

## Filtered state (if present in your KFAS version): att is x_{t|t}
GDP_filt <- rep(NA_real_, Tobs)
if (!is.null(kfs_hat$att)) {
  GDP_filt <- th_hat$mu + as.numeric(drop(kfs_hat$att))
}
Y_dt[, GDP_filt := GDP_filt]

## Innovations and standardized innovations
innov <- as.matrix(kfs_hat$v)
Farr <- kfs_hat$F

if (is.matrix(Farr)) {
  F11 <- as.numeric(Farr[1, ])
  F22 <- as.numeric(Farr[2, ])
} else if (length(dim(Farr)) == 3) {
  F11 <- vapply(1:dim(Farr)[3], function(t) Farr[1, 1, t], numeric(1))
  F22 <- vapply(1:dim(Farr)[3], function(t) Farr[2, 2, t], numeric(1))
} else {
  stop("Unexpected format for kfs_hat$F")
}

Y_dt[, `:=`(
  innov_E = as.numeric(innov[, 1]),
  innov_I = as.numeric(innov[, 2]),
  innov_sdE = sqrt(F11),
  innov_sdI = sqrt(F22)
)]


#############
## Figures ##
#############
p_state <-
  ggplot(Y_dt) +
  geom_line(
    mapping = aes(
      x = observation_date,
      y = GDP_filt,
      color = 'Filtered',
      linetype = 'Filtered'
    )
  ) +
  geom_line(
    mapping = aes(
      x = observation_date,
      y = GDP_smooth,
      color = 'Smoothed',
      linetype = 'Smoothed'
    )
  ) +
  scale_linetype_manual(
    values = c('Filtered' = 'solid', 'Smoothed' = 'dashed')
  ) +
  scale_color_manual(
    values = c('Filtered' = 'black', 'Smoothed' = 'red')
  ) +
  theme_bw() +
  labs(
    title = "Extracted GDP growth",
    x = NULL,
    y = "GDP_t",
    color = NULL,
    linetype = NULL
  ) +
  theme(plot.title.position = 'plot', legend.position = 'bottom')

p_state

p_inE <-
  ggplot(Y_dt, aes(observation_date, innov_E / innov_sdE)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    title = "Standardized innovations: GDPE",
    x = NULL,
    y = "vE_t / sqrt(FE_t)"
  )

p_inE

p_inI <-
  ggplot(Y_dt, aes(observation_date, innov_I / innov_sdI)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    title = "Standardized innovations: GDPI",
    x = NULL,
    y = "vI_t / sqrt(FI_t)"
  )

p_inI

ggsave(
  plot = p_state,
  filename = "blue extracted_GDP.png",
  height = 4,
  width = 7,
  dpi = 300
)
ggsave(
  plot = p_inE,
  filename = "blue standardized_innovations_GDPE.png",
  height = 4,
  width = 7,
  dpi = 300
)
ggsave(
  plot = p_inI,
  filename = "blue standardized_innovations_GDPI.png",
  height = 4,
  width = 7,
  dpi = 300
)

## Check values:
t.test(Y_dt$innov_E / Y_dt$innov_sdE)
t.test(Y_dt$innov_I / Y_dt$innov_sdI)
var(Y_dt$innov_E / Y_dt$innov_sdE)
var(Y_dt$innov_I / Y_dt$innov_sdI)


###############################################################################
## 2(b) Bayesian posterior simulation: RW Metropolis-Hastings
###############################################################################
prior <- list(m0 = 0, s0 = 5, a = 2.5, b = 0.10)

log_ig <- function(s, a, b) {
  if (!is.finite(s) || s <= 0) {
    return(-Inf)
  }
  a * log(b) - lgamma(a) - (a + 1) * log(s) - b / s
}

logprior <- function(par) {
  th <- unpack_par(par)

  lp_mu <- dnorm(th$mu, mean = prior$m0, sd = prior$s0, log = TRUE)
  lp_rho <- log(1 - th$rho^2) # Uniform(-1,1) on rho + Jacobian for tanh

  lp_g <- log_ig(th$sigG2, prior$a, prior$b) + log(th$sigG2)
  lp_e <- log_ig(th$sigE2, prior$a, prior$b) + log(th$sigE2)
  lp_i <- log_ig(th$sigI2, prior$a, prior$b) + log(th$sigI2)

  lp_mu + lp_rho + lp_g + lp_e + lp_i
}

logpost <- function(par, y) {
  ll <- loglik_kf(par, y)
  lp <- logprior(par)
  if (!is.finite(ll) || !is.finite(lp)) {
    return(-Inf)
  }
  ll + lp
}

run_mh <- function(
  y,
  par_start,
  n_iter = 30000,
  burn = 8000,
  thin = 10,
  step = c(0.05, 0.06, 0.10, 0.10, 0.10),
  seed = 42
) {
  set.seed(seed)
  p <- length(par_start)
  draws <- matrix(NA_real_, nrow = n_iter, ncol = p)
  colnames(draws) <- c("mu", "p_rho", "log_sigG2", "log_sigE2", "log_sigI2")

  cur <- par_start
  cur_lp <- logpost(cur, y)
  acc <- 0L

  for (i in 1:n_iter) {
    prop <- cur + rnorm(p, 0, step)
    prop_lp <- logpost(prop, y)

    if (is.finite(prop_lp) && log(runif(1)) < (prop_lp - cur_lp)) {
      cur <- prop
      cur_lp <- prop_lp
      acc <- acc + 1L
    }
    draws[i, ] <- cur

    if (i %% 5000 == 0) {
      cat("MH iter", i, "acc rate so far:", round(acc / i, 3), "\n")
    }
  }

  keep <- seq(burn + 1, n_iter, by = thin)
  list(
    draws_raw = draws,
    draws = draws[keep, , drop = FALSE],
    acc_rate = acc / n_iter
  )
}

mh_out <- run_mh(y_ts, par_start = mle$par)
cat("\nFinal MH acceptance rate:", round(mh_out$acc_rate, 3), "\n")

## Transform posterior draws
draws_u <- mh_out$draws
theta_draws <- rbindlist(lapply(1:nrow(draws_u), function(i) {
  th <- unpack_par(draws_u[i, ])
  data.table(
    mu = th$mu,
    rho = th$rho,
    sigG2 = th$sigG2,
    sigE2 = th$sigE2,
    sigI2 = th$sigI2
  )
}))

post_summ <- theta_draws[, .(
  mu_p05 = quantile(mu, 0.05),
  mu_p50 = quantile(mu, 0.50),
  mu_p95 = quantile(mu, 0.95),
  rho_p05 = quantile(rho, 0.05),
  rho_p50 = quantile(rho, 0.50),
  rho_p95 = quantile(rho, 0.95),
  g_p05 = quantile(sigG2, 0.05),
  g_p50 = quantile(sigG2, 0.50),
  g_p95 = quantile(sigG2, 0.95),
  e_p05 = quantile(sigE2, 0.05),
  e_p50 = quantile(sigE2, 0.50),
  e_p95 = quantile(sigE2, 0.95),
  i_p05 = quantile(sigI2, 0.05),
  i_p50 = quantile(sigI2, 0.50),
  i_p95 = quantile(sigI2, 0.95)
)]
cat("\nPosterior quantiles (5/50/95):\n")
print(post_summ)

## Trace plots
trace_dt <- cbind(data.table(iter = 1:nrow(theta_draws)), theta_draws)

p_mu <- ggplot(trace_dt, aes(iter, mu)) +
  geom_line() +
  theme_bw() +
  labs(title = "Trace: mu", x = NULL)
p_rho <- ggplot(trace_dt, aes(iter, rho)) +
  geom_line() +
  theme_bw() +
  labs(title = "Trace: rho", x = NULL)
p_g <- ggplot(trace_dt, aes(iter, sigG2)) +
  geom_line() +
  theme_bw() +
  labs(title = "Trace: sigG2", x = NULL)
p_e <- ggplot(trace_dt, aes(iter, sigE2)) +
  geom_line() +
  theme_bw() +
  labs(title = "Trace: sigE2", x = NULL)
p_i <- ggplot(trace_dt, aes(iter, sigI2)) +
  geom_line() +
  theme_bw() +
  labs(title = "Trace: sigI2", x = NULL)

print(p_mu)
print(p_rho)
print(p_g)
print(p_e)
print(p_i)

ggsave(
  plot = p_mu,
  filename = "blue trace_mu.png",
  height = 4,
  width = 7,
  dpi = 300
)
ggsave(
  plot = p_rho,
  filename = "blue trace_rho.png",
  height = 4,
  width = 7,
  dpi = 300
)
ggsave(
  plot = p_g,
  filename = "blue trace_sigG2.png",
  height = 4,
  width = 7,
  dpi = 300
)
ggsave(
  plot = p_e,
  filename = "blue trace_sigE2.png",
  height = 4,
  width = 7,
  dpi = 300
)
ggsave(
  plot = p_i,
  filename = "blue trace_sigI2.png",
  height = 4,
  width = 7,
  dpi = 300
)

###############################################################################
## Posterior bands for extracted GDP_t = mu + x_t via simulateSSM
###############################################################################
extract_state_draw <- function(sim_obj, Tobs) {
  x <- sim_obj
  if (is.list(sim_obj)) {
    if (!is.null(sim_obj$alphasim)) {
      x <- sim_obj$alphasim
    } else if (!is.null(sim_obj$states)) {
      x <- sim_obj$states
    }
  }
  v <- as.numeric(x)
  n <- length(v)

  if (n == Tobs) {
    return(v)
  }
  if (n == Tobs + 1) {
    return(tail(v, Tobs))
  }
  if (n %% Tobs == 0) {
    return(v[(n - Tobs + 1):n])
  }
  if (n %% (Tobs + 1) == 0) {
    block <- v[(n - (Tobs + 1) + 1):n]
    return(tail(block, Tobs))
  }
  stop(paste0("simulateSSM length ", n, " not compatible with Tobs=", Tobs))
}

idx <- sample.int(nrow(theta_draws), size = min(NSIM_STATES, nrow(theta_draws)))
x_draws_dev <- matrix(NA_real_, nrow = Tobs, ncol = length(idx))
mu_draws <- theta_draws$mu[idx]

for (j in seq_along(idx)) {
  th <- theta_draws[idx[j]]
  mod_j <- build_model(y_ts, th$mu, th$rho, th$sigG2, th$sigE2, th$sigI2)
  sim_j <- KFAS::simulateSSM(mod_j, type = "states", nsim = 1)
  x_draws_dev[, j] <- extract_state_draw(sim_j, Tobs)
  if (j %% 50 == 0) cat("Simulated states:", j, "/", length(idx), "\n")
}

GDP_draws <- sweep(x_draws_dev, 2, mu_draws, FUN = "+")
Y_dt[, `:=`(
  GDP_post_mean = rowMeans(GDP_draws),
  GDP_post_p05 = apply(GDP_draws, 1, quantile, 0.05),
  GDP_post_p95 = apply(GDP_draws, 1, quantile, 0.95)
)]

p_band <- ggplot(Y_dt, aes(observation_date)) +
  geom_ribbon(aes(ymin = GDP_post_p05, ymax = GDP_post_p95), alpha = 0.2) +
  geom_line(aes(y = GDP_post_mean)) +
  theme_bw() +
  labs(
    title = "Posterior mean and 90% band for extracted GDP growth",
    x = NULL,
    y = "GDP_t"
  )

print(p_band)
ggsave(
  plot = p_band,
  filename = "blue posterior_state_band.png",
  height = 4,
  width = 7,
  dpi = 300
)
