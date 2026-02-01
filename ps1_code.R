######################################
##
## Problem Set 1
## Question 2
##
## Author: Martin Sielfeld
## Last edition: 01/27/2026
##
######################################

## Settings:
rm(list = ls())
options(scipen = 999)

## Load packages:
library(data.table) ## data manipulation
library(htmltools) ## To display html tables
library(ggplot2) ## figures
library(texreg) ## To format regression tables
library(vars) ## bivariate autocorrelation

## Load data:
houst = fread('houst_detrended.csv')
computsa = fread('COMPUTSA_detrended.csv')


## Sort ddbb:
houst = houst[order(observation_date)]
computsa = computsa[order(observation_date)]

## Create lags:
computsa[, lag1 := shift(COMPUTSA, n = 1, type = 'lag')]
computsa[, lag2 := shift(COMPUTSA, n = 2, type = 'lag')]
computsa[, lag3 := shift(COMPUTSA, n = 3, type = 'lag')]

## Drop newest 6 months:
houst_short = houst[-((nrow(houst) - 5):nrow(houst)), ]
computsa_short = computsa[-((nrow(computsa) - 5):nrow(computsa)), ]

########################################################################################
## Question 1: Fit a univariate AR(3) to completions. Display and discuss the implied ##
##             univariate autocorrelation                                             ##
########################################################################################

## Run model:
lm(COMPUTSA ~ lag1, data = computsa_short)
lm(COMPUTSA ~ lag1 + lag2, data = computsa_short)
lm(COMPUTSA ~ lag1 + lag2 + lag3, data = computsa_short)

## AR(3):
ar3 = lm(COMPUTSA ~ lag1 + lag2 + lag3, data = computsa_short)
summary(ar3)

## Print in HTML:
html_print(HTML(htmlreg(ar3, custom.model.names = c('AR(3)'))))

## The implied ACF is a function of these coefficients only. So we extract
## the coefficiets.
coeffs = coef(ar3)[c('lag1', 'lag2', 'lag3')]

## Now we compute the theoretical ACF implied by an AR moded:
acf = ARMAacf(ar = coeffs, lag.max = 36) ## 12 months * 3 years

## We plot the theoretical ACF:
plot01 =
  ggplot(data.table(x = 0:36, y = acf), aes(x, y)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  coord_cartesian(ylim = c(0.5, 1), expand = F) +
  theme(plot.title.position = 'plot') +
  labs(
    title = "Implied ACF from fitted AR(3) (Completions)",
    x = "Lag (months)",
    y = "Autocorrelation"
  )

plot01

ggsave(
  plot = plot01,
  filename = 'Implied ACF from fitted AR3.png',
  width = 6,
  height = 4,
  dpi = 300
)


##########################################################################################
## Question 2: Use the fitted AR(3) for Wiener-Kolmogorov-Wold 1-step-ahead forecasting ##
##             of the six hold-out completions observations, and assess accuracy.       ##
##########################################################################################

## Predict for the last 6 perdiods
computsa[, predict := predict.lm(ar3, newdata = computsa)]

## Just for plotting reasons, Replace N - 6 with true data:
computsa[-((nrow(computsa) - 5):nrow(computsa)), predict := COMPUTSA]

## Plot:
plot02 =
  ggplot(tail(computsa, 48), aes(x = observation_date, y = predict)) +
  geom_line(aes(
    x = observation_date,
    y = predict,
    linetype = 'AR(3)',
    color = 'AR(3)'
  )) +
  geom_line(aes(
    x = observation_date,
    y = COMPUTSA,
    linetype = 'True value',
    color = 'True value'
  )) +
  scale_linetype_manual(
    values = c('AR(3)' = 'dashed', 'True value' = 'solid')
  ) +
  scale_color_manual(values = c('AR(3)' = 'red', 'True value' = 'black')) +
  theme_bw() +
  theme(plot.title.position = 'plot', legend.position = 'bottom') +
  labs(
    title = "Fitted 6 newest periods using AR(3) model",
    x = "months",
    y = "Total completed houses",
    linetype = NULL,
    color = NULL
  )

plot02

ggsave(
  plot = plot02,
  filename = 'Fitted 6 newest periods using AR3.png',
  width = 6,
  height = 4,
  dpi = 300
)

## Assess accuracy:
test = computsa[((nrow(computsa) - 5):nrow(computsa))]
test[, ar3_err := COMPUTSA - predict]

## Accuracy metrics on holdout
RMSE_ar3 <- sqrt(mean(test$ar3_err^2, na.rm = TRUE))
MAE_ar3 <- mean(abs(test$ar3_err), na.rm = TRUE)
ME_ar3 <- mean(test$ar3_err, na.rm = TRUE)

data.table(model = "AR(3)", RMSE = RMSE_ar3, MAE = MAE_ar3, ME = ME_ar3)


###########################################################################################
## Question 3: Fit a bivariate VAR(3) to starts and completions. Display and discuss the ##
##             implied bivariate autocorrelation                                         ##
###########################################################################################

## Check the ddbb end in the same period:
last(houst$observation_date)
last(computsa$observation_date)

## Prep data for bivariate model:
data = merge(
  computsa[, .(observation_date, COMPUTSA)],
  houst[, .(observation_date, HOUST)],
  by = 'observation_date'
) ## Will drop non-matching dates

## Sort:
data = data[order(observation_date)]

## Drop newest 6 months:
data_short = data[-((nrow(data) - 5):nrow(data)), ]

## VAR(3):
var3 = VAR(data_short[, .(COMPUTSA, HOUST)], p = 3, type = "const")
summary(var3)
roots(var3) # stable if all Mod(roots) > 1

## Print in HTML:
html_print(HTML(htmlreg(
  list(ar3, var3$varresult$COMPUTSA, var3$varresult$HOUST),
  custom.model.names = c(
    'AR(3)  = COMPUTSA',
    'VAR(3) = COMPUTSA',
    'VAR(3) = HOUST'
  ),
  custom.coef.names = c(
    '(Intercept)',
    'lag1' = 'COMPUTSA.l1',
    'lag2' = 'COMPUTSA.l2',
    'lag3' = 'COMPUTSA.l3',
    'COMPUTSA.l1',
    'HOUST.l1',
    'COMPUTSA.l2',
    'HOUST.l2',
    'COMPUTSA.l3',
    'HOUST.l3',
    'const' = '(Intercept)'
  )
)))

## Now we extract the innovations and Wold matrices:
SigmaU = summary(var3)$covres
M = 300
Psi_arr = Psi(var3, nstep = M)

## Now we compute the implied autocovariance matrices:
Kmax = 36 ## 12 monts * 3 years
Gamma_list = vector("list", Kmax + 1)

for (k in 0:Kmax) {
  Gk = matrix(0, 2, 2)
  for (j in 0:(M - k)) {
    Psi_jk = Psi_arr[,, j + k + 1] # Ψ_{j+k}
    Psi_j = Psi_arr[,, j + 1] # Ψ_j
    Gk = Gk + Psi_jk %*% SigmaU %*% t(Psi_j)
  }
  Gamma_list[[k + 1]] <- Gk
}

## Convert Γ(k) to implied autocorrelation matrices R(k):
Gamma0 = Gamma_list[[1]]
Dinvhalf = diag(1 / sqrt(diag(Gamma0)))
R_list = lapply(Gamma_list, function(G) Dinvhalf %*% G %*% Dinvhalf)

acf_dt = rbindlist(lapply(0:Kmax, function(k) {
  Rk = R_list[[k + 1]]
  data.table(
    lag = k,
    `Completions - Completions` = Rk[1, 1],
    `Starts - Starts` = Rk[2, 2],
    `Completions with lagged starts` = Rk[1, 2],
    `Starts with lagged completions` = Rk[2, 1]
  )
}))

acf_dt[, type := "Implied (VAR)"]

## Sample correlations (data)
xC = data_short$COMPUTSA
xH = data_short$HOUST

## Sample ACFs (include lag 0):
acf_C = as.numeric(acf(xC, plot = FALSE, lag.max = Kmax)$acf)
acf_H = as.numeric(acf(xH, plot = FALSE, lag.max = Kmax)$acf)

## Sample cross-correlations:
get_neg_lags = function(ccf_obj, Kmax) {
  lags = as.numeric(ccf_obj$lag)
  vals = as.numeric(ccf_obj$acf)
  out = vals[lags <= 0] # lags: -Kmax,...,0
  out = rev(out) # reorder to 0,1,...,Kmax
  out[1:(Kmax + 1)]
}

cc_CH = ccf(xC, xH, plot = FALSE, lag.max = Kmax) # for Corr(COMPUTSA_t, HOUST_{t-k})
cc_HC = ccf(xH, xC, plot = FALSE, lag.max = Kmax) # for Corr(HOUST_t, COMPUTSA_{t-k})

sample_CH = get_neg_lags(cc_CH, Kmax)
sample_HC = get_neg_lags(cc_HC, Kmax)

acf_sample_dt = data.table(
  lag = 0:Kmax,
  `Completions - Completions` = acf_C,
  `Starts - Starts` = acf_H,
  `Completions with lagged starts` = sample_CH,
  `Starts with lagged completions` = sample_HC
)
acf_sample_dt[, type := "Sample"]

## Combine implied + sample, melt, plot:
acf_both = rbind(acf_dt, acf_sample_dt, fill = TRUE)

acf_long = melt(
  acf_both,
  id.vars = c("lag", "type"),
  variable.name = "pair",
  value.name = "corr"
)

plot03 =
  ggplot(acf_long, aes(lag, corr, linetype = type, color = type)) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c('Implied (VAR)' = 'black', 'Sample' = 'red')) +
  facet_wrap(~pair, ncol = 2) +
  theme_bw() +
  theme(plot.title.position = 'plot', legend.position = 'bottom') +
  labs(
    title = "Bivariate autocorrelations: VAR-implied vs sample",
    x = "Lag (months)",
    y = "Correlation",
    linetype = NULL,
    color = NULL
  )

plot03

ggsave(
  plot = plot03,
  filename = 'Implied bivariate autocorrelations from VAR3.png',
  width = 6,
  height = 4,
  dpi = 300
)

########################################################################################
## Question 4: Perform a full Granger-causality analysis                              ##
########################################################################################

## Test "HOUST Granger-causes COMPUTSA":
gc_H_to_C = causality(var3, cause = "HOUST")
gc_H_to_C$Granger
gc_H_to_C$Instant ## Include “instantaneous causality”

## Test "COMPUTSA Granger-causes HOUST":
gc_C_to_H = causality(var3, cause = "COMPUTSA")
gc_C_to_H$Granger
gc_C_to_H$Instant ## Include “instantaneous causality”


#########################################################################################
## Question 5: Using Cholesky-factor (CF) identification, calculate and graph the full ##
##             set of impulse-response functions. Why is CF identification especially  ##
##             reasonable in a starts/completions VAR?                                 ##
#########################################################################################

set.seed(123) # reproducible bootstrap bands

## Rerun VAR(3) with correct order:
var3_v2 = VAR(data_short[, .(HOUST, COMPUTSA)], p = 3, type = "const")

# Impulse: COMPLETIONS
irf_C = irf(
  var3_v2,
  impulse = "COMPUTSA",
  response = c("HOUST", "COMPUTSA"),
  n.ahead = 36,
  ortho = TRUE,
  boot = TRUE,
  runs = 1000,
  ci = 0.95
)

png("IRF_CF_impulse_COMPUTSA.png", width = 1600, height = 1000, res = 200)
plot(irf_C)
dev.off()

# Impulse: STARTS
irf_H = irf(
  var3_v2,
  impulse = "HOUST",
  response = c("HOUST", "COMPUTSA"),
  n.ahead = 36,
  ortho = TRUE,
  boot = TRUE,
  runs = 1000,
  ci = 0.95
)

png("IRF_CF_impulse_HOUST.png", width = 1600, height = 1000, res = 200)
plot(irf_H)
dev.off()


###########################################################################################
## Question 6: Use the fitted VAR(3) for Wiener-Kolmogorov-Wold 1-step-ahead forecasting ##
##             of the six hold-out completions observations, and assess accuracy. How    ##
##             does the multivariate accuracy compare to the earlier-assessed univariate ##
##             accuracy?                                                                 ##
###########################################################################################

## Produce 1-step-ahead forecasts for the next 6 periods:
pred_var = predict(var3, n.ahead = 6, ci = 0.95)

# Extract point forecasts:
fc_H = as.numeric(pred_var$fcst$HOUST[, "fcst"])
fc_C = as.numeric(pred_var$fcst$COMPUTSA[, "fcst"])

## Keep last 6 periods:
data_test = data[(nrow(data) - 5):nrow(data), ]

## Add forcast to true values:
data_test[, var_fcst := fc_C]
data_test[, var_err := COMPUTSA - var_fcst]

## Assess VAR forecast accuracy on completions:
RMSE_var = sqrt(mean(data_test$var_err^2))
MAE_var = mean(abs(data_test$var_err))
ME_var = mean(data_test$var_err)

data.table(model = "VAR(3)", RMSE = RMSE_var, MAE = MAE_var, ME = ME_var)

## 1-step-ahead forecasts for each of the 6 holdout points using actual lags:
computsa[, ar_err := COMPUTSA - predict]

## Putting VAR and AR forecasts side-by-side:
RMSE_ar = sqrt(mean(computsa[(nrow(computsa) - 5):nrow(computsa)]$ar_err^2))
MAE_ar = mean(abs(computsa[(nrow(computsa) - 5):nrow(computsa)]$ar_err))
ME_ar = mean(computsa[(nrow(computsa) - 5):nrow(computsa)]$ar_err)

## Error measure comparison:
cmp = rbind(
  data.table(model = "AR(3)", RMSE = RMSE_ar, MAE = MAE_ar, ME = ME_ar),
  data.table(model = "VAR(3)", RMSE = RMSE_var, MAE = MAE_var, ME = ME_var)
)
cmp

## Figure:
plot_dt = merge(
  data_test[, .(observation_date, var_fcst)],
  computsa[, .(observation_date, actual = COMPUTSA, predict)],
  by = "observation_date",
  all.y = T
)
plot_dt[nrow(plot_dt) - 6, var_fcst := actual] ## Just for plotting purpuses

plot06 =
  ggplot(tail(plot_dt, 36)) +
  geom_line(aes(
    observation_date,
    var_fcst,
    linetype = 'VAR(3)',
    color = 'VAR(3)'
  )) +
  geom_line(aes(
    observation_date,
    predict,
    linetype = 'AR(3)',
    color = 'AR(3)'
  )) +
  geom_line(aes(
    observation_date,
    actual,
    linetype = 'True value',
    color = 'True value'
  )) +
  scale_linetype_manual(
    values = c(
      'AR(3)' = 'dashed',
      'VAR(3)' = 'dotted',
      'True value' = 'solid'
    )
  ) +
  scale_color_manual(
    values = c(
      'AR(3)' = 'red',
      'VAR(3)' = 'blue',
      'True value' = 'black'
    )
  ) +
  theme_bw() +
  theme(legend.position = "bottom", plot.title.position = 'plot') +
  labs(
    title = "6-month holdout: AR(3) vs VAR(3) forecasts for completions",
    x = NULL,
    y = "COMPUTSA (detrended)",
    linetype = NULL,
    color = NULL
  )

plot06

ggsave(
  plot = plot06,
  filename = '6-month holdout AR3 vs VAR3 forecasts for completions.png',
  width = 6,
  height = 4,
  dpi = 300
)
