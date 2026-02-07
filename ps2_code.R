######################################
##
## Problem Set 2
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

#########################################################################################
## Question 1: Fit a univariate AR(3) to completions. Display and discuss the spectral ##
##             density functions.                                                      ##
#########################################################################################

sp1 = spec.ar(computsa_short$COMPUTSA, order = 3, plot = TRUE)
sp1 = data.table(
  omega = sp1$freq * 2 * pi,
  spec = sp1$spec[, 1] / (2 * pi), # convert density-per-cycle to density-per-radian
  model = "AR(3) completions"
)

plot01 =
  ggplot(sp1, aes(x = omega, y = log(spec))) +
  geom_line(color = 'black') +
  theme_bw() +
  theme(plot.title.position = 'plot') +
  labs(
    x = expression(omega ~ "(radians)"),
    y = 'Log(Spectrum)',
    title = 'Completions spectrum in radians: AR(3)'
  )

plot01

ggsave(
  plot = plot01,
  filename = 'blue Completions spectrum in radians AR3.png',
  width = 6,
  height = 4,
  dpi = 300
)

#######################################################################################
## Question 2: Fit a bivariate VAR(3) to starts and completions. Display and discuss ##
##             the implied bivariate autocorrelation and spectral density functions. ##
#######################################################################################

set.seed(123) # reproducible bootstrap bands

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

# Fit VAR(3) to starts & completions
var3 = VAR(data_short[, .(HOUST, COMPUTSA)], p = 3, type = "const")

A <- Acoef(var3)
U <- resid(var3)
SigU <- crossprod(U) / nrow(U)

spec_var <- function(omega, A, SigU) {
  k <- nrow(A[[1]])
  z <- exp(-1i * omega)
  Az <- diag(k) - A[[1]] * z - A[[2]] * z^2 - A[[3]] * z^3
  H <- solve(Az)
  (1 / (2 * pi)) * H %*% SigU %*% Conj(t(H))
}

omega <- seq(0, pi, length.out = 1024)
k <- ncol(U)
Fomega <- array(NA_complex_, dim = c(k, k, length(omega)))
for (j in seq_along(omega)) {
  Fomega[,, j] <- spec_var(omega[j], A, SigU)
}

## Total completions spectrum:
sp2 <- data.table(
  omega = omega,
  spec = Re(Fomega[2, 2, ]), # completions = 2nd variable
  model = "VAR(3) implied completions"
)

## Total starts spectrum:
sp3 = data.table(
  omega = omega,
  spec = Re(Fomega[1, 1, ]), # completions = 2nd variable
  model = "VAR(3) implied starts"
)

## Append
sp = rbind(sp1, sp2, sp3)

## Compare f11 and f22 VAR(3):
plot02 =
  ggplot(
    sp[model %in% c('VAR(3) implied starts', 'VAR(3) implied completions')],
    aes(x = omega, y = log(spec), color = model)
  ) +
  geom_line() +
  scale_color_manual(
    values = c(
      'VAR(3) implied starts' = 'green',
      'VAR(3) implied completions' = 'red'
    )
  ) +
  theme_bw() +
  theme(plot.title.position = 'plot', legend.position = 'bottom') +
  labs(
    x = expression(omega ~ "(radians)"),
    y = "Log(Spectrum)",
    title = "VAR(3) implied completions and starts spectrum in radians",
    color = NULL
  )

plot02

ggsave(
  plot = plot02,
  filename = 'blue VAR3 implied completions and starts spectrum in radians.png',
  width = 6,
  height = 4,
  dpi = 300
)

## Compare AR(3) to VAR(3) spectral:
plot03 =
  ggplot(
    sp[model %in% c('AR(3) completions', 'VAR(3) implied completions')],
    aes(x = omega, y = log(spec), color = model)
  ) +
  geom_line() +
  scale_color_manual(
    values = c(
      'AR(3) completions' = 'black',
      'VAR(3) implied completions' = 'red'
    )
  ) +
  theme_bw() +
  theme(plot.title.position = 'plot', legend.position = 'bottom') +
  labs(
    x = expression(omega ~ "(radians)"),
    y = "Log(Spectrum)",
    title = "Completions spectrum in radians: AR(3) vs VAR(3)-implied",
    color = NULL
  )

plot03

ggsave(
  plot = plot03,
  filename = 'blue Completions spectrum in radians AR3 vs VAR3-implied.png',
  width = 6,
  height = 4,
  dpi = 300
)

f11 <- Re(Fomega[1, 1, ])
f22 <- Re(Fomega[2, 2, ])
f12 <- Fomega[1, 2, ]

coh12 <- (Mod(f12)^2) / (f11 * f22)

plot04 =
  ggplot(data.table(omega = omega, coh = coh12), aes(x = omega, y = coh)) +
  geom_line(color = 'blue') +
  theme_bw() +
  theme(plot.title.position = 'plot', legend.position = 'bottom') +
  labs(
    x = expression(omega ~ "(radians)"),
    y = "Coherence",
    title = "VAR(3): coherence between starts and completions"
  )

plot04

ggsave(
  plot = plot04,
  filename = 'blue VAR3 coherence between starts and completions.png',
  width = 6,
  height = 4,
  dpi = 300
)
