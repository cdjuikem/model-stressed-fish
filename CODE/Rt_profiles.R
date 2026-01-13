# --- R(t) with constant vs seasonal alpha(t), preserving your existing colours/linetypes ---

rm(list = ls())
suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})


# --- Directories ---
root_dir <- "C:/Users/cdjuk/OneDrive/Pictures/Desktop/Cours UoM/Paper with stress/CODE"
setwd(root_dir)


# Parameters (consistent with your CTMC/ODE blocks)
params <- list(
  Lambda = 10,  
  mu = 0.001,
  beta1  = 0.00005, 
  beta2 = 0.0001,
  gamma1 = 0.10,   
  gamma2 = 0.05,
  d1 = 0.001, 
  d2 = 0.002,
  alpha_max = 0.05,
  Wcrit = 6, 
  kW = 1,
  W0 = 6, 
  A_W = 2
)

nu1 <- with(params, gamma1 + mu + d1)
nu2 <- with(params, gamma2 + mu + d2)

# Water → stress
S_W_of_W  <- function(W, Wcrit, kW) 1 / (1 + exp(-kW * (Wcrit - W)))
alpha_of_W <- function(W, par) par$alpha_max * S_W_of_W(W, par$Wcrit, par$kW)

# Seasonal W(t) and alpha(t)
W_seasonal <- function(t, par) par$W0 + par$A_W * sin(2 * pi * t / 365)
alpha_t    <- function(t, par) alpha_of_W(W_seasonal(t, par), par)

# R0(alpha) for constant alpha
R0_of_alpha <- function(alpha, par) {
  S1_0 <- par$Lambda / (alpha + par$mu)
  S2_0 <- alpha * par$Lambda / (par$mu * (alpha + par$mu))
  (par$beta1 * S1_0) / nu1 + (par$beta2 * S2_0) / nu2
}

# Constant scenarios (wet/medium/dry) via constant W
W_wet <- 8.5; W_medium <- params$Wcrit; W_dry <- 4.5
alpha_wet    <- alpha_of_W(W_wet,    params)
alpha_medium <- alpha_of_W(W_medium, params)
alpha_dry    <- alpha_of_W(W_dry,    params)
R0_wet    <- R0_of_alpha(alpha_wet,    params)
R0_medium <- R0_of_alpha(alpha_medium, params)
R0_dry    <- R0_of_alpha(alpha_dry,    params)

# DFE ODE for seasonal alpha(t): dS1/dt = Λ - (α(t)+μ)S1; dS2/dt = α(t)S1 - μ S2
rhs_dfe <- function(t, y, par) {
  S1 <- y[1]; S2 <- y[2]; a <- alpha_t(t, par)
  with(par, {
    dS1 <- Lambda - (a + mu) * S1
    dS2 <- a * S1 - mu * S2
    list(c(dS1, dS2))
  })
}

# Start seasonal DFE from steady state at t=0 using alpha(0)
alpha0  <- alpha_t(0, params)
S1_init <- params$Lambda / (alpha0 + params$mu)
S2_init <- alpha0 * params$Lambda / (params$mu * (alpha0 + params$mu))
y0 <- c(S1 = S1_init, S2 = S2_init)

times <- seq(0, 365, by = 1)
sol <- as.data.frame(ode(y = y0, times = times, func = rhs_dfe, parms = params, method = "lsoda"))

# Seasonal R(t)
Rt_seasonal <- with(sol, (params$beta1 * S1) / nu1 + (params$beta2 * S2) / nu2)
df_seasonal <- tibble(time = sol$time, scenario = "seasonal", Rt = Rt_seasonal)

# Constant scenarios as flat series (keep your palette/linetypes by NOT redefining scales)
df_const <- tibble(
  time = rep(times, 3),
  scenario = rep(c("wet", "medium", "dry"), each = length(times)),
  Rt = c(rep(R0_wet, length(times)),
         rep(R0_medium, length(times)),
         rep(R0_dry, length(times)))
)

df_all <- bind_rows(df_const, df_seasonal)

#ymax <- ceiling(max(df_const$value, na.rm = TRUE))

# Plot (no manual scales -> preserves your existing colours/linetypes)
p_Rt <- ggplot(df_all, aes(x = time, y = Rt, colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 1) +
  coord_cartesian(ylim = c(15, 19)) +
  labs(x = "Time (days)", y = expression(R(t)), colour = "Scenario", linetype = "Scenario") +
  theme_bw()

print(p_Rt)

# Save without touching your style
if (!dir.exists("fig")) dir.create("fig", recursive = TRUE)
ggsave("fig/Rt_profiles.pdf", p_Rt, width = 8.5, height = 4.8, dpi = 600)
