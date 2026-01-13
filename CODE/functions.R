# ============================================================
# File: functions.R
# Purpose: One-way stress SIR model with:
#   (i) constant alpha
#   (ii) time-dependent alpha(t) driven by dissolved oxygen W(t)
# ============================================================

# ------------------------------------------------------------
# 1. Force of infection: λ(I1, I2) = beta1 * I1 + beta2 * I2
# ------------------------------------------------------------
lambda_fun <- function(I1, I2, pars) {
  with(as.list(pars), {
    beta1 * I1 + beta2 * I2
  })
}

# ------------------------------------------------------------
# 2. Water and stress index for time-dependent alpha(t)
# ------------------------------------------------------------

# 2.1 Water availability W(t)
# Here we use a simple seasonal pattern around a baseline W0:
#   W(t) = W0 + A_W * sin(2π t / 365)
# You can set A_W = 0 for constant water.
water_fun <- function(t, pars) {
  with(as.list(pars), {
    W <- W0 + A_W * sin(2 * pi * t / 365)
    # Ensure non-negative DO
    max(W, 0)
  })
}

# 2.2 Water-stress index S_W(t) in [0,1]
#   S_W(t) = 1 / (1 + exp(-kW * (Wcrit - W(t))))
stress_index_fun <- function(W, pars) {
  with(as.list(pars), {
    1 / (1 + exp(-kW * (Wcrit - W)))
  })
}

# 2.3 Time-dependent stress rate alpha(t) = alpha_max * S_W(t)
alpha_time_fun <- function(t, state, pars) {
  W_t  <- water_fun(t, pars)
  S_Wt <- stress_index_fun(W_t, pars)
  with(as.list(pars), {
    alpha_max * S_Wt
  })
}

# ------------------------------------------------------------
# 3. One-way stress model with constant alpha
# ------------------------------------------------------------
# State variables:
#   S1 : normal susceptible
#   S2 : stressed susceptible
#   I1 : infected from S1 (non-stressed)
#   I2 : infected from S2 (stressed)
#   R  : recovered
#
# Parameters (for constant-alpha model):
#   Lambda : recruitment (birth) rate
#   alpha  : constant stress rate S1 -> S2
#   mu     : natural death rate
#   beta1  : transmission coefficient from I1
#   beta2  : transmission coefficient from I2
#   gamma1 : recovery rate of I1
#   gamma2 : recovery rate of I2
#   d1     : disease-induced death rate of I1
#   d2     : disease-induced death rate of I2
sir_stress_oneway_const <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    # Force of infection
    lambda <- lambda_fun(I1, I2, pars)
    
    dS1 <- Lambda - alpha * S1 - lambda * S1 - mu * S1
    dS2 <- alpha * S1        - lambda * S2 - mu * S2
    dI1 <- lambda * S1 - (gamma1 + mu + d1) * I1
    dI2 <- lambda * S2 - (gamma2 + mu + d2) * I2
    dR  <- gamma1 * I1 + gamma2 * I2 - mu * R
    
    list(c(dS1, dS2, dI1, dI2, dR))
  })
}

# ------------------------------------------------------------
# 4. One-way stress model with time-dependent alpha(t)
# ------------------------------------------------------------
# Same state variables and infection parameters, but:
#   alpha(t) = alpha_time_fun(t, state, pars)
# Parameters (extra, for time-dependent alpha):
#   alpha_max : maximal stress rate
#   W0        : baseline dissolved oxygen
#   A_W       : amplitude of seasonal DO variation
#   Wcrit     : critical DO level
#   kW        : sensitivity of stress to oxygen deficit
sir_stress_oneway_time <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    # Time-dependent stress rate alpha(t)
    alpha_t <- alpha_time_fun(t, state, pars)
    
    # Force of infection
    lambda <- lambda_fun(I1, I2, pars)
    
    dS1 <- Lambda - alpha_t * S1 - lambda * S1 - mu * S1
    dS2 <- alpha_t * S1         - lambda * S2 - mu * S2
    dI1 <- lambda * S1 - (gamma1 + mu + d1) * I1
    dI2 <- lambda * S2 - (gamma2 + mu + d2) * I2
    dR  <- gamma1 * I1 + gamma2 * I2 - mu * R
    
    list(c(dS1, dS2, dI1, dI2, dR))
  })
}

