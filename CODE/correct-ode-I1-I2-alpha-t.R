# ============================================================
# ODE with time-varying alpha(t) for 4 scenarios
#   wet / medium / dry / seasonal
#   Regime: R0_min > 1  (supercritical)
#   We plot I1 (non-stressed infected) and I2 (stressed infected)
# ============================================================

rm(list = ls())
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------- 1) Parameters (chosen so that R0_min > 1) ----------------
pars <- list(
  Lambda = 10,
  mu     = 0.005,
  beta1  = 0.00005,
  beta2  = 0.0001,
  gamma1 = 0.10,
  gamma2 = 0.1,
  d1     = 0.01,
  d2     = 0.01,
  # water-stress link
  alpha_max = 0.05,  # max stress rate (day^-1)
  Wcrit     = 6,
  kW        = 1,
  W0        = 6,
  A_W       = 2
)

# Helper: R0(alpha) for the autonomous model
R0_alpha <- function(alpha, p) {
  nu1 <- p$gamma1 + p$mu + p$d1
  nu2 <- p$gamma2 + p$mu + p$d2
  p$Lambda / (alpha + p$mu) *
    ( p$beta1/nu1 + alpha * p$beta2/(p$mu * nu2) )
}

alpha_min <- 0.01
alpha_max <- pars$alpha_max

R0_min <- R0_alpha(alpha_min, pars)
R0_max <- R0_alpha(alpha_max, pars)

cat("R0_min (alpha = alpha_min)      =", R0_min, "\n")
cat("R0_max (alpha = alpha_max)      =", R0_max, "\n")

# ---------------- 2) Initial conditions and time grid ----------------
y0 <- c(S1 = pars$Lambda/pars$mu, S2 = 0, I1 = 1, I2 = 0, R = 0)

t_end <- 365*4
times <- seq(0, t_end, by = 0.5)

scenarios <- c("wet","medium","dry","seasonal")

# ---------------- 3) W(t) and alpha(t) ----------------
W_of_t <- function(t, scen, p) {
  if (scen == "wet")      return(8.5)
  if (scen == "medium")   return(p$Wcrit)
  if (scen == "dry")      return(4.5)
  if (scen == "seasonal") return(p$W0 + p$A_W * sin(2*pi*t/365))
  p$Wcrit
}

alpha_of_t <- function(t, scen, p) {
  Wt  <- W_of_t(t, scen, p)
  SWt <- 1 / (1 + exp(-p$kW * (p$Wcrit - Wt)))  # in [0,1]
  p$alpha_max * SWt
}

# ---------------- 4) RHS of ODE system ----------------
rhs <- function(t, y, p, scen) {
  with(as.list(c(y, p)), {
    a_t   <- alpha_of_t(t, scen, p)
    lambda <- beta1 * I1 + beta2 * I2
    
    dS1 <- Lambda - a_t*S1 - lambda*S1 - mu*S1
    dS2 <- a_t*S1   - lambda*S2 - mu*S2
    dI1 <- lambda*S1 - (gamma1 + mu + d1)*I1
    dI2 <- lambda*S2 - (gamma2 + mu + d2)*I2
    dR  <- gamma1*I1 + gamma2*I2 - mu*R
    
    list(c(dS1, dS2, dI1, dI2, dR))
  })
}

# ---------------- 5) Solve for each scenario ----------------
solve_one <- function(scen) {
  out <- ode(
    y     = y0,
    times = times,
    func  = function(t, y, parms) rhs(t, y, parms, scen),
    parms = pars,
    method = "lsoda"
  )
  as.data.frame(out) |>
    mutate(scenario = scen)
}

sol_all <- do.call(rbind, lapply(scenarios, solve_one))

# --- Set factor order + colours consistent with other figures ----
sol_all$scenario <- factor(
  sol_all$scenario,
  levels = c("dry", "medium", "seasonal", "wet")
)

scenario_cols <- c(
  "dry"      = "#D55E00",   # rouge/orange
  "medium"   = "#009E73",   # vert
  "seasonal" = "#0072B2",   # bleu
  "wet"      = "#CC79A7"    # violet
)


# ---------------- 6) Plot I1 and I2 only ----------------
long_I <- sol_all |>
  select(time, I1, I2, scenario) |>
  pivot_longer(
    cols      = c(I1, I2),
    names_to  = "variable",
    values_to = "value"
  ) |>
  mutate(
    variable = recode(
      variable,
      I1 = "'Non-stressed infected'~(I[1])",
      I2 = "'Stressed infected'~(I[2])"
    )
  )

ymax <- ceiling(max(long_I$value, na.rm = TRUE))

p_I <- ggplot(long_I,
              aes(x = time, y = value, colour = scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scenario_cols, name = "Scenario") +
  facet_wrap(
    ~ variable,
    ncol = 2,
    scales = "free_y",
    labeller = labeller(variable = label_parsed)
  ) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(
    x = "Time (days)",
    y = "Number of infected fish"
  ) +
  theme_bw()

print(p_I)

# ---------------- 7) Save figure ----------------
if (!dir.exists("fig")) dir.create("fig")
ggsave("fig/ODE_I1_I2_alpha_t_min_gt1.pdf",
       p_I, width = 9, height = 4, dpi = 300)
