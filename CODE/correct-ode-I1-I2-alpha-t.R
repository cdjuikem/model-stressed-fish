# ============================================================
# ODE with time-varying alpha(t) for 4 scenarios
#   wet / medium / dry / seasonal
#   Regime: R0_min > 1  (supercritical)
#   We plot IN (non-stressed infected) and IS (stressed infected)
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
  betaN  = 0.00005,
  betaS  = 0.0001,
  gammaN = 0.10,
  gamma2 = 0.1,
  dN     = 0.005,
  dS     = 0.005,
  # water-stress link
  alpha_max = 0.05,  # max stress rate (day^-1)
  Wcrit     = 6,
  kW        = 1,
  W0        = 6,
  A_W       = 2
)

# Helper: R0(alpha) for the autonomous model
R0_alpha <- function(alpha, p) {
  nuN <- p$gammaN + p$mu + p$dN
  nuS <- p$gamma2 + p$mu + p$dS
  p$Lambda / (alpha + p$mu) *
    ( p$betaN/nuN + alpha * p$betaS/(p$mu * nuS) )
}

alpha_min <- 0.01
alpha_max <- pars$alpha_max

R0_min <- R0_alpha(alpha_min, pars)
R0_max <- R0_alpha(alpha_max, pars)

cat("R0_min (alpha = alpha_min)      =", R0_min, "\n")
cat("R0_max (alpha = alpha_max)      =", R0_max, "\n")

# ---------------- 2) Initial conditions and time grid ----------------
y0 <- c(SN = pars$Lambda/pars$mu, SS = 0, IN = 1, IS = 0, R = 0)

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
    lambda <- betaN * IN + betaS * IS
    
    dSN <- Lambda - a_t*SN - lambda*SN - mu*SN
    dSS <- a_t*SN   - lambda*SS - mu*SS
    dIN <- lambda*SN - (gammaN + mu + dN)*IN
    dIS <- lambda*SS - (gamma2 + mu + dS)*IS
    dR  <- gammaN*IN + gamma2*IS - mu*R
    
    list(c(dSN, dSS, dIN, dIS, dR))
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


# ---------------- 6) Plot IN and IS only ----------------
long_I <- sol_all |>
  select(time, IN, IS, scenario) |>
  pivot_longer(
    cols      = c(IN, IS),
    names_to  = "variable",
    values_to = "value"
  ) |>
  mutate(
    variable = recode(
      variable,
      IN = "'Non-stressed infected'~(I[N])",
      IS = "'Stressed infected'~(I[S])"
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
ggsave("fig/ODE_IN_IS_alpha_t_min_gt1.pdf",
       p_I, width = 9, height = 4, dpi = 300)
