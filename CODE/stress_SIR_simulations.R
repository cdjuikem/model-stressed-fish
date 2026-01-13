library(deSolve)
library(ggplot2)
library(patchwork)

# Load the ODE models and alpha(t) helpers
source("functions.R")

# ============================================================
# 1. Time grid and initial conditions
# ============================================================
times <- seq(0, 365, by = 0.5)   # one year

init_state <- c(
  S1 = 300,   # normal susceptible fish
  S2 = 10,   # stressed susceptible fish
  I1 = 1,     # initially infected (non-stressed)
  I2 = 0,     # initially infected (stressed)
  R  = 0      # recovered
)

# ============================================================
# 2. Common epidemiological parameters (fish-inspired)
#    (orders of magnitude consistent with fish diseases)
# ============================================================
pars_common <- c(
  Lambda = 10,      # recruitment (fish/day)
  mu     = 0.001,   # natural mortality (day^-1)
  beta1  = 0.0005,  # transmission from I1
  beta2  = 0.001,   # transmission from I2 (higher if stress ↑ susceptibility)
  gamma1 = 0.10,    # recovery from I1 (day^-1)
  gamma2 = 0.05,    # recovery from I2 (slower recovery)
  d1     = 0.001,   # disease-induced death I1 (day^-1)
  d2     = 0.002    # disease-induced death I2 (day^-1)
)

# ============================================================
# 3. Parameters for CONSTANT alpha
# ============================================================
pars_const <- c(
  pars_common,
  alpha = 0.02      # constant stress rate S1 -> S2 (day^-1)
)

# ============================================================
# 4. Parameters for TIME-DEPENDENT alpha(t)
#    alpha(t) = alpha_max * S_W(t),
#    W(t) seasonal around W0 = Wcrit
# ============================================================
pars_time <- c(
  pars_common,
  alpha_max = 0.05,  # maximal stress rate (day^-1)
  W0        = 6,     # baseline DO (mg/L)
  A_W       = 2,     # amplitude of seasonal variation (mg/L)
  Wcrit     = 6,     # critical DO for stress (mg/L)
  kW        = 1      # sensitivity (1 / (mg/L))
)

# ============================================================
# 5. Run the ODE simulations
# ============================================================
out_const <- ode(
  y     = init_state,
  times = times,
  func  = sir_stress_oneway_const,
  parms = pars_const
)
out_const <- as.data.frame(out_const)
out_const$scenario <- "constant alpha"

out_time <- ode(
  y     = init_state,
  times = times,
  func  = sir_stress_oneway_time,
  parms = pars_time
)
out_time <- as.data.frame(out_time)
out_time$scenario <- "time-dependent alpha(t)"

# Combine trajectories
out_all <- rbind(out_const, out_time)


y_lim_S <- range(c(out_all$S1, out_all$S2), na.rm = TRUE)
y_lim_I <- range(c(out_all$I1, out_all$I2), na.rm = TRUE)

# ============================================================
# 6. Helper: month ticks for nicer x-axis
# ============================================================
month_breaks <- seq(0, 360, by = 30)
month_labels <- month.abb  # "Jan", "Feb", ...

# ============================================================
# 7. Build the four panels S1, S2, I1, I2  (x in DAYS)
#    (each with both scenarios)
# ============================================================
make_panel <- function(var, ylab) {
  ggplot(out_all, aes(x = time, y = .data[[var]], colour = scenario)) +
    geom_line(size = 1) +
    labs(
      x = "Time (days)",
      y = ylab
    ) +
    theme_minimal() +
    coord_cartesian(ylim = y_lim_S) +
    theme(legend.position = "none")
}

p_S1 <- make_panel("S1", "S1")
p_S2 <- make_panel("S2", "S2")
p_I1 <- make_panel("I1", "I1")
p_I2 <- make_panel("I2", "I2")

# Arrange as 4 subplots: S1 S2 on top, I1 I2 on bottom
fig_states <- (p_S1 | p_S2) / (p_I1 | p_I2)+
  plot_layout(guides = "collect") &   # collect the guides (legend)
  theme(legend.position = "top")   # single legend at bottom

print(fig_states)

# ============================================================
# 8. Save the figure in the 'fig' folder
# ============================================================
if (!dir.exists("fig")) dir.create("fig")

ggsave(
  filename = "fig/fig_states_const_vs_alpha_t.png",
  plot     = fig_states,
  width    = 9,
  height   = 6,
  dpi      = 300
)
