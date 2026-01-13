library(ggplot2)
library(patchwork)

# ============================================================
# PARAMETRES COMMUNS
# ============================================================
alpha_max <- 0.05   # max stress rate (day^-1)
W_crit    <- 6      # critical DO (mg/L)
k_W       <- 1      # sensitivity (1 / (mg/L))

times <- seq(0, 365, by = 1)  # 1 year, in days

# --- Directories ---
root_dir <- "C:/Users/cdjuk/OneDrive/Pictures/Desktop/Cours UoM/Paper with stress/CODE"
setwd(root_dir)

# ------------------------------------------------------------
# Fonction index de stress S_W(W)
# ------------------------------------------------------------
stress_index <- function(W, W_crit, k_W) {
  1 / (1 + exp(-k_W * (W_crit - W)))
}

# ============================================================
# FIGURE 2 : 4 SCENARIOS (wet / medium / dry / seasonal)
# Panel 1 : W(t) ; Panel 2 : alpha(t)
# ============================================================

water_wet <- function(t) {
  8.5   # > W_crit : bonne oxygénation
}

water_medium <- function(t) {
  6     # = W_crit : limite
}

water_dry <- function(t) {
  4.5   # < W_crit : conditions pauvres
}

water_seasonal <- function(t) {
  6 + 2 * sin(2 * pi * t / 365)  # entre ~4 et 8 mg/L
}

# Construire le data.frame pour les 4 scénarios
df_scen <- data.frame(
  time     = rep(times, times = 4),
  scenario = rep(c("wet", "medium", "dry", "seasonal"), each = length(times)),
  W        = NA,
  alpha    = NA
)

for (i in seq_along(times)) {
  t <- times[i]
  df_scen$W[df_scen$scenario == "wet"      & df_scen$time == t] <- water_wet(t)
  df_scen$W[df_scen$scenario == "medium"   & df_scen$time == t] <- water_medium(t)
  df_scen$W[df_scen$scenario == "dry"      & df_scen$time == t] <- water_dry(t)
  df_scen$W[df_scen$scenario == "seasonal" & df_scen$time == t] <- water_seasonal(t)
}

df_scen$S_W   <- stress_index(df_scen$W, W_crit, k_W)
df_scen$alpha <- alpha_max * df_scen$S_W

# --- Palette commune pour toutes les figures ---
df_scen$scenario <- factor(
  df_scen$scenario,
  levels = c("dry", "medium", "seasonal", "wet")
)

scenario_cols <- c(
  "dry"      = "#D55E00",   # rouge/orange
  "medium"   = "#009E73",   # vert
  "seasonal" = "#0072B2",   # bleu
  "wet"      = "#CC79A7"    # violet
)


# Panel du haut : W(t) pour chaque scénario
p_W_scen <- ggplot(df_scen, aes(x = time, y = W, colour = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = scenario_cols, name = "Stress scenario") +
  labs(
    x = "Time (days)",
    y = expression(W(t)),
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Panel du bas : alpha(t) pour chaque scénario
p_alpha_scen <- ggplot(df_scen, aes(x = time, y = alpha, colour = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = scenario_cols, name = "Stress scenario") +
  labs(
    x = "Time (days)",
    y = expression(alpha(t)),
    title = ""
  ) +
  theme_minimal()

# Figure 2 : 2 sous-plots (ici côte à côte ; mets / si tu veux empiler)
fig2 <- p_W_scen | p_alpha_scen
print(fig2)

ggsave(
  filename = "fig/fig2_W_alpha_scenarios.pdf",
  plot     = fig2,
  width    = 10,
  height   = 3,
  dpi      = 300
)
