## ============================================================
## Endemic equilibrium vs alpha (ODE), with scenario colors
## ============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

## --------------------------
## Parameters (match CTMC set)
## --------------------------
par <- list(
  Lambda = 10,
  mu     = 0.001,
  beta1  = 0.0005,
  beta2  = 0.001,
  gamma1 = 0.10,
  gamma2 = 0.05,
  d1     = 0.001,
  d2     = 0.002,
  # stress mapping W -> alpha
  alpha_max = 0.5,
  Wcrit     = 6,
  kW        = 1,
  W0        = 6,
  A_W       = 2
)

## Palette (lock to your choices)
SCEN_LEVELS <- c("dry","medium","seasonal","wet")
PAL_SCEN <- c(
  "dry"      = "#D55E00",   # rouge/orange
  "medium"   = "#009E73",   # vert
  "seasonal" = "#0072B2",   # bleu
  "wet"      = "#CC79A7"    # violet
)

## --------------------------
## Helpers: W(t) and alpha(t)
## --------------------------
S_W_of_W <- function(W, Wcrit, kW) 1/(1 + exp(-kW*(Wcrit - W)))
alpha_of_W <- function(W, alpha_max, Wcrit, kW) alpha_max * S_W_of_W(W, Wcrit, kW)

W_const <- c(
  dry    = 4.5,
  medium = par$Wcrit,
  wet    = 8.5
)

alpha_const <- vapply(names(W_const), function(sc)
  alpha_of_W(W_const[[sc]], par$alpha_max, par$Wcrit, par$kW),
  numeric(1))

## Seasonal alpha(t): min/max over one year
t_year <- seq(0, 365, by = 0.5)
W_seas <- par$W0 + par$A_W * sin(2*pi*t_year/365)
alpha_seasonal <- alpha_of_W(W_seas, par$alpha_max, par$Wcrit, par$kW)
alpha_seasonal_min <- min(alpha_seasonal)
alpha_seasonal_max <- max(alpha_seasonal)

## --------------------------
## ODE system with constant α
## --------------------------
ode_rhs <- function(t, y, p) {
  with(as.list(c(y, p)), {
    lambda <- beta1*I1 + beta2*I2
    dS1 <- Lambda - alpha*S1 - lambda*S1 - mu*S1
    dS2 <- alpha*S1 - lambda*S2 - mu*S2
    dI1 <- lambda*S1 - (gamma1 + mu + d1)*I1
    dI2 <- lambda*S2 - (gamma2 + mu + d2)*I2
    dR  <- gamma1*I1 + gamma2*I2 - mu*R
    list(c(dS1,dS2,dI1,dI2,dR))
  })
}

## Integrate to steady state for a given alpha
solve_to_ss <- function(alpha_val) {
  y0 <- c(S1 = 300, S2 = 10, I1 = 1, I2 = 0, R = 0)
  p  <- c(par, alpha = alpha_val)
  tt <- seq(0, 4000, by = 1)  # long horizon to reach equilibrium
  sol <- ode(y = y0, times = tt, func = ode_rhs, parms = p, atol = 1e-8, rtol = 1e-8)
  tail(as.data.frame(sol), 1)
}

## Sweep α over a grid
alpha_grid <- seq(0, par$alpha_max, length.out = 60)
ss_df <- map_dfr(alpha_grid, function(a) {
  ss <- solve_to_ss(a)
  with(ss, tibble::tibble(
    alpha = a,
    I_total = I1 + I2,
    prop_stressed = (S2 + I2) / (S1 + S2 + I1 + I2 + R),
    mort_eq = par$d1*I1 + par$d2*I2
  ))
})

## --------------------------
## Build ggplot
## --------------------------
if (!dir.exists("FIGS_clo")) dir.create("FIGS_clo", recursive = TRUE)

# Base curves
long <- ss_df |>
  tidyr::pivot_longer(
    cols = c(I_total, prop_stressed, mort_eq),
    names_to = "quantity", values_to = "value"
  ) |>
  dplyr::mutate(
    quantity = factor(quantity,
                      levels = c("I_total","prop_stressed","mort_eq"),
                      labels = c(
                        "Total endemic prevalence  (I1* + I2*)",
                        "Proportion stressed at equilibrium  ((S2*+I2*)/N*)",
                        "Disease-induced mortality at equilibrium  (d1*I1* + d2*I2*)"
                      ))
  )

p <- ggplot(long, aes(x = alpha, y = value)) +
  geom_line(color = "grey30", linewidth = 1) +
  facet_wrap(~ quantity, scales = "free_y", ncol = 1) +
  labs(x = expression(alpha), y = NULL) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

# Add constant-scenario vertical lines
add_vline <- function(p, a, col) {
  p + geom_vline(xintercept = a, color = col, linewidth = 1, alpha = 0.9)
}
p <- add_vline(p, alpha_const["dry"],    PAL_SCEN["dry"])
p <- add_vline(p, alpha_const["medium"], PAL_SCEN["medium"])
p <- add_vline(p, alpha_const["wet"],    PAL_SCEN["wet"])

# Add seasonal alpha(t) span as a horizontal segment on x-axis (cyan)
p <- p +
  annotate("segment",
           x = alpha_seasonal_min, xend = alpha_seasonal_max,
           y = -Inf, yend = -Inf,
           colour = PAL_SCEN["seasonal"], linewidth = 3)

# Legend guide (manual)
leg_df <- tibble::tibble(
  x = c(alpha_const["dry"], alpha_const["medium"], alpha_const["wet"],
        (alpha_seasonal_min + alpha_seasonal_max)/2),
  y = Inf,
  lab = c("dry","medium","wet","seasonal")
)
p <- p +
  geom_point(data = leg_df[1:3,], aes(x = x, y = y, colour = lab),
             inherit.aes = FALSE, shape = 15, size = 3) +
  scale_colour_manual(values = PAL_SCEN, breaks = SCEN_LEVELS, name = "Scenario") +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 4))) +
  theme(legend.position = "right")

print(p)
ggsave("fig/endemic_vs_alpha.png", p, width = 6.5, height = 8.5, dpi = 300)

cat("Saved: fig/endemic_vs_alpha.png\n")
