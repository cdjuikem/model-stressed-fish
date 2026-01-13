# =====================================================
# File: stress_model_simulation_R0_legend.R
# Purpose: Simulate temporal dynamics for 3 regimes
#          and show R0, R01 values in figure legend
# =====================================================

library(deSolve)
library(ggplot2)
library(dplyr)
library(patchwork)

# --- Parameters (biological constants) ---
Lambda <- 10
alpha  <- 0.2
gamma1 <- 0.2
gamma2 <- 0.15
mu     <- 0.01
d1     <- 0.02
d2     <- 0.025

# --- Equilibrium susceptibles at DFE ---
S1_0 <- Lambda / (alpha + mu)
S2_0 <- alpha * S1_0 / mu

# --- Define R0 targets for three scenarios ---
R01_target <- c(0.7, 0.7, 1.3)
R0_target  <- c(0.8, 2, 2.0)

# --- Compute betas ---
cases <- tibble::tibble(
  Case = c("Case 1: Subcritical (R0 < 1)",
           "Case 2: Stress-driven persistence",
           "Case 3: Fully endemic regime"),
  R01 = R01_target,
  R0  = R0_target
) %>%
  mutate(
    R02 = pmax(R0 - R01, 0),
    beta1 = R01 * (gamma1 + mu + d1) / S1_0,
    beta2 = R02 * (gamma2 + mu + d2) / S2_0
  )

print(cases)

# --- ODE System ---
stress_model <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    lambda <- beta1 * I1 + beta2 * I2
    dS1 <- Lambda - alpha*S1 - (lambda + mu)*S1
    dS2 <- alpha*S1 - (lambda + mu)*S2
    dI1 <- lambda*S1 - (gamma1 + mu + d1)*I1
    dI2 <- lambda*S2 - (gamma2 + mu + d2)*I2
    dR  <- gamma1*I1 + gamma2*I2 - mu*R
    list(c(dS1, dS2, dI1, dI2, dR))
  })
}

# --- Initial conditions ---
y0 <- c(S1 = 100, S2 = 20, I1 = 1, I2 = 0.5, R = 0)
times <- seq(0, 1000, by = 1)

# --- Simulations ---
sim_results <- lapply(1:nrow(cases), function(i){
  p <- as.list(cases[i, ])
  parms <- c(Lambda, alpha, gamma1, gamma2, mu, d1, d2,
             beta1 = p$beta1, beta2 = p$beta2)
  
  out <- ode(y = y0, times = times, func = stress_model, parms = parms)
  as.data.frame(out) %>%
    mutate(Case = p$Case,
           R01 = p$R01,
           R0  = p$R0)
}) %>% bind_rows()

# --- Plot function for each case ---
plot_case <- function(df_case, color_I1, color_I2) {
  R01_val <- unique(df_case$R01)
  R0_val  <- unique(df_case$R0)
  
  label_text <- paste0("R[0][1]==", R01_val, "*','~R[0]==", R0_val)
  
  ggplot(df_case, aes(x = time)) +
    geom_line(aes(y = I1), color = color_I1, linewidth = 1.2) +
    geom_line(aes(y = I2), color = color_I2, linewidth = 1.2) +
    annotate("text",
             x = max(df_case$time) * 0.7,
             y = max(df_case$I1 + df_case$I2, na.rm = TRUE) * 0.85,
             label = label_text,
             parse = TRUE,
             color = "black", size = 4.2) +
    labs(
      title = unique(df_case$Case),
      x = "Time (days)",
      y = "Infected"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}


# --- Individual plots ---
#p1 <- plot_case(sim_results %>% filter(Case == cases$Case[1]),
 #               "#0072B2", "#D55E00")
p2 <- plot_case(sim_results %>% filter(Case == cases$Case[2]),
                "#0072B2", "#D55E00")
p3 <- plot_case(sim_results %>% filter(Case == cases$Case[3]),
                "#0072B2", "#D55E00")

# --- Combine in 2x2 layout ---
blank_panel <- ggplot() + theme_void()

final_plot <- p2 + p3 
  plot_annotation(
    title = "Temporal evolution of infected classes under different R₀ regimes",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
    legend.position = "top"  # single legend at bottom
     # plot.subtitle = element_text(size = 12)
    )
  )

# --- Display and Save ---
print(final_plot)

#ggsave("stress_model_three_cases_legend.pdf", final_plot, width = 10, height = 8, dpi = 300)
#ggsave("stress_model_three_cases_legend.png", final_plot, width = 10, height = 8, dpi = 400)
