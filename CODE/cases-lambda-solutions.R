# =====================================================
# File: three_cases_lambda_plot.R
# Purpose: Visualize λ* vs R0 for 3 biologically possible regimes
# =====================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# --- Parameters ---
mu     <- 0.005
alpha  <- 0.1
gamma1 <- 0.2
gamma2 <- 0.25
d1     <- 0.02
d2     <- 0.025

# --- Constants ---
C1 <- (gamma1 + mu + d1) * (gamma2 + mu + d2)
C2 <- C1 * (alpha + 2 * mu)
C3 <- C1 * mu * (alpha + mu)

# --- Function for λ* ---
lambda_star <- function(R0, R01) {
  A2 <- C1
  A1 <- C2 * (1 - R01)
  A0 <- C3 * (1 - R0)
  disc <- A1^2 - 4 * A2 * A0
  if (disc <= 0) return(0)
  roots <- (-A1 + c(1, -1) * sqrt(disc)) / (2 * A2)
  pos_roots <- roots[roots > 0 & is.finite(roots)]
  if (length(pos_roots) == 0) return(0)
  max(pos_roots)
}

# --- R0 range ---
R0_vals <- seq(0, 2.5, length.out = 300)

# --- Define 3 biologically consistent cases ---
cases <- tibble::tribble(
  ~Case, ~R01, ~Color,
  "Case 1: R01 < 1, R0 < 1", 0.7, "#0072B2",
  "Case 2: R01 < 1, R0 > 1", 0.7, "#D55E00",
  "Case 3: R01 > 1, R0 > 1", 1.2, "#CC79A7"
)

# --- Compute λ* ---
df <- cases %>%
  group_by(Case, R01, Color) %>%
  do({
    tibble(
      R0 = R0_vals,
      lambda = sapply(R0_vals, lambda_star, R01 = .$R01)
    )
  }) %>%
  ungroup()

# --- Plot template ---
plot_case <- function(data, title, col){
  R01_val <- unique(data$R01)
  
  ggplot(data, aes(x = R0, y = lambda)) +
    annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf,
             fill = "gray90", alpha = 0.4) +
    geom_line(color = col, linewidth = 1.3) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.9) +
    annotate("text",
             x = max(data$R0) * 0.2,
             y = max(data$lambda, na.rm = TRUE) * 0.9,
             label = paste0("R[0][1]==", R01_val),
             parse = TRUE,
             color = "black", size = 4.5, fontface = "bold") +
    labs(
      title = title,
      x = expression(R[0]),
      y = expression(lambda^"*")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", color = col, size = 14),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# --- Generate subplots ---
p1 <- plot_case(df %>% filter(Case == "Case 1: R01 < 1, R0 < 1"),
                "Case 1: Subcritical", "#0072B2")
p2 <- plot_case(df %>% filter(Case == "Case 2: R01 < 1, R0 > 1"),
                "Case 2: Stress-driven EE", "#D55E00")
p3 <- plot_case(df %>% filter(Case == "Case 3: R01 > 1, R0 > 1"),
                "Case 3: Fully EE", "#CC79A7")

# --- Combine all three panels neatly ---
final_plot <- (p1 + p2 + p3) +
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16)
    )
  )

# --- Display ---
print(final_plot)
ggsave("lambda_star_three_cases.png", final_plot, width = 20, height = 5, dpi = 300)
