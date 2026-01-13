## Packages ----
library(deSolve)
library(dplyr)
library(ggplot2)


# --- Directories ---
root_dir <- "C:/Users/cdjuk/OneDrive/Pictures/Desktop/Cours UoM/Paper with stress/CODE"
setwd(root_dir)


## Parameters (example values – replace with yours) ----
Lambda  <- 10      # recruitment
mu      <- 0.005     # natural death
beta1   <- 0.00005 
beta2   <- 0.0001
gamma1  <- 0.1
gamma2  <- 0.1
d1      <- 0.01
d2      <- 0.010

nu1 <- gamma1 + mu + d1
nu2 <- gamma2 + mu + d2

alpha_min <- 0.01
alpha_max <- 0.05

## R0(alpha) for the autonomous model ----
R0_alpha <- function(alpha) {
  Lambda/(alpha + mu) * ( beta1/nu1 + alpha*beta2/(mu*nu2) )
}

R0_min <- R0_alpha(alpha_min)
R0_max <- R0_alpha(alpha_max)

## Time-varying stress profiles alpha(t) ----
alpha_fun_list <- list(
  wet = function(t) rep(alpha_min, length(t)),
  dry = function(t) rep(alpha_max, length(t)),
  medium = function(t) rep(0.5*(alpha_min + alpha_max), length(t)),
  seasonal = function(t) {
    # oscillates between alpha_min and alpha_max over a 1-year period
    m  <- 0.5*(alpha_min + alpha_max)
    a  <- 0.5*(alpha_max - alpha_min)
    m + a * sin(2*pi*t/365)
  }
)

scenario_cols <- c(
      "dry"      = "#D55E00",   # rouge/orange
      "medium"   = "#009E73",   # vert
      "seasonal" = "#0072B2",   # bleu
      "wet"      = "#CC79A7"    # violet
    )


## ODE for the disease-free trajectory (S1^0, S2^0) ----
dS_dt <- function(t, state, parms) {
  S1 <- state["S1"]
  S2 <- state["S2"]
  
  alpha_t <- parms$alpha_fun(t)
  
  with(parms, {
    dS1 <- Lambda - (alpha_t + mu)*S1
    dS2 <- alpha_t*S1 - mu*S2
    list(c(S1 = dS1, S2 = dS2))
  })
}

## Time grid ----
tmax  <- 365*5       # 3 years, say
times <- seq(0, tmax, by = 0.1)

## Simulate R(t) for each alpha(t) scenario ----
Rt_df <- lapply(names(alpha_fun_list), function(scen) {
  pars <- list(
    Lambda   = Lambda,
    mu       = mu,
    alpha_fun = alpha_fun_list[[scen]]
  )
  
  # some initial condition; after a transient it doesn't matter much
  y0 <- c(S1 = Lambda/mu, S2 = 0)
  
  out <- ode(
    y     = y0,
    times = times,
    func  = dS_dt,
    parms = pars
  )
  
  out <- as.data.frame(out)
  out$scenario <- scen
  
  # instantaneous reproduction number along the DFE trajectory
  out$R_t <- beta1 * out$S1 / nu1 + beta2 * out$S2 / nu2
  out
}) |> bind_rows()

# Range for padding
y_range <- max(Rt_df$R_t) - min(Rt_df$R_t)
x_max   <- max(Rt_df$time)

p <- ggplot(Rt_df, aes(x = time, y = R_t, colour = scenario)) +
  geom_hline(yintercept = R0_min, linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = R0_max, linetype = "dashed", linewidth = 1) +
  geom_line(linewidth = 1) +
  scale_colour_manual(
    values = scenario_cols,
    name   = "Stress scenario"
  ) +
  annotate(
    "text",
    x = x_max * 0.98,
    y = R0_min - 0.08 * y_range,
    label = paste0("R[0](alpha[min])==", round(R0_min, 2)),
    hjust = 1, vjust = 0,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = x_max * 0.98,
    y = R0_max + 0.01 * y_range,
    label = paste0("R[0](alpha[max])==", round(R0_max, 2)),
    hjust = 1, vjust = 0,
    parse = TRUE
  ) +
  labs(
    x = "Time (days)",
    y = expression(R(t)),
    colour = "Stress scenario"
  ) +
  theme_bw()

print(p)

ggsave("fig/Rt_profiles_ED0.pdf", p, width = 9, height = 4, dpi = 300)
