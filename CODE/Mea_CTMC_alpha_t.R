rm(list = ls())

library(adaptivetau)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

## ============================================================
## 1. Common epidemiological + stress parameters
## ============================================================

params_base <- list(
  # epidemic part
  Lambda = 10,
  mu     = 0.001,
  beta1  = 0.0005,
  beta2  = 0.001,
  gamma1 = 0.10,
  gamma2 = 0.05,
  d1     = 0.001,
  d2     = 0.002,
  # water / stress part
  alpha_max = 0.05,   # maximal stress rate (day^-1)
  Wcrit     = 6,      # critical DO (mg/L)
  kW        = 1,      # sensitivity in S_W(t)
  W0        = 6,      # baseline DO for seasonal case
  A_W       = 2,      # amplitude for seasonal case
  scenario  = NA      # will be set to "wet","medium","dry","seasonal"
)

# initial state
init_state <- c(
  S1 = 300,
  S2 = 10,
  I1 = 1,
  I2 = 0,
  R  = 0
)

t_final    <- 365
dt_plot    <- 1
times_grid <- seq(0, t_final, by = dt_plot)

n_sims     <- 10000   # simulations per scenario
scenarios  <- c("wet", "medium", "dry", "seasonal")

## ============================================================
## 2. Water profile W(t) and alpha(t)
## ============================================================

W_of_t <- function(t, params) {
  with(params, {
    if (scenario == "wet") {
      8.5  # well-oxygenated
    } else if (scenario == "medium") {
      Wcrit  # borderline
    } else if (scenario == "dry") {
      4.5  # low oxygen
    } else if (scenario == "seasonal") {
      # oscillation between ~ Wcrit - A_W and Wcrit + A_W
      W0 + A_W * sin(2 * pi * t / 365)
    } else {
      Wcrit
    }
  })
}

## ============================================================
## 3. CTMC transitions and rate function with alpha(t)
## ============================================================

# State = (S1, S2, I1, I2, R)
transitions_stress <- list(
  c(S1 = +1),                        # birth -> S1
  c(S1 = -1),                        # natural death S1
  c(S2 = -1),                        # natural death S2
  c(I1 = -1),                        # natural death I1
  c(I2 = -1),                        # natural death I2
  c(R  = -1),                        # natural death R
  c(S1 = -1, S2 = +1),               # stress S1 -> S2
  c(S1 = -1, I1 = +1),               # infection S1 -> I1
  c(S2 = -1, I2 = +1),               # infection S2 -> I2
  c(I1 = -1, R  = +1),               # recovery I1 -> R
  c(I2 = -1, R  = +1),               # recovery I2 -> R
  c(I1 = -1),                        # disease death I1
  c(I2 = -1)                         # disease death I2
)

lvrates_stress_time <- function(x, params, t) {
  with(as.list(c(x, params)), {
    # water profile and stress index
    W_t <- W_of_t(t, params)
    S_W <- 1 / (1 + exp(-kW * (Wcrit - W_t)))
    alpha_t <- alpha_max * S_W
    
    lambda_inf <- beta1 * I1 + beta2 * I2
    
    c(
      Lambda,
      mu * S1,
      mu * S2,
      mu * I1,
      mu * I2,
      mu * R,
      alpha_t * S1,    # time-dependent stress rate
      lambda_inf * S1,
      lambda_inf * S2,
      gamma1 * I1,
      gamma2 * I2,
      d1 * I1,
      d2 * I2
    )
  })
}

## ============================================================
## 4. One CTMC trajectory projected on the time grid
## ============================================================

simulate_one_traj_on_grid <- function(params_local) {
  sol <- ssa.exact(
    init.values = init_state,
    transitions = transitions_stress,
    rateFunc    = lvrates_stress_time,
    params      = params_local,
    tf          = t_final
  )
  sol <- as.data.frame(sol)  # columns: time, S1,S2,I1,I2,R
  
  nT <- length(times_grid)
  res <- matrix(NA, nrow = nT, ncol = 5)
  colnames(res) <- c("S1","S2","I1","I2","R")
  
  idx <- 1
  for (j in seq_len(nT)) {
    tgj <- times_grid[j]
    while (idx < nrow(sol) && sol$time[idx + 1] <= tgj) {
      idx <- idx + 1
    }
    res[j, ] <- as.numeric(sol[idx, c("S1","S2","I1","I2","R")])
  }
  
  data.frame(
    time = times_grid,
    S1 = res[,"S1"],
    S2 = res[,"S2"],
    I1 = res[,"I1"],
    I2 = res[,"I2"],
    R  = res[,"R"]
  )
}

## ============================================================
## 5. Mean trajectories (I1, I2) for one scenario
## ============================================================

run_mean_for_scenario <- function(scen) {
  params_local <- params_base
  params_local$scenario <- scen
  
  sum_I1 <- rep(0, length(times_grid))
  sum_I2 <- rep(0, length(times_grid))
  
  for (k in seq_len(n_sims)) {
    df <- simulate_one_traj_on_grid(params_local)
    sum_I1 <- sum_I1 + df$I1
    sum_I2 <- sum_I2 + df$I2
  }
  
  data.frame(
    time    = times_grid,
    I1_mean = sum_I1 / n_sims,
    I2_mean = sum_I2 / n_sims,
    scenario = scen
  )
}

## ============================================================
## 6. Parallel computation over scenarios
## ============================================================

set.seed(123)

n_cores_total <- detectCores()
n_cores_use   <- max(1, floor(0.8 * n_cores_total))

cat("Total cores :", n_cores_total, "\n")
cat("Using cores :", n_cores_use,   "\n")

cl <- makeCluster(n_cores_use)

clusterEvalQ(cl, {
  library(adaptivetau)
})

clusterExport(
  cl,
  varlist = c("params_base","init_state","t_final","times_grid",
              "n_sims","transitions_stress","lvrates_stress_time",
              "W_of_t","simulate_one_traj_on_grid","run_mean_for_scenario"),
  envir = environment()
)

mean_list <- parLapply(cl, scenarios, run_mean_for_scenario)

stopCluster(cl)

mean_all <- bind_rows(mean_list)

## ============================================================
## 7. Plot: only I1 and I2, mean trajectories by scenario
## ============================================================

mean_long <- mean_all %>%
  pivot_longer(
    cols      = c("I1_mean","I2_mean"),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = recode(variable,
                      I1_mean = "Non-stressed infected - I1",
                      I2_mean = "Stressed infected - I2")
  )

# borne max commune (arrondie vers le haut)
ymax <- ceiling(max(mean_long$value, na.rm = TRUE))
p_mean_I <- ggplot(mean_long,
                   aes(x = time,
                       y = value,
                       colour = scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, ncol = 2, scales = "free_y") +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(
    x = "Time (days)",
    y = "Mean number of infected fish",
    colour = "Water scenario",
    title = ""
    #title = "Mean CTMC trajectories under time-dependent stress rate α(t)"
  ) +
  theme_bw()

print(p_mean_I)

# ============================================================
# 8. Save the figure in the 'fig' folder
# ============================================================

ggsave(
  filename = "fig/Mean_ctmc_alpha_t.pdf",
  plot     = p_mean_I,
  width    = 9,
  height   = 4,
  dpi      = 600
)
