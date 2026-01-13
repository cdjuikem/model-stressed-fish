rm(list = ls())

library(adaptivetau)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

## ============================================================
## 1. Paramètres + état initial
## ============================================================

params_base <- list(
  Lambda = 3,
  mu     = 0.001,
  beta1  = 0.0005,
  beta2  = 0.001,
  gamma1 = 0.10,
  gamma2 = 0.05,
  d1     = 0.001,
  d2     = 0.002,
  alpha  = NA      # sera fixé pour chaque alpha
  # (on ajoutera plus tard alpha_max, etc. si on veut alpha(t))
)

init_state <- c(
  S1 = 300,
  S2 = 10,
  I1 = 1,
  I2 = 0,
  R  = 0
)

t_final    <- 365
dt_plot    <- 1        # pas de temps pour la moyenne
times_grid <- seq(0, t_final, by = dt_plot)

n_sims <- 50          # nombre de trajectoires par alpha
alpha_values <- c(0, 0.01, 0.02, 0.04)   # exemple

## ============================================================
## 2. Transitions CTMC et fonction de taux (alpha constant)
## ============================================================

transitions_stress <- list(
  c(S1 = +1),                        # naissance -> S1
  c(S1 = -1),                        # mort nat. S1
  c(S2 = -1),                        # mort nat. S2
  c(I1 = -1),                        # mort nat. I1
  c(I2 = -1),                        # mort nat. I2
  c(R  = -1),                        # mort nat. R
  c(S1 = -1, S2 = +1),               # stress S1 -> S2
  c(S1 = -1, I1 = +1),               # infection S1 -> I1
  c(S2 = -1, I2 = +1),               # infection S2 -> I2
  c(I1 = -1, R  = +1),               # guérison I1
  c(I2 = -1, R  = +1),               # guérison I2
  c(I1 = -1),                        # décès maladie I1
  c(I2 = -1)                         # décès maladie I2
)

lvrates_stress <- function(x, params, t) {
  with(as.list(c(x, params)), {
    lambda_inf <- beta1 * I1 + beta2 * I2
    
    c(
      Lambda,
      mu * S1,
      mu * S2,
      mu * I1,
      mu * I2,
      mu * R,
      alpha * S1,        # ici alpha est CONSTANT pour l’instant
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
## 3. Une trajectoire CTMC, projetée sur la grille de temps
## ============================================================

simulate_one_traj_on_grid <- function(params_local) {
  sol <- ssa.exact(
    init.values = init_state,
    transitions = transitions_stress,
    rateFunc    = lvrates_stress,
    params      = params_local,
    tf          = t_final
  )
  sol <- as.data.frame(sol)   # colonnes: time, S1, S2, I1, I2, R
  
  # Processus en escalier : pour chaque t_grid, on prend le dernier état connu
  nT <- length(times_grid)
  res <- matrix(NA, nrow = nT, ncol = 5)
  colnames(res) <- c("S1","S2","I1","I2","R")
  
  idx <- 1
  for (j in seq_len(nT)) {
    tgj <- times_grid[j]
    # avancer dans la trajectoire tant que le temps suivant <= tgj
    while (idx < nrow(sol) && sol$time[idx + 1] <= tgj) {
      idx <- idx + 1
    }
    res[j, ] <- as.numeric(sol[idx, c("S1","S2","I1","I2","R")])
  }
  
  df <- data.frame(
    time = times_grid,
    S1 = res[,"S1"],
    S2 = res[,"S2"],
    I1 = res[,"I1"],
    I2 = res[,"I2"],
    R  = res[,"R"]
  )
  df
}

## ============================================================
## 4. Moyenne des trajectoires pour une valeur donnée de alpha
## ============================================================

run_mean_for_alpha <- function(alpha_val) {
  params_local <- params_base
  params_local$alpha <- alpha_val
  
  # matrice pour accumuler la somme
  sum_mat <- matrix(0, nrow = length(times_grid), ncol = 5)
  colnames(sum_mat) <- c("S1","S2","I1","I2","R")
  
  for (k in seq_len(n_sims)) {
    df <- simulate_one_traj_on_grid(params_local)
    sum_mat <- sum_mat + as.matrix(df[, c("S1","S2","I1","I2","R")])
  }
  
  mean_mat <- sum_mat / n_sims
  
  data.frame(
    time = times_grid,
    S1   = mean_mat[,"S1"],
    S2   = mean_mat[,"S2"],
    I1   = mean_mat[,"I1"],
    I2   = mean_mat[,"I2"],
    R    = mean_mat[,"R"],
    alpha = alpha_val
  )
}

## ============================================================
## 5. Lancer en parallèle sur les différentes valeurs de alpha
## ============================================================

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
              "n_sims","transitions_stress","lvrates_stress",
              "simulate_one_traj_on_grid","run_mean_for_alpha"),
  envir = environment()
)

mean_list <- parLapply(cl, alpha_values, run_mean_for_alpha)

stopCluster(cl)

mean_all <- bind_rows(mean_list)

## ============================================================
## 6. Plot : trajectoires moyennes pour S1, S2, I1, I2
## ============================================================

mean_long <- mean_all %>%
  pivot_longer(
    cols = c("I1","I2"),
    names_to = "variable",
    values_to = "value"
  )

p_mean <- ggplot(mean_long,
                 aes(x = time,
                     y = value,
                     colour = factor(alpha))) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  labs(
    x = "Time (days)",
    y = "Mean number of fish",
    colour = expression(alpha),
    title = ""
    #title = "Mean CTMC trajectories for I1, I2\nfor different stress rates α"
  ) +
  theme_bw()

print(p_mean)
