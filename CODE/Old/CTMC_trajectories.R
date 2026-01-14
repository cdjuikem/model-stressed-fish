rm(list = ls())

library(adaptivetau)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

## ============================================================
## 1. Paramètres + état initial (modèle stress S1, S2, I1, I2, R)
## ============================================================

params <- list(
  Lambda = 5,      # recrutement (poissons / jour)
  mu     = 0.001,   # mortalité naturelle
  beta1  = 0.0005,  # transmission depuis I1
  beta2  = 0.001,   # transmission depuis I2
  gamma1 = 0.10,    # guérison I1
  gamma2 = 0.05,    # guérison I2
  d1     = 0.001,   # décès maladie I1
  d2     = 0.002,   # décès maladie I2
  alpha  = 0.002     # taux de stress S1 -> S2
)

init_state <- c(
  S1 = 300,   # susceptibles non stressés
  S2 = 10,    # susceptibles stressés
  I1 = 1,     # infectés non stressés (introduction initiale)
  I2 = 0,     # infectés stressés
  R  = 0      # guéris
)

t_final  <- 100      # horizon de simulation (jours)
n_traj   <- 50       # nombre de trajectoires CTMC à tracer

## ============================================================
## 2. Transitions CTMC et fonction de taux
## ============================================================

# Etat = (S1, S2, I1, I2, R)
transitions_stress <- list(
  c(S1 = +1),                        # 1) naissance -> S1
  c(S1 = -1),                        # 2) mort naturelle S1
  c(S2 = -1),                        # 3) mort naturelle S2
  c(I1 = -1),                        # 4) mort naturelle I1
  c(I2 = -1),                        # 5) mort naturelle I2
  c(R  = -1),                        # 6) mort naturelle R
  c(S1 = -1, S2 = +1),               # 7) stress S1 -> S2
  c(S1 = -1, I1 = +1),               # 8) infection S1 -> I1
  c(S2 = -1, I2 = +1),               # 9) infection S2 -> I2
  c(I1 = -1, R  = +1),               # 10) guérison I1 -> R
  c(I2 = -1, R  = +1),               # 11) guérison I2 -> R
  c(I1 = -1),                        # 12) décès maladie I1
  c(I2 = -1)                         # 13) décès maladie I2
)

lvrates_stress <- function(x, params, t) {
  with(as.list(c(x, params)), {
    lambda_inf <- beta1 * I1 + beta2 * I2  # force d'infection commune
    
    c(
      Lambda,            # 1) naissance
      mu * S1,           # 2) mort nat. S1
      mu * S2,           # 3) mort nat. S2
      mu * I1,           # 4) mort nat. I1
      mu * I2,           # 5) mort nat. I2
      mu * R,            # 6) mort nat. R
      alpha * S1,        # 7) stress S1->S2
      lambda_inf * S1,   # 8) infection S1
      lambda_inf * S2,   # 9) infection S2
      gamma1 * I1,       # 10) guérison I1
      gamma2 * I2,       # 11) guérison I2
      d1 * I1,           # 12) décès maladie I1
      d2 * I2            # 13) décès maladie I2
    )
  })
}

## ============================================================
## 3. Une trajectoire CTMC (sera appelée en parallèle)
## ============================================================

simulate_one_traj <- function(traj_id) {
  sol <- ssa.exact(
    init.values = init_state,
    transitions = transitions_stress,
    rateFunc    = lvrates_stress,
    params      = params,
    tf          = t_final
  )
  df <- as.data.frame(sol)
  df$traj <- traj_id
  df
}

set.seed(123)  # pour reproductibilité

## ============================================================
## 4. Simulation en parallèle des trajectoires
## ============================================================

# Nombre de cœurs : 80 % des cœurs dispo
n_cores_total <- detectCores()
n_cores_use   <- max(1, floor(0.8 * n_cores_total))

cat("Total cores :", n_cores_total, "\n")
cat("Using cores :", n_cores_use,   "\n")

cl <- makeCluster(n_cores_use)

# Charger adaptivetau sur chaque nœud
clusterEvalQ(cl, {
  library(adaptivetau)
})

# Exporter les objets nécessaires aux workers
clusterExport(
  cl,
  varlist = c("init_state",
              "params",
              "t_final",
              "transitions_stress",
              "lvrates_stress",
              "simulate_one_traj"),
  envir = environment()
)

traj_ids <- seq_len(n_traj)

# Simuler en parallèle
all_traj_list <- parLapply(cl, traj_ids, simulate_one_traj)

stopCluster(cl)

all_traj <- bind_rows(all_traj_list)

## ============================================================
## 5. Format long + 4 sous-graphes S1, S2, I1, I2
## ============================================================

long_traj <- all_traj %>%
  pivot_longer(
    cols      = c("I1", "I2"),
    names_to  = "variable",
    values_to = "value"
  )

p_ctmc <- ggplot(long_traj,
                 aes(x = time,
                     y = value,
                     group = traj)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~ variable,
             scales = "free_y",
             ncol = 2) +
  labs(
    x = "Time (days)",
    y = "Number of fish",
    title = "Parallel CTMC trajectories for S1, S2, I1, I2"
  ) +
  theme_bw()

print(p_ctmc)
