## ============================================================
## CTMC with time-varying alpha(t), first-introduction times
## and extinction probabilities across scenarios
## ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(adaptivetau)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(parallel)
  library(forcats)
})

## --------------------------
## 0) Parameters & scenarios
## --------------------------
params <- list(
  # Epidemic part
  Lambda = 10,
  mu     = 0.001,
  beta1  = 0.00005,  # transmission from I1
  beta2  = 0.0001,   # transmission from I2
  gamma1 = 0.10,    # recovery I1
  gamma2 = 0.1,    # recovery I2
  d1     = 0.01,   # disease death I1
  d2     = 0.01,   # disease death I2
  # Water / stress mapping
  alpha_max = 0.05,  # max stress rate (day^-1)
  Wcrit     = 6,     # mg/L
  kW        = 1,     # sensitivity in S_W(W)
  # Seasonal
  W0        = 6,     # baseline DO
  A_W       = 2      # amplitude
)

# canonical scenario order and palette (lock legend/colors everywhere)
SCEN_LEVELS <- c("dry","medium","seasonal","wet")
PAL_SCEN <- c(
  "dry"      = "#D55E00",  
  "medium"   = "#009E73",   
  "seasonal" = "#0072B2",   
  "wet"      = "#CC79A7"    
)
scale_scenario_colour <- function(...) {
  ggplot2::scale_colour_manual(values = PAL_SCEN, breaks = SCEN_LEVELS,
                               name = "Scenario", drop = FALSE, ...)
}
scale_scenario_fill <- function(...) {
  ggplot2::scale_fill_manual(values = PAL_SCEN, breaks = SCEN_LEVELS,
                             name = "Scenario", drop = FALSE, ...)
}

# Time horizon & sims
t_final  <- 25
n_sims   <- 50000                  # increase when ready (e.g., 2000+)
scenarios <- SCEN_LEVELS       # use the canonical order

# initializers (we'll seed one infected later)
init_S1 <- params$Lambda/(params$mu)
init_S2 <- 0
init_R  <- 0

## ----------------------------------------------
## 1) Stress mapping: W(t) -> S_W(W) -> alpha(t)
## ----------------------------------------------
S_W_of_W <- function(W, Wcrit, kW){
  1 / (1 + exp(-kW * (Wcrit - W)))
}
alpha_of_W <- function(W, alpha_max, Wcrit, kW){
  alpha_max * S_W_of_W(W, Wcrit, kW)
}
W_of_t <- function(t, scenario, par){
  with(par, {
    switch(
      scenario,
      "wet"      = 8.5,                 # well-oxygenated (constant)
      "medium"   = Wcrit,               # borderline (constant)
      "dry"      = 4.5,                 # low oxygen (constant)
      "seasonal" = W0 + A_W * sin(2*pi*t/365),
      Wcrit
    )
  })
}
alpha_t <- function(t, scenario, par){
  Wt <- W_of_t(t, scenario, par)
  alpha_of_W(W = Wt, alpha_max = par$alpha_max, Wcrit = par$Wcrit, kW = par$kW)
}

## ----------------------------------------------
## 2) CTMC transitions & rates (time-dependent)
## ----------------------------------------------
# State = (S1, S2, I1, I2, R)
transitions <- list(
  c(S1 = +1),                        # 1) birth -> S1
  c(S1 = -1),                        # 2) nat. death S1
  c(S2 = -1),                        # 3) nat. death S2
  c(I1 = -1),                        # 4) nat. death I1
  c(I2 = -1),                        # 5) nat. death I2
  c(R  = -1),                        # 6) nat. death R
  c(S1 = -1, S2 = +1),               # 7) stress S1 -> S2
  c(S1 = -1, I1 = +1),               # 8) infection S1 -> I1
  c(S2 = -1, I2 = +1),               # 9) infection S2 -> I2
  c(I1 = -1, R  = +1),               # 10) recovery I1 -> R
  c(I2 = -1, R  = +1),               # 11) recovery I2 -> R
  c(I1 = -1),                        # 12) disease death I1
  c(I2 = -1)                         # 13) disease death I2
)

rate_fun <- function(x, par, t){
  with(as.list(c(x, par)), {
    # time-varying stress
    a_t <- alpha_t(t, scenario = scenario, par = par)
    lambda <- beta1 * I1 + beta2 * I2
    c(
      Lambda,
      mu * S1,
      mu * S2,
      mu * I1,
      mu * I2,
      mu * R,
      a_t * S1,     # stress S1->S2
      lambda * S1,  # inf S1
      lambda * S2,  # inf S2
      gamma1 * I1,
      gamma2 * I2,
      d1 * I1,
      d2 * I2
    )
  })
}
rate_fun_wrapped <- function(x, params, t) rate_fun(x, params, t)

## ---------------------------------------------------------
## 3) One trajectory: seed in I1 or I2, track first-crossing
## ---------------------------------------------------------
simulate_once <- function(scenario, seed_in = c("I1","I2"), par, t_final, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  seed_in <- match.arg(seed_in)
  
  # initial state
  if (seed_in == "I1"){
    x0 <- c(S1 = init_S1, S2 = init_S2, I1 = 1, I2 = 0, R = init_R)
  } else {
    x0 <- c(S1 = init_S1, S2 = init_S2, I1 = 0, I2 = 1, R = init_R)
  }
  
  # attach scenario to params (so rate_fun sees it)
  par$scenario <- scenario
  
  sol <- ssa.exact(
    init.values = x0,
    transitions = transitions,
    rateFunc    = rate_fun_wrapped,
    params      = par,
    tf          = t_final
  )
  sol <- as.data.frame(sol)  # columns: time, S1,S2,I1,I2,R
  
  # First introduction (cross to other class)
  if (seed_in == "I1"){
    idx <- which(sol$I2 > 0)[1]
    tau_cross <- if (is.na(idx)) NA_real_ else sol$time[idx]
  } else {
    idx <- which(sol$I1 > 0)[1]
    tau_cross <- if (is.na(idx)) NA_real_ else sol$time[idx]
  }
  
  # Extinction at final time?
  extinct <- (tail(sol$I1 + sol$I2, 1) == 0)
  
  c(tau_cross = tau_cross, extinct = as.numeric(extinct))
}

## ---------------------------------------------------------
## 4) Parallel batch runner (exports all helpers)
## ---------------------------------------------------------
run_batch <- function(scenario, seed_in = c("I1","I2"), n_sims = 2000L){
  seed_in  <- match.arg(seed_in)
  nb_cores <- max(1, floor(0.8 * parallel::detectCores()))
  cl <- parallel::makeCluster(nb_cores, type = "PSOCK")
  
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(adaptivetau))
    NULL
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "params","t_final","scenario","seed_in","n_sims",
      "simulate_once","transitions",
      "rate_fun_wrapped","rate_fun",
      "alpha_t","W_of_t","alpha_of_W","S_W_of_W",
      "init_S1","init_S2","init_R"
    ),
    envir = environment()
  )
  
  # give each job its own RNG seed for reproducibility
  jobs <- lapply(seq_len(n_sims), function(i) list(i = i))
  out_list <- parallel::parLapply(
    cl, jobs,
    function(job)
      simulate_once(scenario = scenario, seed_in = seed_in,
                    par = params, t_final = t_final, seed = 100000L + job$i)
  )
  parallel::stopCluster(cl)
  
  df <- do.call(rbind, lapply(out_list, function(v) as.data.frame(as.list(v))))
  df$scenario <- scenario
  df$seed_in  <- seed_in
  rownames(df) <- NULL
  tibble::as_tibble(df)
}

## ---------------------------------------------------------
## 5) Run all experiments
## ---------------------------------------------------------
set.seed(123)

message(">>> Running batches... (this may take a few minutes)")
res_list <- list()
for (sc in scenarios){
  res_list[[paste0(sc,"_I1")]] <- run_batch(sc, "I1", n_sims = n_sims)
  res_list[[paste0(sc,"_I2")]] <- run_batch(sc, "I2", n_sims = n_sims)
}
res <- bind_rows(res_list) %>%
  mutate(
    scenario = factor(
      ifelse(as.character(scenario) == "wel", "wet", as.character(scenario)),
      levels = SCEN_LEVELS
    ),
    seed_in = factor(seed_in, levels = c("I1","I2"))
  )

# Quick sanity: how many NA (no cross by t_final)?
summary_NA <- res %>%
  group_by(scenario, seed_in) %>%
  summarise(no_cross = sum(is.na(tau_cross)),
            frac_no_cross = mean(is.na(tau_cross)),
            .groups = "drop")
print(summary_NA)

## ---------------------------------------------------------
## 6) Figures
## ---------------------------------------------------------
if (!dir.exists("fig")) dir.create("fig")

# 6a) First-introduction times (only finite tau_cross)
dens_df <- res %>%
  filter(is.finite(tau_cross)) %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS))

# Base plot with histogram only
p_base <- ggplot(dens_df,
                 aes(x = tau_cross, fill = scenario, colour = scenario)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, alpha = 0.15, position = "identity")

# Build ggplot object and extract computed densities
gb         <- ggplot_build(p_base)
max_density <- max(gb$data[[1]]$density, na.rm = TRUE)*1.9

seed_lab <- c(
  I1 = "'Introduction in'~I[N]",
  I2 = "'Introduction in'~I[S]"
)

# Now make the faceted plot and force the same y-limit in both panels
p_tau <- p_base +
  facet_wrap(~ seed_in, ncol = 1,
             labeller = as_labeller(seed_lab, label_parsed)) +
  scale_scenario_colour() + scale_scenario_fill() +
  labs(x = "First introduction time (days)", y = "Density") +
  coord_cartesian(ylim = c(0, max_density)) +
  theme_bw()

print(p_tau)

ggsave("fig/first_intro_times.pdf", p_tau, width = 9, height = 6, dpi = 300)

# 6b) Extinction probabilities at t_final
ext_df <- res %>%
  group_by(scenario, seed_in) %>%
  summarise(
    p_ext = mean(as.numeric(extinct) == 1),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = SCEN_LEVELS)
  )

p_ext <- ggplot(ext_df, aes(x = scenario, y = p_ext, fill = scenario)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100*p_ext)),
            vjust = -0.6, size = 3.5) +
  facet_wrap(~ seed_in, nrow = 1,
                labeller = as_labeller(seed_lab, label_parsed)) +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_scenario_fill() +
  labs(x = NULL, y = "Extinction probability") +
  theme_bw()
print(p_ext)

ggsave("fig/extinction_bars.pdf", p_ext, width = 9, height = 4, dpi = 300)

message(">>> (Optional) Save enabled lines are commented out. Increase n_sims for publication plots.")

## ---------------------------------------------------------
## 7) (Optional) Compact textual observations in console
## ---------------------------------------------------------
obs <- dens_df %>%
  group_by(scenario, seed_in) %>%
  summarise(
    median_tau = median(tau_cross),
    p25 = quantile(tau_cross, 0.25),
    p75 = quantile(tau_cross, 0.75),
    .groups = "drop"
  ) %>%
  arrange(seed_in, median_tau)

print(obs)
