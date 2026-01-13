## ============================================================
## Violin plots of extinction probability vs alpha
## P_ext(alpha, i10, i20) = (u1_star(alpha))^i10 * (u2_star(alpha))^i20
## Two seeding scenarios: seed in I1 only, seed in I2 only
## ============================================================

rm(list = ls())

library(dplyr)
library(ggplot2)

## --- Parameters (your values) -------------------------------
Lambda  <- 10        # recruitment
mu      <- 0.005     # natural death
beta1   <- 0.00005   # transmission from I1
beta2   <- 0.0001    # transmission from I2
gamma1  <- 0.1       # recovery I1
gamma2  <- 0.1       # recovery I2
d1      <- 0.01      # disease death I1
d2      <- 0.01      # disease death I2

nu1 <- gamma1 + mu + d1
nu2 <- gamma2 + mu + d2

## --- Choose alpha values (those you used: 0.01, 0.255, 0.5) --
alpha_vals <- c(0.001, 0.005, 0.01, 0.05)

## --- R0(alpha) for the autonomous model ---------------------
R0_alpha <- function(alpha) {
  Lambda/(alpha + mu) * ( beta1/nu1 + alpha*beta2/(mu*nu2) )
}

## --- DFE for constant alpha (autonomous model) --------------
S_DFE <- function(alpha) {
  S1_0 <- Lambda / (alpha + mu)
  S2_0 <- alpha * Lambda / (mu * (alpha + mu))
  c(S1_0 = S1_0, S2_0 = S2_0)
}

## --- Branching PGF map f(u) for given alpha -----------------
## u = (u1, u2)
f_map <- function(u, alpha) {
  S  <- S_DFE(alpha)
  S1_0 <- S["S1_0"]; S2_0 <- S["S2_0"]
  
  u1 <- u[1]; u2 <- u[2]
  
  denom1 <- nu1 + beta1 * (S1_0 + S2_0)
  denom2 <- nu2 + beta2 * (S1_0 + S2_0)
  
  f1 <- (beta1 * S1_0 * u1^2 + beta1 * S2_0 * u1 * u2 + nu1) / denom1
  f2 <- (beta2 * S2_0 * u2^2 + beta2 * S1_0 * u1 * u2 + nu2) / denom2
  
  # clamp to [0,1] to avoid tiny numeric excursions
  f1 <- min(max(f1, 0), 1)
  f2 <- min(max(f2, 0), 1)
  
  c(f1, f2)
}

## --- Compute u* for a given alpha (fixed point of f) --------
compute_u_star <- function(alpha, tol = 1e-8, max_iter = 2000) {
  R0 <- R0_alpha(alpha)
  if (R0 <= 1) return(c(1, 1))  # sure extinction
  
  u <- c(0, 0)
  for (k in 1:max_iter) {
    u_new <- f_map(u, alpha)
    
    if (any(is.na(u_new))) {
      # sécurité : si ça explose numériquement, on retourne extinction certaine
      return(c(1, 1))
    }
    
    if (max(abs(u_new - u)) < tol) {
      u <- u_new
      break
    }
    u <- u_new
  }
  u
}

## --- Precompute u*(alpha) WITHOUT rowwise/list complications --
u_mat <- sapply(alpha_vals, compute_u_star)  # matrix 2 x length(alpha_vals)

u_star_tbl <- tibble(
  alpha       = alpha_vals,
  u1          = u_mat[1, ],
  u2          = u_mat[2, ],
  R0          = R0_alpha(alpha_vals),
  alpha_label = paste0("alpha = ", alpha_vals)
)

print(u_star_tbl)
# tu dois maintenant voir u1 ET u2 non-NA, ~ (0.79, 0.66), etc.

## --- Range of initial infected ------------------------------
i10_max <- 20
i20_max <- 20

## 1) Seed in I1 only: i10 = 1..20, i20 = 0 --------------------
df_I1 <- expand.grid(
  alpha = alpha_vals,
  i10   = 1:i10_max   # >= 1
) %>%
  left_join(u_star_tbl, by = "alpha") %>%
  mutate(
    i20  = 0,
    Pext = (u1^i10) * (u2^i20),
    seed = "Seed in I1"
  )

## 2) Seed in I2 only: i20 = 1..20, i10 = 0 --------------------
df_I2 <- expand.grid(
  alpha = alpha_vals,
  i20   = 1:i20_max
) %>%
  left_join(u_star_tbl, by = "alpha") %>%
  mutate(
    i10  = 0,
    Pext = (u1^i10) * (u2^i20),
    seed = "Seed in I2"
  )

## Combine both scenarios --------------------------------------
df_all <- bind_rows(df_I1, df_I2) %>%
  mutate(
    seed = factor(seed, levels = c("Seed in I1", "Seed in I2")),
    # map alpha -> stress scenario, to reuse same colours as before
    scenario = case_when(
      alpha == 0.001 ~ "wet",
      alpha == 0.005 ~ "medium",
      alpha == 0.01  ~ "seasonal",
      alpha == 0.05  ~ "dry",
      TRUE           ~ "other"
    ),
    scenario = factor(scenario, levels = c("dry", "medium", "seasonal", "wet"))
  )

## --- Violin plot with colour by stress scenario --------------
p_violin <- ggplot(
  df_all,
  aes(x = factor(alpha), y = Pext, fill = scenario)
) +
  geom_violin(colour = "grey20", alpha = 0.7) +
  geom_jitter(width = 0.05, height = 0, alpha = 0.4, size = 1) +
  facet_wrap(~ seed) +
  scale_y_continuous(limits = c(0, 1)) +
  # labels de l’axe x en maths : alpha == 0.01, etc.
  scale_x_discrete(
    labels = function(x) parse(text = paste0("alpha==", x))
  ) +
  # même palette que pour tes autres figures (à adapter si besoin)
  scale_fill_manual(
    name   = "Scenario",
    values = c(
      "dry"      = "#D55E00",   # rouge/orange
      "medium"   = "#009E73",   # vert
      "seasonal" = "#0072B2",   # bleu
      "wet"      = "#CC79A7"    # violet
    )
  ) +
  labs(
    x = expression(alpha),
    y = expression(P[ext])
    # title = expression("Distribution of extinction probability" ~ P[ext] ~
    #                    " over initial conditions, for different " * alpha)
  ) +
  theme_bw()

print(p_violin)

ggsave("fig/violin_P_ext.pdf", p_violin, width = 9, height = 4, dpi = 300)
