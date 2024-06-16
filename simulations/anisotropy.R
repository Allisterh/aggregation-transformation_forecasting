# Required packages
packages <- c("RandomFields", "ggplot2", "gsl", "R.utils", "BSDA", "latex2exp", "doParallel", "foreach", "mev")
lapply(packages, require, ch = TRUE)
source("../utils/anisotropy_SRs.R")
source("../utils/general.R")

# RamdomFields options
RFoptions(seed = NA)
RFoptions(spConform = TRUE)

# General setup
n_obs <- 500
n_rep <- 10
n_tot <- n_obs * n_rep
n_sample <- 100
d1 <- 20
d2 <- 20
my_grid <- list(x = seq(0, d1 - 1, by = 1), y = seq(0, d2 - 1, by = 1), ind.x = 1:d1, ind.y = 1:d2)

# Observations
sigma_0 <- 1
lambda_0 <- 3
beta_0 <- 1

angle_0 <- pi / 4
ratio_0 <- 1 / 2
A_0 <- RMangle(angle = angle_0, ratio = ratio_0)

model <- RMstable(alpha = beta_0, scale = lambda_0, var = sigma_0, Aniso = A_0)
plot(model, dim = 2)
image(model)

simu <- RFsimulate(model, my_grid$x, my_grid$y, grid = TRUE, spConform = FALSE, n = n_obs * n_rep, cPrintlevel = 0)

# Weights for the Variogram score
weights_d <- 1 / as.matrix(dist(coords(my_grid), diag = T, upper = T))
diag(weights_d) <- rep(0, nrow(weights_d))
weights_d <- weights_d / sum(weights_d)

weights_dA <- 1 / mev::distg(loc = coords(my_grid), scale = 1 / ratio_0, rho = -angle_0)
diag(weights_dA) <- rep(0, nrow(weights_dA))
weights_dA <- weights_dA / sum(weights_dA)

# --------------------------------- FORECASTS -------------------------------- #

# Ideal
mu_ideal <- matrix(rep(0, d1 * d2), nrow = d1)
sigma_ideal <- sigma_0
lambda_ideal <- lambda_0
beta_ideal <- beta_0
A_ideal <- A_0

# Angle 1 : small-angle
mu_angle1 <- mu_ideal
sigma_angle1 <- sigma_ideal
lambda_angle1 <- lambda_0
beta_angle1 <- beta_0
A_angle1 <- RMangle(angle = 0, ratio = ratio_0)

# Angle 2 : large-angle
mu_angle2 <- mu_ideal
sigma_angle2 <- sigma_ideal
lambda_angle2 <- lambda_0
beta_angle2 <- beta_0
A_angle2 <- RMangle(angle = pi / 2, ratio = ratio_0)

# Ratio 1 : isotropic
mu_ratio1 <- mu_ideal
sigma_ratio1 <- sigma_ideal
lambda_ratio1 <- lambda_0
beta_ratio1 <- beta_0
A_ratio1 <- RMangle(angle = angle_0, ratio = 1)

# Ratio 2 : over-anisotropic
mu_ratio2 <- mu_ideal
sigma_ratio2 <- sigma_ideal
lambda_ratio2 <- lambda_0
beta_ratio2 <- beta_0
A_ratio2 <- RMangle(angle = angle_0, ratio = 1 / 3)

A_list <- list("ideal" = A_ideal, "angle1" = A_angle1, "angle2" = A_angle2, "ratio1" = A_ratio1, "ratio2" = A_ratio2)


# ---------------------------------------------------------------------------- #
#                       COMPUTATION OF THE SCORING RULES                       #
# ---------------------------------------------------------------------------- #

p_list <- c(.5) # List of orders for the Variogram score
h_list <- c(1, 2, 3, 4, 5) # List of scale for the Anisotropic score

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

df_exp <- foreach(m = c("ideal", "angle1", "angle2", "ratio1", "ratio2"), .combine = rbind, .packages = c('gsl', 'scoringRules', 'RandomFields')) %dopar% {
    df <- data.frame(model = character(), score_name = character(), score_value = numeric())
    mu <- mu_ideal
    sigma <- sigma_0
    lambda <- lambda_0
    beta <- beta_0
    A <- A_list[[m]]

    # Variogram score (naive weights)
    for (p in p_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                            score_name = rep(paste0("vs_", p), n_tot),
                            score_value = apply(simu, 3, function(x) { variogram_score(mu = mu, sigma = sigma, lambda = lambda, beta = beta, y = x, p = p, weights = weights_d, my_grid = my_grid, A = A) })
        )
        df <- rbind(df, df_temp)
    }

    # Variogram score (informed weights)
    for (p in p_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                            score_name = rep(paste0("vs-i_", p), n_tot),
                            score_value = apply(simu, 3, function(x) { variogram_score(mu = mu, sigma = sigma, lambda = lambda, beta = beta, y = x, p = p, weights = weights_dA, my_grid = my_grid, A = A) })
        )
        df <- rbind(df, df_temp)
    }

    # ----------------------------- Ensemble forecast ---------------------------- #
    forecast <- RMstable(alpha = beta, scale = lambda, var = sigma, Aniso = A)
    ens_temp <- RFsimulate(forecast, my_grid$x, my_grid$y, grid = TRUE, spConform = FALSE, n = n_sample * n_tot, cPrintlevel = 0)
    ensemble <- ens_temp

    # Anisotropic score (old = preprint of Allen et al.)
    for (h in h_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                            score_name = rep(paste0("as_old_h", h), n_tot),
                            score_value = sapply(1:n_tot, function(i) { SE_Tiso_old(ens = ensemble[,, ((i - 1) * n_sample + 1):(i * n_sample)], y = simu[,, i], h = h) })
      )
        df <- rbind(df, df_temp)
    }

    # Anisotropic score (a1 = axes corresponding to x and y)
    for (h in h_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                            score_name = rep(paste0("as_a1_h", h), n_tot),
                            score_value = sapply(1:n_tot, function(i) { SE_Tiso_a1(ens = ensemble[,, ((i - 1) * n_sample + 1):(i * n_sample)], y = simu[,, i], h = h) })
      )
        df <- rbind(df, df_temp)
    }

    # Anisotropic score (a2 = axes corresponding to bisectors)
    for (h in h_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                            score_name = rep(paste0("as_a2_h", h), n_tot),
                            score_value = sapply(1:n_tot, function(i) { SE_Tiso_a2(ens = ensemble[,, ((i - 1) * n_sample + 1):(i * n_sample)], y = simu[,, i], h = h) })
      )
        df <- rbind(df, df_temp)
    }
    df
}
stopCluster(cl)


# Aggregated anisotropic score (w_h=1/h)

# List of models and scoring rules
list_models <- unique(df_exp$model)
list_sr <- unique(df_exp$score_name)

for (m in list_models) {
    anisotropic_wh <- rep(0, n_tot)
    for (h in h_list) {
        anisotropic_wh = anisotropic_wh + 1 / h * subset(df_exp, model == m & score_name == paste0("as_a2_h", h))$score_value
    }
    df_temp <- data.frame(model = rep(m, n_tot),
                            score_name = rep("as_a2_wh", n_tot),
                            score_value = anisotropic_wh
      )
    df_exp <- rbind(df_exp, df_temp)
}

# ---------------------------------------------------------------------------- #
#                             DIEBOLD-MARIANO TEST                             #
# ---------------------------------------------------------------------------- #

# List of models and scoring rules
list_models <- unique(df_exp$model)
list_sr <- unique(df_exp$score_name)

# Initialize list of p-values and test statistics 
mat_empty <- matrix(, nrow = length(list_models), ncol = length(list_models))
row.names(mat_empty) <- list_models
colnames(mat_empty) <- list_models
p_val <- list(mat_empty, mat_empty, mat_empty, mat_empty, mat_empty,
              mat_empty, mat_empty, mat_empty, mat_empty, mat_empty,
              mat_empty, mat_empty, mat_empty, mat_empty, mat_empty,
              mat_empty, mat_empty, mat_empty)
names(p_val) <- list_sr
test_val <- p_val

for (sr in list_sr) {
    for (m_ref in list_models) {
        s2 <- subset(df_exp, model == m_ref & score_name == sr)$score_value
        for (m in list_models) {
            if (m != m_ref) {
                s1 <- subset(df_exp, model == m & score_name == sr)$score_value
                x <- (s1 - s2) / sd(s1 - s2)

                z_test = z.test(x, sigma.x = 1)
                p_val[[sr]][m_ref, m] = z_test$p.value
                test_val[[sr]][m_ref, m] = z_test$estimate
            }
        }
    }
}

# ----------------------- DM tests used in the article ----------------------- #

m = "ideal"

# Variogram score (standard)
for (p in p_list) {
    print(paste0("vs_", p))
    print(rbind(p_val[[paste0("vs_", p)]][m,] < 0.05))
    print(test_val[[paste0("vs_", p)]][m,])
}

# Variogram score (informed)
for (p in p_list) {
    print(paste0("vs-i_", p))
    print(rbind(p_val[[paste0("vs-i_", p)]][m,] < 0.05))
    print(test_val[[paste0("vs-i_", p)]][m,])
}

# Anisotropic score
for (h in h_list) {
    print(paste0("as_a2_h", h))
    print(rbind(p_val[[paste0("as_a2_h", h)]][m,] < 0.05))
    print(test_val[[paste0("as_a2_h", h)]][m,])
}

# Aggregated anisotropic score
print("as_a2_wh")
print(rbind(p_val[["as_a2_wh"]][m,] < 0.05))
print(test_val[["as_a2_wh"]][m,])


# ---------------------------------------------------------------------------- #
#                                     PLOTS                                    #
# ---------------------------------------------------------------------------- #

stride <- n_obs
chunk_size <- n_obs

# Averaging over repetitions
df_exp_avg <- data.frame(model = character(), score_name = character(), score_value = numeric())
for (sr in list_sr) {
    for (m in list_models) {
        x <- subset(df_exp, model == m & score_name == sr)$score_value
        n_chunk <- length(seq(1, length(x) - chunk_size + 1, stride))
        df_temp <- data.frame(model = rep(m, n_chunk),
                              score_name = rep(sr, n_chunk),
                              score_value = mean_chunk(x, chunk_size = chunk_size, stride = stride))
        df_exp_avg <- rbind(df_exp_avg, df_temp)
    }
}

# Rescaling scoring rules with respect to the mean of the ideal forecast
df_exp_rescaled <- data.frame(model = character(), score_name = character(), score_value = numeric())
for (sr in list_sr) {
    x_ideal <- mean(subset(df_exp_avg, model == "ideal" & score_name == sr)$score_value)
    for (m in list_models) {
        x <- subset(df_exp_avg, model == m & score_name == sr)$score_value / x_ideal
        df_temp <- data.frame(model = rep(m, length(x)),
                              score_name = rep(sr, length(x)),
                              score_value = x)
        df_exp_rescaled <- rbind(df_exp_rescaled, df_temp)
    }
}

# Labels
sr.labs <- c(paste0("Variogram Score (p=", p_list, ")"), paste0("Informed VS (p=", p_list, ")"), paste0("Anisotropic score (h=", h_list, ")"), paste0("Anisotropic score (h=", h_list, ")"), paste0("Anisotropic score (h=", h_list, ")"), "Anisotropic score (w_h=1/h)")
names(sr.labs) <- c(paste0("vs_", p_list), paste0("vs-i_", p_list), paste0("as_old_h", h_list), paste0("as_a1_h", h_list), paste0("as_a2_h", h_list), "as_a2_wh")


# ------------------------------ Variogram score ----------------------------- #

dff <- subset(df_exp_rescaled, score_name %in% c("vs_0.5", "vs-i_0.5"))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf, alpha = rep(c(1, 1, 1, 1, 1), 2)) +
  geom_boxplot(outlier.shape = NA, coef = Inf, alpha = rep(c(1, 1, 1, 1, 1), 2)) +
  facet_wrap(~score_name, nrow = 1, labeller = labeller(score_name = sr.labs)) +
  scale_x_discrete(limits = c("ideal", "angle1", "angle2", "ratio1", "ratio2")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("cornflowerblue", "darkblue", "violet", "tan1", "darkorange"))
ggsave("../anisotropy_vs.png", width = 6, height = 4)


# ------------ Anisotropic score (old) [not shown in the article] ------------ #
dff <- subset(df_exp_rescaled, score_name %in% paste0("as_old_h", h_list))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf) +
  geom_boxplot(outlier.shape = NA, coef = Inf) +
  facet_wrap(~score_name, scale = "free", nrow = 1, labeller = labeller(score_name = sr.labs)) +
  scale_x_discrete(limits = c("ideal", "angle1", "angle2", "ratio1", "ratio2")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("cornflowerblue", "darkblue", "violet", "tan1", "darkorange"))
ggsave("../anisotropy_as_old.png", width = 15, height = 4)


# ------------- Anisotropic score (a1) [not shown in the article] ------------ #

dff <- subset(df_exp_rescaled, score_name %in% paste0("as_a1_h", h_list))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf) +
  geom_boxplot(outlier.shape = NA, coef = Inf) +
  facet_wrap(~score_name, scale = "free", nrow = 1, labeller = labeller(score_name = sr.labs)) +
  scale_x_discrete(limits = c("ideal", "angle1", "angle2", "ratio1", "ratio2")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("cornflowerblue", "darkblue", "violet", "tan1", "darkorange"))
ggsave("../anisotropy_as_a1.png", width = 15, height = 4)


# -------------------------- Anisotropic score (a2) -------------------------- #
dff <- subset(df_exp_rescaled, score_name %in% c(paste0("as_a2_h", h_list), "as_a2_wh"))
dff$score_name <- factor(dff$score_name,
                         levels = c(paste0("as_a2_h", h_list), "as_a2_wh"),
                         labels = TeX(c(sprintf(r'(AS ($h = %f$))', h_list), r'(AS ($w_h = 1/h$))')))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf) +
  geom_boxplot(outlier.shape = NA, coef = Inf) +
  facet_wrap(~score_name, scale = "free", nrow = 1, labeller = label_parsed) +
  scale_x_discrete(limits = c("ideal", "angle1", "angle2", "ratio1", "ratio2")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("cornflowerblue", "darkblue", "violet", "tan1", "darkorange"))
ggsave("../anisotropy_as_a2.png", width = 11, height = 4)





