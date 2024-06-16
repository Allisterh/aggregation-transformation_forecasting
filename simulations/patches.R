# Required packages
packages <- c("RandomFields", "ggplot2", "gsl", "R.utils", "BSDA", "doParallel", "foreach")
lapply(packages, require, ch = TRUE)
source("../utils/patched_SRs.R")
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

weights_d <- 1 / as.matrix(dist(coords(my_grid), diag = T, upper = T))
diag(weights_d) <- rep(0, nrow(weights_d))
weights_d <- weights_d / sum(weights_d)

# Observations
sigma_0 <- 1
lambda_0 <- 3
beta_0 <- 1

model <- RMstable(alpha = beta_0, scale = lambda_0, var = sigma_0)
plot(model)

simu <- RFsimulate(model, my_grid$x, my_grid$y, grid = TRUE, spConform = FALSE, n = n_tot, cPrintlevel = 0)

# --------------------------------- FORECASTS -------------------------------- #

lambdas <- c(1, 5)
betas <- c(.5, 2)

# Ideal
mu_ideal <- matrix(rep(0, d1 * d2), nrow = d1)
sigma_ideal <- sigma_0
lambda_ideal <- lambda_0
beta_ideal <- beta_0

mu_list <- list("ideal" = mu_ideal)
sigma_list <- list("ideal" = sigma_ideal)
lambda_list <- list("ideal" = lambda_ideal)
beta_list <- list("ideal" = beta_ideal)

# Misspecification of lambda
for (lbda in lambdas) {
    mu_f <- mu_ideal
    sigma_f <- sigma_ideal
    lambda_f <- lbda
    beta_f <- beta_ideal

    mu_temp <- list(mu_f)
    names(mu_temp) <- c(paste0("lambda_", lbda))
    sigma_temp <- list(sigma_f)
    names(sigma_temp) <- c(paste0("lambda_", lbda))
    lambda_temp <- list(lambda_f)
    names(lambda_temp) <- c(paste0("lambda_", lbda))
    beta_temp <- list(beta_f)
    names(beta_temp) <- c(paste0("lambda_", lbda))

    mu_list <- c(mu_list, mu_temp)
    sigma_list <- c(sigma_list, sigma_temp)
    lambda_list <- c(lambda_list, lambda_temp)
    beta_list <- c(beta_list, beta_temp)
}

# Misspecification of beta
for (beta in betas) {
    mu_f <- mu_ideal
    sigma_f <- sigma_ideal
    lambda_f <- lambda_ideal
    beta_f <- beta

    mu_temp <- list(mu_f)
    names(mu_temp) <- c(paste0("beta_", beta))
    sigma_temp <- list(sigma_f)
    names(sigma_temp) <- c(paste0("beta_", beta))
    lambda_temp <- list(lambda_f)
    names(lambda_temp) <- c(paste0("beta_", beta))
    beta_temp <- list(beta_f)
    names(beta_temp) <- c(paste0("beta_", beta))

    mu_list <- c(mu_list, mu_temp)
    sigma_list <- c(sigma_list, sigma_temp)
    lambda_list <- c(lambda_list, lambda_temp)
    beta_list <- c(beta_list, beta_temp)
}


# ---------------------------------------------------------------------------- #
#                       COMPUTATION OF THE SCORING RULES                       #
# ---------------------------------------------------------------------------- #

p_vs_list <- c(0.5, 1, 2) # List of orders for the Variogram score
p_pvs_list <- c(0.5, 1, 2) # List of power for p-Variation score
patch_list <- c(2, 3, 4, 5, 10) # List of patch size

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)


df_exp <- foreach(m = c("ideal", paste0('lambda_', lambdas), paste0('beta_', betas)), .combine = rbind, .packages = c('gsl', 'scoringRules', 'RandomFields')) %dopar% {
    df <- data.frame(model = character(), score_name = character(), score_value = numeric())
    mu <- mu_list[[m]]
    sigma <- sigma_list[[m]]
    lambda <- lambda_list[[m]]
    beta <- beta_list[[m]]

    # Variogram score
    for (p in p_vs_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                              score_name = rep(paste0("vs_", p), n_tot),
                              score_value = apply(simu, 3, function(x) { variogram_score(mu = mu, sigma = sigma, lambda = lambda, beta = beta, y = x, p = p, weights = weights_d, my_grid = my_grid) })
        )
        df <- rbind(df, df_temp)
    }

    # p-Variation score
    for (p in p_pvs_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                              score_name = rep(paste0("pvs_", p), n_tot),
                              score_value = apply(simu, 3, function(x) { p_variation_score(sigma = sigma, lambda = lambda, beta = beta, y = x, p = p, my_grid = my_grid) })
        )
        df <- rbind(df, df_temp)
    }

    # Aggregated CRPS
    df_temp <- data.frame(model = rep(m, n_tot),
                          score_name = rep("agg_crps", n_tot),
                          score_value = apply(simu, 3, function(x) { agg_CRPS_norm(mu = mu, sigma = sigma, y = x) })
    )
    df <- rbind(df, df_temp)

    # ----------------------------- Ensemble forecast ---------------------------- #
    forecast <- RMstable(alpha = beta, scale = lambda, var = sigma)
    ens_temp <- RFsimulate(forecast, my_grid$x, my_grid$y, grid = TRUE, spConform = FALSE, n = n_sample * n_tot, cPrintlevel = 0)

    ensemble <- ens_temp


    # Patched Energy score
    for (patch_size in patch_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                            score_name = rep(paste0("es_p", patch_size), n_tot),
                            score_value = sapply(1:n_tot, function(i) { patched_energy_score(ens = ensemble[,, ((i - 1) * n_sample + 1):(i * n_sample)], y = simu[,, i], p1 = patch_size, p2 = patch_size, my_grid = my_grid) })
      )
        df <- rbind(df, df_temp)
    }


    # Energy score
    df_temp <- data.frame(model = rep(m, n_tot),
                            score_name = rep("es", n_tot),
                            score_value = sapply(1:n_tot, function(i) { energy_score(ens = ensemble[,, ((i - 1) * n_sample + 1):(i * n_sample)], y = simu[,, i], my_grid = my_grid) })
      )
    df <- rbind(df, df_temp)

    df
}
stopCluster(cl)


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
              mat_empty, mat_empty, mat_empty)
names(p_val) <- list_sr
test_val <- p_val

for (sr in list_sr) {
    # Ignore agg. CRPS because same values for all forecasts 
    if (!(sr %in% c("agg_crps"))) {
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
}

# ----------------------- DM tests used in the article ----------------------- #

m = "ideal"

# Variogram score
for (p in p_vs_list) {
    print(paste0("vs_", p))
    print(rbind(p_val[[paste0("vs_", p)]][m,] < 0.05))
    print(test_val[[paste0("vs_", p)]][m,])
}

# p-Variation score
for (p in p_pvs_list) {
    print(paste0("pvs_", p))
    print(rbind(p_val[[paste0("pvs_", p)]][m,] < 0.05))
    print(test_val[[paste0("pvs_", p)]][m,])
}

# Patched Energy score
for (sr in c("agg_crps", paste0("es_p", patch_list), "es")) {
    print(sr)
    print(rbind(p_val[[sr]][m,] < 0.05))
    print(test_val[[sr]][m,])
}


# ---------------------------------------------------------------------------- #
#                                     PLOTS                                    #
# ---------------------------------------------------------------------------- #

list_models <- unique(df_exp$model)
list_sr <- unique(df_exp$score_name)

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
sr.labs <- c("Agg. CRPS", paste0("Patched ES (s=", patch_list, ")"), "Energy Score", paste0("Variogram Score (p=", p_vs_list, ")"), paste0("p-Variation Score (p=", p_pvs_list, ")"))
names(sr.labs) <- c("agg_crps", paste0("es_p", patch_list), "es", paste0("vs_", p_vs_list), paste0("pvs_", p_pvs_list))


# ------------------------------ Variogram score ----------------------------- #

dff <- subset(df_exp_rescaled, score_name %in% c(paste0("vs_", p_vs_list)))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf, alpha = rep(c(1, 1, 1, 1, 1), length(p_vs_list))) +
  geom_boxplot(outlier.shape = NA, coef = Inf, alpha = rep(c(1, 1, 1, 1, 1), length(p_vs_list))) +
  facet_wrap(~score_name, nrow = 1, labeller = labeller(score_name = sr.labs)) +
  scale_x_discrete(limits = c("ideal", paste0('lambda_', lambdas), paste0('beta_', betas))) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("tan1", "darkorange", "violet", "cornflowerblue", "darkblue"))
ggsave("../patches_vs.png", width = 5.5, height = 4)


# ----------------------------- p-Variation score ---------------------------- #

dff <- subset(df_exp_rescaled, score_name %in% c(paste0("pvs_", p_pvs_list)))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf) +
  geom_boxplot(outlier.shape = NA, coef = Inf) +
  facet_wrap(~score_name, nrow = 1, labeller = labeller(score_name = sr.labs)) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(limits = c("ideal", paste0('lambda_', lambdas), paste0('beta_', betas))) +
  scale_fill_manual(values = c("tan1", "darkorange", "violet", "cornflowerblue", "darkblue"))
ggsave("../patches_pvs.png", width = 5.5, height = 4)


# --------------------------- Patched Energy score --------------------------- #

dff <- subset(df_exp_rescaled, score_name %in% c("agg_crps", paste0("es_p", patch_list), "es"))
dff$score_name <- factor(dff$score_name,
                         levels = c("agg_crps", paste0("es_p", patch_list), "es"))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf) +
  geom_boxplot(outlier.shape = NA, coef = Inf) +
  facet_wrap(~score_name, nrow = 1, labeller = labeller(score_name = sr.labs)) +
  scale_x_discrete(limits = c("ideal", paste0('lambda_', lambdas), paste0('beta_', betas))) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("tan1", "darkorange", "violet", "cornflowerblue", "darkblue"))
ggsave("../patches_es.png", width = 11, height = 4)

