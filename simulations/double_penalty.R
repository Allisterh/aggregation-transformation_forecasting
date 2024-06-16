# Required packages
packages <- c("RandomFields", "ggplot2", "gsl", "R.utils", "BSDA", "doParallel", "foreach")
lapply(packages, require, ch = TRUE)
source("../utils/double_penalty_SRs.R")
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

model <- RMstable(alpha = beta_0, scale = lambda_0, var = sigma_0)
plot(model)

simu <- RFsimulate(model, my_grid$x, my_grid$y, grid = TRUE, spConform = FALSE, n = n_tot, cPrintlevel = 0)


# --------------------------------- FORECASTS -------------------------------- #

eps_list <- c(.1, .25, .5) # List of range of the noises

# Ideal
mu_ideal <- matrix(rep(0, d1 * d2), nrow = d1)
sigma_ideal <- sigma_0 * matrix(rep(1, d1 * d2), nrow = d1)
lambda_ideal <- lambda_0
beta_ideal <- beta_0

mu_list <- list("ideal" = mu_ideal)
sigma_list <- list("ideal" = sigma_ideal)
lambda_list <- list("ideal" = lambda_ideal)
beta_list <- list("ideal" = beta_ideal)

for (eps in eps_list) {
    # Additive noise
    mu_bias <- matrix(runif(n = d1 * d2, min = -eps, max = eps), nrow = d1)
    sigma_bias <- sigma_ideal
    lambda_bias <- lambda_0
    beta_bias <- beta_0

    # Multiplicative noise
    mu_variance <- mu_ideal
    sigma_variance <- matrix(runif(n = d1 * d2, min = 1 - eps, max = 1 + eps), nrow = d1)
    lambda_variance <- lambda_0
    beta_variance <- beta_0

    # Storing parameters
    mu_temp <- list(mu_bias, mu_variance)
    names(mu_temp) <- c(paste0("bias_eps", eps), paste0("variance_eps", eps))
    mu_list <- c(mu_list, mu_temp)
    sigma_temp <- list(sigma_bias, sigma_variance)
    names(sigma_temp) <- c(paste0("bias_eps", eps), paste0("variance_eps", eps))
    sigma_list <- c(sigma_list, sigma_temp)
    lambda_temp <- list(lambda_bias, lambda_variance)
    names(lambda_temp) <- c(paste0("bias_eps", eps), paste0("variance_eps", eps))
    lambda_list <- c(lambda_list, lambda_temp)
    beta_temp <- list(beta_bias, beta_variance)
    names(beta_temp) <- c(paste0("bias_eps", eps), paste0("variance_eps", eps))
    beta_list <- c(beta_list, beta_temp)
}


# ---------------------------------------------------------------------------- #
#                       COMPUTATION OF THE SCORING RULES                       #
# ---------------------------------------------------------------------------- #

patch_list <- c(2, 3, 5, 10) # List of patch sizes
t_list <- c(1) # List of threshold for the Brier score and Squared error of FTE

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

df_exp <- foreach(m = c("ideal", paste0("bias_eps", eps_list), paste0("variance_eps", eps_list)), .combine = rbind, .packages = c('gsl', 'scoringRules', 'RandomFields')) %dopar% {
    df <- data.frame(model = character(), score_name = character(), score_value = numeric())
    mu <- mu_list[[m]]
    sigma <- sigma_list[[m]]
    lambda <- lambda_list[[m]]
    beta <- beta_list[[m]]

    # Aggregated CRPS
    df_temp <- data.frame(model = rep(m, n_obs),
                        score_name = rep("agg_crps", n_obs),
                        score_value = apply(simu, 3, function(x) { agg_CRPS_norm(mu = mu, sigma = sigma, y = x) })
  )
    df <- rbind(df, df_temp)

    # Aggregated CRPS of spatial mean
    for (patch_size in patch_list) {
        df_temp <- data.frame(model = rep(m, n_obs),
                        score_name = rep(paste0("agg_crps_sm_p", patch_size), n_obs),
                        score_value = apply(simu, 3, function(x) { agg_CRPS_mean(mu = mu, sigma = sigma, lambda = lambda, beta = beta, y = x, p1 = patch_size, p2 = patch_size, my_grid = my_grid) })
  )
        df <- rbind(df, df_temp)
    }

    # Aggregated Brier score
    for (t in t_list) {
        df_temp <- data.frame(model = rep(m, n_obs),
                          score_name = rep(paste0("agg_bs_t", t), n_obs),
                          score_value = apply(simu, 3, function(x) { agg_BS_norm(mu = mu, sigma = sigma, y = x, threshold = t) })
    )
        df <- rbind(df, df_temp)
    }

    # ----------------------------- Ensemble forecast ---------------------------- #
    forecast <- RMstable(alpha = beta, scale = lambda, var = sigma_0)
    ens_temp <- RFsimulate(forecast, my_grid$x, my_grid$y, grid = TRUE, spConform = FALSE, n = n_sample * n_tot, cPrintlevel = 0)
    ensemble <- array(, dim = dim(ens_temp))
    for (i in 1:dim(ens_temp)[3]) {
        ensemble[,, i] <- sigma * ens_temp[,, i] + mu
    }

    # Aggregated Brier score of FTE
    for (patch_size in patch_list) {
        for (t in t_list) {
            df_temp <- data.frame(model = rep(m, n_tot),
                            score_name = rep(paste0("agg_bs_fte_t", t, "_p", patch_size), n_tot),
                            score_value = sapply(1:n_tot, function(i) { agg_BS_fte(ens = ensemble[,, ((i - 1) * n_sample + 1):(i * n_sample)], y = simu[,, i], p1 = patch_size, p2 = patch_size, my_grid = my_grid, t = t) })
      )
            df <- rbind(df, df_temp)
        }
    }
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
              mat_empty, mat_empty, mat_empty, mat_empty, mat_empty)
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

# Aggregated CRPS and Aggregated CRPS of spatial mean
for (sr in c("agg_crps", paste0("agg_crps_sm_p", patch_list))) {
    print(sr)
    print(rbind(p_val[[sr]][m,] < 0.05))
    print(test_val[[sr]][m,])
}

# Aggregated Brier score and Aggregated Squared error of FTE
for (sr in c(paste0("agg_bs_t", t_list), paste0("agg_bs_fte_t", t_list, "_p", patch_list))) {
    print(sr)
    print(rbind(p_val[[sr]][m,] < 0.05))
    print(test_val[[sr]][m,])
}


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
sr.labs <- c("Aggregated CRPS", paste0("Aggregated CRPS of\n spatial mean (s=", patch_list, ")"), paste0("Aggregated BS (t=", t_list, ")"), paste0("Aggregated SE of\n FTE (t=", t_list, ", s=", patch_list, ")"))
names(sr.labs) <- c("agg_crps", paste0("agg_crps_sm_p", patch_list), paste0("agg_bs_t", t_list), paste0("agg_bs_fte_t", t_list, "_p", patch_list))


# -------------------- Continuous Ranked Probability Score ------------------- #

dff <- subset(df_exp_rescaled, score_name %in% c(paste0("agg_crps_sm_p", patch_list), "agg_crps"))
dff$score_name <- factor(dff$score_name,
                         levels = c("agg_crps", paste0("agg_crps_sm_p", patch_list)))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf, alpha = rep(rep(1, 7), length(patch_list) + 1)) +
  geom_boxplot(outlier.shape = NA, coef = Inf, alpha = rep(rep(1, 7), length(patch_list) + 1)) +
  facet_wrap(~score_name, nrow = 1, labeller = labeller(score_name = sr.labs)) +
  scale_x_discrete(limits = c("ideal", paste0("bias_eps", eps_list), paste0("variance_eps", eps_list))) +
  theme(legend.position = "none",
        strip.text = element_text(size = 15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("cornflowerblue", "blue", "navy", "violet", "tan1", "darkorange", "darkorange4"))
ggsave("../double_penalty_crps.png", width = 11, height = 4)


# -------------------------------- Brier score ------------------------------- #
dff <- subset(df_exp_rescaled, score_name %in% c(paste0("agg_bs_fte_t", t_list, "_p", patch_list), paste0("agg_bs_t", t_list)))
dff$score_name <- factor(dff$score_name,
                         levels = c(paste0("agg_bs_t", t_list), paste0("agg_bs_fte_t", t_list, "_p", patch_list)))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf) +
geom_boxplot(outlier.shape = NA, coef = Inf) +
facet_wrap(~score_name, nrow = 1, labeller = labeller(score_name = sr.labs)) +
  scale_x_discrete(limits = c("ideal", paste0("bias_eps", eps_list), paste0("variance_eps", eps_list))) +
  theme(legend.position = "none",
        strip.text = element_text(size = 15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("cornflowerblue", "blue", "navy", "violet", "tan1", "darkorange", "darkorange4"))
ggsave("../double_penalty_bs.png", width = 11, height = 4)

