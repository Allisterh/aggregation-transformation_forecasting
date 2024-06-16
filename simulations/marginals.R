# Required packages
packages <- c("RandomFields", "ggplot2", "gsl", "R.utils", "BSDA", "latex2exp", "doParallel", "foreach")
lapply(packages, require, ch = TRUE)
source("../utils/agg_univariate_SRs.R")
source("../utils/general.R")

# RamdomFields options
RFoptions(seed = NA)
RFoptions(spConform = TRUE)

# General setup
n_obs <- 500
n_rep <- 10
n_tot <- n_obs * n_rep
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

# Ideal
mu_ideal <- matrix(rep(0, d1 * d2), nrow = d1)
sigma_ideal <- sigma_0
lambda_ideal <- lambda_0
beta_ideal <- beta_0

# Biased
mu_bias <- 0.255
sigma_bias <- sigma_ideal
lambda_bias <- lambda_0
beta_bias <- beta_0

# Underdispersed
mu_under <- mu_ideal
sigma_under <- 2 / 3 * sigma_0
lambda_under <- lambda_0
beta_under <- beta_0

# Overdispersed
mu_over <- mu_ideal
sigma_over <- 1.4 * sigma_0
lambda_over <- lambda_0
beta_over <- beta_0

mu_list <- list("ideal" = mu_ideal, "biased" = mu_bias, "underdispersed" = mu_under, "overdispersed" = mu_over)
sigma_list <- list("ideal" = sigma_ideal, "biased" = sigma_bias, "underdispersed" = sigma_under, "overdispersed" = sigma_over)
lambda_list <- list("ideal" = lambda_ideal, "biased" = lambda_bias, "underdispersed" = lambda_under, "overdispersed" = lambda_over)
beta_list <- list("ideal" = beta_ideal, "biased" = beta_bias, "underdispersed" = beta_under, "overdispersed" = beta_over)


# ---------------------------------------------------------------------------- #
#                       COMPUTATION OF THE SCORING RULES                       #
# ---------------------------------------------------------------------------- #

t_list <- c(0.5, 1) # List of thresholds for the Brier score
alpha_list <- c(.5, .75, .95) # List of quantiles for the Quantile score

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)


df_exp <- foreach(m = c("ideal", "biased", "underdispersed", "overdispersed"), .combine = rbind, .packages = c('gsl', 'scoringRules', 'RandomFields')) %dopar% {
    df <- data.frame(model = character(), score_name = character(), score_value = numeric())
    mu <- mu_list[[m]]
    sigma <- sigma_list[[m]]
    lambda <- lambda_list[[m]]
    beta <- beta_list[[m]]

    # Aggregated CRPS
    df_temp <- data.frame(model = rep(m, n_tot),
                          score_name = rep("agg_crps", n_tot),
                          score_value = apply(simu, 3, function(x) { agg_CRPS_norm(mu = mu, sigma = sigma, y = x) })
    )
    df <- rbind(df, df_temp)


    # Aggregated Squared Error
    df_temp <- data.frame(model = rep(m, n_tot),
                          score_name = rep("agg_se", n_tot),
                          score_value = apply(simu, 3, function(x) { agg_SE(mu = mu, y = x) })
    )
    df <- rbind(df, df_temp)


    # Aggregated Dawid-Sebastiani Score
    df_temp <- data.frame(model = rep(m, n_tot),
                          score_name = rep("agg_dss", n_tot),
                          score_value = apply(simu, 3, function(x) { agg_DSS(mu = mu, sigma = sigma, y = x) })
    )
    df <- rbind(df, df_temp)

    # Aggregated Brier Score
    for (t in t_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                              score_name = rep(paste0("agg_bs_t", t), n_tot),
                              score_value = apply(simu, 3, function(x) { agg_BS_norm(mu = mu, sigma = sigma, y = x, threshold = t) })
      )
        df <- rbind(df, df_temp)
    }


    # Aggregated Quantile Score
    for (alpha in alpha_list) {
        df_temp <- data.frame(model = rep(m, n_tot),
                              score_name = rep(paste0("agg_qs_a", alpha), n_tot),
                              score_value = apply(simu, 3, function(x) { agg_QS_norm(mu = mu, sigma = sigma, y = x, alpha = alpha) })
        )
        df <- rbind(df, df_temp)
    }
    df

}
stopCluster(cl)

# ------------------------------- Student model ------------------------------ #
m = "student"
deg_freedom <- 5
tau <- 0.745 / sqrt(deg_freedom / (deg_freedom - 2))
sigma <- tau * sqrt(deg_freedom / (deg_freedom - 2))
mu <- mu_ideal

# Aggregated CRPS
df_temp <- data.frame(model = rep(m, n_tot),
                      score_name = rep("agg_crps", n_tot),
                      score_value = apply(simu, 3, function(x) { agg_CRPS_t(df = deg_freedom, y = x, tau = tau) })
)
df_exp <- rbind(df_exp, df_temp)


# Aggregated Squared Error
df_temp <- data.frame(model = rep(m, n_tot),
                      score_name = rep("agg_se", n_tot),
                      score_value = apply(simu, 3, function(x) { agg_SE(mu = mu, y = x) })
)
df_exp <- rbind(df_exp, df_temp)


# Aggregated Dawid-Sebastiani Score
df_temp <- data.frame(model = rep(m, n_tot),
                      score_name = rep("agg_dss", n_tot),
                      score_value = apply(simu, 3, function(x) { agg_DSS(mu = mu, sigma = sigma, y = x) })
)
df_exp <- rbind(df_exp, df_temp)

# Aggregated Brier Score
for (t in t_list) {
    df_temp <- data.frame(model = rep(m, n_tot),
                          score_name = rep(paste0("agg_bs_t", t), n_tot),
                          score_value = apply(simu, 3, function(x) { agg_BS_t(df = deg_freedom, y = x, threshold = t, tau = tau) })
  )
    df_exp <- rbind(df_exp, df_temp)
}


# Aggregated Quantile Score
for (alpha in alpha_list) {
    df_temp <- data.frame(model = rep(m, n_tot),
                          score_name = rep(paste0("agg_qs_a", alpha), n_tot),
                          score_value = apply(simu, 3, function(x) { agg_QS_t(df = deg_freedom, y = x, alpha = alpha) })
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
p_val <- list(mat_empty, mat_empty, mat_empty, mat_empty, mat_empty, mat_empty, mat_empty, mat_empty)
names(p_val) <- list_sr
test_val <- p_val

for (sr in list_sr) {
    # Ignore scoring rules leading to same values for multiple forecasts 
    if (!(sr %in% c("agg_se", "agg_qs_a0.5"))) { 
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

# Aggregated CRPS
for (i in 1:length(list_models)) {
    m = list_models[i]
    print(m)
    print(rbind(p_val$agg_crps[i, - i] < 0.05))
    print(test_val$agg_crps[i, - i])
}

# Student - Quantile score
m = "student"
print(m)
print(rbind(p_val$agg_qs_a0.75[m,] < 0.05))
print(test_val$agg_qs_a0.75[m,])

# Biased forecast - Quantile score/Brier score/Dawid-Sebastiani score
m = "biased"
print(m)
print(rbind(p_val$agg_qs_a0.95[m,] < 0.05))
print(test_val$agg_qs_a0.95[m,])

print(rbind(p_val$agg_bs_t1[m,] < 0.05))
print(test_val$agg_bs_t1[m,])

print(rbind(p_val$agg_dss[m,] < 0.05))
print(test_val$agg_dss[m,])

# Overdispersed forecast - Brier score
m = "overdispersed"
print(m)
print(rbind(p_val$agg_bs_t0.5[m,] < 0.05))
print(test_val$agg_bs_t0.5[m,])

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
sr.labs <- c("Aggregated CRPS", "Aggregated SE", "Aggregated DSS", paste0("Aggregated BS (t=", t_list, ")"), paste0("Aggregated QS (alpha=", alpha_list, ")"))
names(sr.labs) <- c("agg_crps", "agg_se", "agg_dss", paste0("agg_bs_t", t_list), paste0("agg_qs_a", alpha_list))


# ----------------------------------- CRPS ----------------------------------- #

ggplot(subset(df_exp_rescaled, score_name == "agg_crps"), aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf, alpha = c(.25, 1, 1, 1, 1)) +
  geom_boxplot(outlier.shape = NA, coef = Inf, alpha = c(.25, 1, 1, 1, 1)) +
  facet_wrap(~score_name, scale = "free", nrow = 1, labeller = labeller(score_name = sr.labs)) +
  scale_x_discrete(limits = c("ideal", "biased", "underdispersed", "overdispersed", "student")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("darkorange", "violet", "darkblue", "springgreen4", "cornflowerblue"))
ggsave("../marginals_crps.png", width = 3, height = 4)


# -------------------------------- Brier score ------------------------------- #

dff <- subset(df_exp_rescaled, score_name %in% paste0("agg_bs_t", t_list))
dff$score_name <- factor(dff$score_name,
                         levels = paste0("agg_bs_t", t_list),
                         labels = TeX(sprintf(r'(Aggregated BS ($t = %f$))', t_list)))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf, alpha = rep(c(.25, 1, 1, 1, 1), 2)) +
  geom_boxplot(outlier.shape = NA, coef = Inf, alpha = rep(c(.25, 1, 1, 1, 1), 2)) +
  facet_wrap(~score_name, nrow = 1, labeller = label_parsed) +
  scale_x_discrete(limits = c("ideal", "biased", "underdispersed", "overdispersed", "student")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("darkorange", "violet", "darkblue", "springgreen4", "cornflowerblue"))
ggsave("../marginals_bs.png", width = 5, height = 4)


# ------------------------------ Quantile score ------------------------------ #

dff <- subset(df_exp_rescaled, score_name %in% paste("agg_qs_a", alpha_list, sep = ""))
dff$score_name <- factor(dff$score_name,
                         levels = paste0("agg_qs_a", alpha_list),
                         labels = TeX(sprintf(r'(Aggregated QS ($\alpha = %f$))', alpha_list)))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf, alpha = rep(c(.25, 1, 1, 1, 1), length(alpha_list))) +
  geom_boxplot(outlier.shape = NA, coef = Inf, alpha = rep(c(.25, 1, 1, 1, 1), length(alpha_list))) +
  facet_wrap(~score_name, scale = "free", nrow = 1, labeller = label_parsed) +
  scale_x_discrete(limits = c("ideal", "biased", "underdispersed", "overdispersed", "student")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("darkorange", "violet", "darkblue", "springgreen4", "cornflowerblue"))
ggsave("../marginals_qs.png", width = 7.5, height = 4)


# --------------------------- Moments-based scores --------------------------- #

dff <- subset(df_exp_rescaled, score_name %in% c("agg_se", "agg_dss"))
dff$score_name <- factor(dff$score_name,
                         levels = c("agg_se", "agg_dss"))

ggplot(dff, aes(x = model, y = score_value, fill = model)) +
  stat_boxplot(geom = 'errorbar', coef = Inf, alpha = rep(c(.25, 1, 1, 1, 1), 2)) +
  geom_boxplot(outlier.shape = NA, coef = Inf, alpha = rep(c(.25, 1, 1, 1, 1), 2)) +
  facet_wrap(~score_name, scale = "free", nrow = 1, labeller = labeller(score_name = sr.labs)) +
  scale_x_discrete(limits = c("ideal", "biased", "underdispersed", "overdispersed", "student")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("darkorange", "violet", "darkblue", "springgreen4", "cornflowerblue"))
ggsave("../marginals_moments.png", width = 5.5, height = 4)
