# ---------------------------------------------------------------------------- #
#                      AGGREGATED UNIVARIATE SCORING RULES                     #
# ---------------------------------------------------------------------------- #

# ---------------------------- REQUIRED LIBRARIES ---------------------------- #

require(scoringRules)
require(extraDistr)

# --------------------------------- FUNCTIONS -------------------------------- #
agg_CRPS_norm <- function(mu, sigma, y) {
    return(mean(crps_norm(y = y, mean = mu, sd = sigma)))
}

agg_CRPS_t <- function(df, y, tau = 1) {
    return(mean(crps_t(y = y, df = df, scale = tau)))
}

agg_SE <- function(mu, y) {
    return(mean((y - mu) ^ 2))
}
agg_DSS <- function(mu, sigma, y) {
    return(mean(2 * log(sigma) + ((y - mu) / sigma) ^ 2))
}

agg_BS_norm <- function(mu, sigma, y, threshold) {
    return(mean((pnorm(q = t, mean = mu, sd = sigma) - (y <= t)) ^ 2))
}

agg_BS_t <- function(df, y, threshold, tau = 1) {
    return(mean((plst(q = t, df = df, sigma = tau) - (y <= t)) ^ 2))
}

agg_QS_norm <- function(mu, sigma, y, alpha) {
    q <- qnorm(p = alpha, mean = mu, sd = sigma)
    return(mean(((y <= q) - alpha) * (q - y)))
}

agg_QS_t <- function(df, y, alpha, tau = 1) {
    q <- qlst(p = alpha, df = df, sigma = tau)
    return(mean(((y <= q) - alpha) * (q - y)))
}