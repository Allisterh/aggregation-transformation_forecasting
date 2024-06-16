# ---------------------------------------------------------------------------- #
#                         DOUBLE-PENALTY SCORING RULES                         #
# ---------------------------------------------------------------------------- #

# -------------------- AGGREGATED UNIVARIATE SCORING RULES ------------------- #

agg_CRPS_norm <- function(mu, sigma, y) {
    return(mean(crps_norm(y = y, mean = mu, sd = sigma)))
}

agg_SE <- function(mu, y) {
    return(mean((y - mu) ^ 2))
}

agg_BS_norm <- function(mu, sigma, y, threshold) {
    return(mean((pnorm(q = t, mean = mu, sd = sigma) - (y <= t)) ^ 2))
}


# --------------------- AGGREGATED PATCHED SCORING RULES --------------------- #

agg_CRPS_mean <- function(mu, sigma, lambda, beta, y, p1, p2, my_grid) {
    res <- c()
    for (i in 1:(d1 - p1 + 1)) {
        for (j in 1:(d2 - p2 + 1)) {
            y_P <- mean(y[i:(i + p1 - 1), j:(j + p2 - 1)])
            mu_P <- mean(mu[i:(i + p1 - 1), j:(j + p2 - 1)])
            coord_patch <- patch_coord(my_grid, i, j, p1, p2)
            ind_patch <- patch_ind(my_grid, i, j, p1, p2)

            diag_sigma <- diag(c(sigma[i:(i + p1 - 1), j:(j + p2 - 1)]))
            corr <- exp(-(as.matrix(dist(coord_patch, diag = T, upper = T)) / lambda) ^ beta)

            sigma_P <- sqrt(mean(diag_sigma %*% corr %*% diag_sigma))
            res = c(res, crps_norm(y = y_P, mean = mu_P, sd = sigma_P))
        }
    }
    return(mean(res))
}

agg_BS_fte <- function(ens, y, p1, p2, my_grid, t) {
    res <- c()
    for (i in 1:(d1 - p1 + 1)) {
        for (j in 1:(d2 - p2 + 1)) {
            coord_patch <- patch_coord(my_grid, i, j, p1, p2)
            ind_patch <- patch_ind(my_grid, i, j, p1, p2)

            y_P <- apply(ind_patch, 1, function(x) { y[x[1], x[2]] })
            ens_P <- t(apply(ind_patch, 1, function(x) { ens[x[1], x[2],] }))

            res = c(res, (mean(ens_P >= t) - mean(y_P >= t)) ^ 2)
        }
    }
    return(mean(res))
}