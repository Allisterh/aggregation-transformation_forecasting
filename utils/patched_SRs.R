# ---------------------------------------------------------------------------- #
#                             PATCHED SCORING RULES                            #
# ---------------------------------------------------------------------------- #

# ---------------------------- REQUIRED LIBRARIES ---------------------------- #

require(scoringRules)

# --------------------------------- FUNCTIONS -------------------------------- #

agg_CRPS_norm <- function(mu, sigma, y) {
    return(mean(scoringRules::crps_norm(y = y, mean = mu, sd = sigma)))
}

variogram_score <- function(mu, sigma, lambda, beta, y, p = .5, weights = matrix(rep(1, (d1 * d2) ^ 2), nrow = d1 * d2) / (d1 * d2) ^ 2, my_grid) {
    coord <- coords(my_grid)
    ind <- indices(my_grid)
    v_mu <- apply(ind, 1, function(i) { mu[i[1], i[2]] })
    v_y <- apply(ind, 1, function(i) { y[i[1], i[2]] })

    sigma2 <- 2 * sigma ^ 2 * (1 - exp(-(as.matrix(dist(coord, diag = T, upper = T)) / lambda) ^ beta))
    temp <- -(outer(v_mu, v_mu, '-') ^ 2) / (2 * sigma2)
    diag(temp) <- rep(0, length(v_mu))

    vs <- (sigma2 ^ (p / 2) * 2 ^ (p / 2) * gamma((p + 1) / 2) / sqrt(pi) * gsl::hyperg_1F1(a = -p / 2, b = .5, x = temp) - as.matrix(dist(v_y, diag = T, upper = T)) ^ p) ^ 2
    return(sum(weights * vs))
}

p_variation_score <- function(sigma, lambda, beta, y, p = .5, weights = NA, my_grid) {
    coord <- coords_var(my_grid)
    ind <- indices_var(my_grid)

    pv_y <- apply(ind, 1, function(i) { abs(y[i[1] + 1, i[2] + 1] - y[i[1] + 1, i[2]] - y[i[1], i[2] + 1] + y[i[1], i[2]]) ^ p })

    sigma2 <- 4 * sigma ^ 2 * (1 + exp(-(sqrt(2) / lambda) ^ beta) - 2 * exp(-(1 / lambda) ^ beta))

    pvs <- mean((sigma2 ^ (p / 2) * 2 ^ (p / 2) * gamma((p + 1) / 2) / sqrt(pi) - pv_y) ^ 2)
    return(pvs)
}


patched_energy_score <- function(ens, y, p1, p2, my_grid) {
    res <- c()
    for (i in 1:(d1 - p1 + 1)) {
        for (j in 1:(d2 - p2 + 1)) {
            coord_patch <- patch_coord(my_grid, i, j, p1, p2)
            ind_patch <- patch_ind(my_grid, i, j, p1, p2)

            y_P <- apply(ind_patch, 1, function(x) { y[x[1], x[2]] })
            ens_P <- t(apply(ind_patch, 1, function(x) { ens[x[1], x[2],] }))

            res = c(res, scoringRules::es_sample(y = y_P, dat = ens_P))
        }
    }
    return(mean(res))
}

energy_score <- function(ens, y, my_grid) {
    coord <- coords(my_grid)
    ind <- indices(my_grid)

    y <- apply(ind, 1, function(x) { y[x[1], x[2]] })
    ens <- t(apply(ind, 1, function(x) { ens[x[1], x[2],] }))

    return(es_sample(y = y, dat = ens))
}