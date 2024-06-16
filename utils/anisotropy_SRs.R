# ---------------------------------------------------------------------------- #
#                           ANISOTROPY SCORING RULES                           #
# ---------------------------------------------------------------------------- #

# ----------------------------- ANISOTROPIC SCORE ---------------------------- #

vario <- function(X, h1, h2) {
    d1 <- dim(X)[1]
    d2 <- dim(X)[2]

    X_ij <- X[max(1, 1 - h1):min(d1, d1 - h1), max(1, 1 - h2):min(d2, d2 - h2)]
    X_ijh <- X[max(1 + h1, 1):min(d1, d1 + h1), max(1, 1 + h2):min(d2, d2 + h2)]
    return(
    mean(
      abs(
        X_ij - X_ijh
        )
      )
    )
}

T_iso_old <- function(X, h) {
    v1 <- vario(X, h, 0)
    v2 <- vario(X, 0, h)
    v3 <- vario(X, h, h)
    v4 <- vario(X, - h, h)
    return(
    -(((v1 - v2) / (v1 + v2)) ^ 2 + ((v3 - v4) / (v3 + v4)) ^ 2)
    )
}

T_iso_a1 <- function(X, h) {
    d1 <- dim(X)[1]
    d2 <- dim(X)[2]
    v1 <- vario(X, h, 0)
    v2 <- vario(X, 0, h)
    return(
    -(v1 - v2) ^ 2 / (2 * v1 ^ 2 / ((d1 - h) * d2) + 2 * v2 ^ 2 / (d1 * (d2 - h)))
    )
}

T_iso_a2 <- function(X, h) {
    d1 <- dim(X)[1]
    d2 <- dim(X)[2]
    v1 <- vario(X, h, h)
    v2 <- vario(X, - h, h)
    return(
    -(v1 - v2) ^ 2 / (2 * v1 ^ 2 / ((d1 - h) * (d2 - h)) + 2 * v2 ^ 2 / ((d1 - h) * (d2 - h)))
    )
}

SE_Tiso_old <- function(ens, y, h) {
    mean_ens_Tiso <- mean(apply(ens, 3, function(x) { T_iso_old(X = x, h = h) }))
    return((mean_ens_Tiso - T_iso_old(X = y, h = h)) ^ 2)
}
SE_Tiso_a1 <- function(ens, y, h) {
    mean_ens_Tiso <- mean(apply(ens, 3, function(x) { T_iso_a1(X = x, h = h) }))
    return((mean_ens_Tiso - T_iso_a1(X = y, h = h)) ^ 2)
}
SE_Tiso_a2 <- function(ens, y, h) {
    mean_ens_Tiso <- mean(apply(ens, 3, function(x) { T_iso_a2(X = x, h = h) }))
    return((mean_ens_Tiso - T_iso_a2(X = y, h = h)) ^ 2)
}


# ------------------------------ VARIOGRAM SCORE ----------------------------- #

variogram_score <- function(mu, sigma, lambda, beta, y, p = .5, weights = matrix(rep(1, (d1 * d2) ^ 2), nrow = d1 * d2) / (d1 * d2) ^ 2, my_grid, A = NULL) {
    coord <- coords(my_grid)
    ind <- indices(my_grid)
    v_mu <- apply(ind, 1, function(i) { mu[i[1], i[2]] })
    v_y <- apply(ind, 1, function(i) { y[i[1], i[2]] })
    if (is.null(A)) {
        D <- as.matrix(dist(coord, diag = T, upper = T))
    } else {
        ratio <- A@par.model$ratio
        a <- A@par.model$angle
        D <- mev::distg(loc = coord, scale = 1 / ratio, rho = -a)
    }
    sigma2 <- 2 * sigma ^ 2 * (1 - exp(-(D / lambda) ^ beta))
    temp <- -(outer(v_mu, v_mu, '-') ^ 2) / (2 * sigma2)
    diag(temp) <- rep(0, length(v_mu))

    vs <- (sigma2 ^ (p / 2) * 2 ^ (p / 2) * gamma((p + 1) / 2) / sqrt(pi) * hyperg_1F1(a = -p / 2, b = .5, x = temp) - as.matrix(dist(v_y, diag = T, upper = T)) ^ p) ^ 2
    return(sum(weights * vs))
}