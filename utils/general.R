# ---------------------------------------------------------------------------- #
#                               GENERAL FUNCTIONS                              #
# ---------------------------------------------------------------------------- #

mean_chunk <- function(x, chunk_size, stride) {
    r <- c()
    for (i in seq(1, length(x) - chunk_size + 1, stride)) {
        r <- c(r, mean(x[i:(i + chunk_size - 1)]))
    }
    return(r)
}

coords <- function(my_grid){
  n1 <- length(my_grid$x)
  n2 <- length(my_grid$y)
  return(cbind(rep(my_grid$x,n2),c(sapply(my_grid$y,function(a){rep(a,n1)}))))
}
indices <- function(my_grid){
  n1 <- length(my_grid$ind.x)
  n2 <- length(my_grid$ind.y)
  return(cbind(rep(my_grid$ind.x,n2),c(sapply(my_grid$ind.y,function(a){rep(a,n1)}))))
}

coords_var <- function(my_grid){
  n1 <- length(my_grid$x)-1
  n2 <- length(my_grid$y)-1
  return(cbind(rep(my_grid$x[-length(my_grid$x)],n2),c(sapply(my_grid$y[-length(my_grid$y)],function(a){rep(a,n1)}))))
}
indices_var <- function(my_grid){
  n1 <- length(my_grid$ind.x)-1
  n2 <- length(my_grid$ind.y)-1
  return(cbind(rep(my_grid$ind.x[-length(my_grid$ind.x)],n2),c(sapply(my_grid$ind.y[-length(my_grid$ind.y)],function(a){rep(a,n1)}))))
}

patch_coord <- function(my_grid,i,j,p1,p2){
  return(cbind(rep(my_grid$x[i:(i+p1-1)],p2),c(sapply(my_grid$y[j:(j+p2-1)],function(a){rep(a,p1)}))))
}
patch_ind <- function(my_grid,i,j,p1,p2){
  return(cbind(rep(my_grid$ind.x[i:(i+p1-1)],p2),c(sapply(my_grid$ind.y[j:(j+p2-1)],function(a){rep(a,p1)}))))
}