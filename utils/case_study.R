# ---------------------------------------------------------------------------- #
#                SCORING RULES FOR ENSEMBLE FORECASTS ON A GRID                #
# ---------------------------------------------------------------------------- #

# ---------------------------- REQUIRED LIBRARIES ---------------------------- #

require(scoringRules)

# ------------------------- UNIVARIATE SCORING RULES ------------------------- #

### Agg. Continuous Ranked Probability Score
agg_CRPS_idx <- function(idx,ens,obs,method='edf'){
  rs_ens <- ens[idx,1:p,1:q,1:M]
  rs_obs <- obs[idx,1:p,1:q]
  dim(rs_ens) <- c(p*q,M)
  dim(rs_obs) <- c(p*q)
  rs_obs <- as.vector(rs_obs)
  return(
    mean(crps_sample(y=rs_obs,dat=rs_ens,method=method))
    )
}

agg_CRPS <- function(ens,obs,method="edf"){
  return(
    unlist(lapply(1:n,function(i){agg_CRPS_idx(idx=i,ens=ens,obs=obs,method=method)}))
    )
}

### Agg. Squared Error
agg_SE <- function(ens,obs){
  return(
    apply(( obs - apply(ens,c(1,2,3),mean) )^2, 1, mean)
    )
}

### Agg. Dawid-Sebastiani Score
agg_DSS_idx <- function(idx,ens,obs){
  rs_ens <- ens[idx,1:p,1:q,1:M]
  rs_obs <- obs[idx,1:p,1:q]
  dim(rs_ens) <- c(p*q,M)
  dim(rs_obs) <- c(p*q)
  rs_obs <- as.vector(rs_obs)

  return(
    mean(dss_sample(y=rs_obs,dat=rs_ens))
    )
}

agg_DSS <- function(ens,obs){
  return(
    unlist(lapply(1:n,function(i){agg_DSS_idx(idx=i,ens=ens,obs=obs)}))
    )
}

### Agg. Brier Score
agg_BS <- function(ens,obs,threshold){
  return(
    apply( ( apply(ens,c(1,2,3),function(x){sum(x>threshold)/M}) - (obs>threshold) )^2, 1, mean )
  )
}

### Agg. Quantile Score
agg_QS <- function(ens,obs,alpha){
  q <- apply(ens,c(1,2,3),function(x){quantile(x,probs=alpha)})
  return(
    apply( ((obs<=q)-alpha) * (q-obs), 1, mean )
    )
}

# ------------------------ MULTIVARIATE SCORING RULES ------------------------ #

### Variogram Score (VS)
variogram_score <- function(ens,obs,p_vs,w_vs=matrix(rep(1,(p*q)^2),nrow=p*q)/(p*q)){#w_vs,
  print(dim(w_vs))
  rs_ens <- ens
  rs_obs <- obs
  dim(rs_ens) <- c(n,p*q,M)
  dim(rs_obs) <- c(n,p*q)
  
  return(
    unlist(lapply(1:n, function(i){vs_sample(y=rs_obs[i,1:(p*q)], dat=rs_ens[i,1:(p*q),1:M], p=p_vs, w_vs=w_vs)}))#, w_vs=w_vs
  )
}

### Power-Variation Score (PVS)
pvar_score_sample <- function(ens,obs,p_vs=.5,grille){
  ind <- indices_var(grille)

  pv_obs <- apply(ind,1,function(i){abs(obs[1:n,i[1]+1,i[2]+1] - obs[1:n,i[1]+1,i[2]] - obs[1:n,i[1],i[2]+1] + obs[1:n,i[1],i[2]])^p_vs})
  pv_ens <- apply(ind,
                  1,
                  function(i){
                    apply(
                      abs(ens[1:n,i[1]+1,i[2]+1,1:M] - ens[1:n,i[1]+1,i[2],1:M] - ens[1:n,i[1],i[2]+1,1:M] + ens[1:n,i[1],i[2],1:M])^p_vs,
                      1,
                      mean
                      )
                    }
                  )
  return(
    apply((pv_ens - pv_obs)^2, 1, mean)
  )
}

### Patched ES
patched_es_idx <- function(ens,obs,p1,p2,grille){
  res <- c()
  for (i in 1:(p-p1+1)){
    for (j in 1:(q-p2+1)){
      ind_patch <- patch_ind(grille,i,j,p1,p2)
      
      y_P <- apply(ind_patch,1,function(x){obs[x[1],x[2]]}) #obs[i:(i+p1-1),j:(j+p2-1)]
      ens_P <- t(apply(ind_patch,1,function(x){ens[x[1],x[2],]}))#ens[i:(i+p1-1),j:(j+p2-1),]
  
      res = c(res,scoringRules::es_sample(y=y_P,dat=ens_P))
    }
  }
  return(mean(res))
}

patched_es_sample <- function(ens,obs,p1,p2,grille){
  return(
    unlist(lapply(1:n, function(i){patched_es_idx(ens=ens[i,1:p,1:q,1:M],obs=obs[i,1:p,1:q],p1=p1,p2=p2,grille=grille)}))
  )
}

### Energy Score
energy_score <- function(ens,obs){
  rs_ens <- ens
  rs_obs <- obs
  dim(rs_ens) <- c(n,p*q,M)
  dim(rs_obs) <- c(n,p*q)
  
  return(
    unlist(lapply(1:n, function(i){es_sample(y=rs_obs[i,1:(p*q)], dat=rs_ens[i,1:(p*q),1:M])}))
  )
}

### Anisotropy score (AS)
vario <- function(X,h1,h2){
  X_ij <- X[1:n,max(1,1-h1):min(p,p-h1),max(1,1-h2):min(q,q-h2)]
  X_ijh <- X[1:n,max(1+h1,1):min(p,p+h1),max(1,1+h2):min(q,q+h2)]
  
  return(
    apply(
      abs(X_ij - X_ijh),
      1,
      mean
      )
    )
}

vario_ens <- function(X,h1,h2){
  X_ij <- X[1:n,max(1,1-h1):min(p,p-h1),max(1,1-h2):min(q,q-h2),1:M]
  X_ijh <- X[1:n,max(1+h1,1):min(p,p+h1),max(1,1+h2):min(q,q+h2),1:M]
  
  return(
    apply(
      abs(X_ij - X_ijh),
      c(1,4),
      mean
      )
    )
}

T_iso_a1 <- function(X,h){
  if (length(dim(X))>3){
    v1 <- vario_ens(X,h,0)
    v2 <- vario_ens(X,0,h)
    
  } else {
    v1 <- vario(X,h,0)
    v2 <- vario(X,0,h)
  }
  return(
    - (v1-v2)^2 / (2*v1^2/((p-h)*q) + 2*v2^2/(p*(q-h)))
    )
}

T_iso_a2 <- function(X,h){
  if (length(dim(X))>3){
    v1 <- vario_ens(X,h,h)
    v2 <- vario_ens(X,-h,h)
  } else {
    v1 <- vario(X,h,h)
    v2 <- vario(X,-h,h)
  }
  return(
    - (v1-v2)^2 / (2*v1^2/((p-h)*(q-h)) + 2*v2^2/((p-h)*(q-h)))
    )
}

SE_Tiso_a1 <- function(ens,obs,h){
  mean_ens_Tiso <- apply(T_iso_a1(X=ens,h=h),2,mean)
  return(
    (mean_ens_Tiso - T_iso_a1(X=obs,h=h))^2
    )
}

SE_Tiso_a2 <- function(ens,obs,h){
  mean_ens_Tiso <- apply(T_iso_a2(X=ens,h=h),2,mean)
  return(
    (mean_ens_Tiso - T_iso_a2(X=obs,h=h))^2
    )
}

### Aggregated SE of FTE 
agg_SE_FTE <- function(ens,obs,p1,p2,grille,t){
  res <- c()
  for (i in 1:(p-p1)){
    for (j in 1:(q-p2)){
      coord_patch <- patch_coord(grille,i,j,p1,p2)
      ind_patch <- patch_ind(grille,i,j,p1,p2)
      
      obs_P <- apply(ind_patch,1,function(x){obs[1:n,x[1],x[2]]})
      ens_P <- apply(ind_patch,1,function(x){ens[1:n,x[1],x[2],1:M]})
      
      dim(ens_P) <- c(n,M,p1*p2)
  
      res = cbind(res,(apply(ens_P>=t,1,mean) - apply(obs_P>=t,1,mean))^2)
    }
  }
  # print(dim(res))
  return(
    apply(res,1,mean)
    )
}

### Agg. CRPS of spatial mean
agg_crps_mean_idx <- function(ens_idx,obs_idx,p1,p2,method=method){
  res <- c()
  for (i in 1:(p-p1+1)){
    for (j in 1:(q-p2+1)){
      obs_P <- mean(obs_idx[i:(i+p1-1),j:(j+p2-1)])
      ens_P <- apply(ens_idx[i:(i+p1-1),j:(j+p2-1),1:M],3,mean)

      res = c(res,scoringRules::crps_sample(y=obs_P,dat=ens_P,method=method))
    }
  }
  return(mean(res))
}

agg_CRPS_mean <- function(ens,obs,p1,p2,method='edf'){
  return(
    unlist(lapply(1:n,function(i){agg_crps_mean_idx(ens=ens[i,1:p,1:q,1:M],obs=obs[i,1:p,1:q],p1=p1,p2=p2,method=method)}))
    )
}