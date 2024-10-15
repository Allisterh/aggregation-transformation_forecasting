# Required packages
packages <- c("scoringRules","BSDA","doParallel","foreach","tidyr","dplyr","R.utils")
lapply(packages, require, ch = TRUE)
source("../utils/case_study.R")

## Import data from MultivCalibration package
load("wind_dat.rda") # wind_dat.rda dowloaded from MultivCalibration package : https://github.com/sallen12/MultivCalibration/tree/main/data/wind_dat.rda
list2env(wind_dat, globalenv())
rm(wind_dat)

### Data dimensions
n <- 1045
p <- 33
q <- 32
M <- 11


### Prerequisites

grille <- list(x=unique(coord[,1])[1:p],y=unique(coord[,2])[1:q],
               ind.x=1:p,ind.y=1:q)

coords <- function(grille){
  n1 <- length(grille$x)
  n2 <- length(grille$y)
  return(cbind(rep(grille$x,n2),c(sapply(grille$y,function(a){rep(a,n1)}))))
}

patch_ind <- function(grille,i,j,p1,p2){
  return(cbind(rep(grille$ind.x[i:(i+p1-1)],p2),c(sapply(grille$ind.y[j:(j+p2-1)],function(a){rep(a,p1)}))))
}
patch_coord <- function(grille,i,j,p1,p2){
  return(cbind(rep(grille$x[i:(i+p1-1)],p2),c(sapply(grille$y[j:(j+p2-1)],function(a){rep(a,p1)}))))
}

indices_var <- function(grille){
  n1 <- length(grille$ind.x)-1
  n2 <- length(grille$ind.y)-1
  return(cbind(rep(grille$ind.x[-length(grille$ind.x)],n2),c(sapply(grille$ind.y[-length(grille$ind.y)],function(a){rep(a,n1)}))))
}



weights_d <- 1/as.matrix(dist(coords(grille),diag=T,upper=T))
diag(weights_d) <- rep(0,nrow(weights_d))
weights_d <- weights_d/sum(weights_d)



## Computing scoring rules

ens_forecasts <- list(fc_ifs=fc_ifs[1:n,1:p,1:q,1:M], fc_ecc=fc_ecc[1:n,1:p,1:q,1:M], fc_ss=fc_ss[1:n,1:p,1:q,1:M])
obs_temp <- obs[1:n,1:p,1:q]
weights_d <- weights_d[1:(p*q),1:(p*q)]

t_list <- c(2.5,5,7.5,10)
alpha_list <- c(.7,.8,.9)
p_vs_list <- c(0.5,1,2)
p_pvs_list <- c(0.5,1,2)
patch_list_es <- c(3)
h_list <- c(1,2,3,4,5)
patch_list_dp <- c(3)
t_list_fte <- c(2.5,5,7.5,10)

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

df_exp <- foreach(m = c("fc_ifs","fc_ecc","fc_ss"), .combine = rbind, .packages=c('scoringRules')) %dopar% {
  df <- data.frame(model=character(),score_name=character(),score_value=numeric())
  ens = ens_forecasts[[m]]
  # ---------------------------------------------------------------------------- #
  #                                  UNIVARIATE                                  #
  # ---------------------------------------------------------------------------- #
  
  ### Agg. CRPS
  df_temp <- data.frame(model = rep(m,n),
                        score_name = rep("agg_crps",n),
                        score_value = agg_CRPS(ens=ens,obs=obs_temp,method='edf')
  )
  df <- rbind(df,df_temp)
  
  
  ### Agg. SE
  df_temp <- data.frame(model = rep(m,n),
                        score_name = rep("agg_se",n),
                        score_value = agg_SE(ens=ens,obs=obs_temp)
  )
  df <- rbind(df,df_temp)
  
  
  ### Agg. DSS
  df_temp <- data.frame(model = rep(m,n),
                        score_name = rep("agg_dss",n),
                        score_value = agg_DSS(ens=ens,obs=obs_temp)
  )
  df <- rbind(df,df_temp)
  
  ### Agg. BS
  for (t in t_list){
    df_temp <- data.frame(model = rep(m,n),
                          score_name = rep(paste("agg_bs_t",t,sep=''),n),
                          score_value = agg_BS(ens=ens,obs=obs_temp,threshold=t)
    )
    df <- rbind(df,df_temp)
  }
  
  ### Agg. QS
  for (alpha in alpha_list){
    df_temp <- data.frame(model = rep(m,n),
                          score_name = rep(paste("agg_qs_a",alpha,sep=''),n),
                          score_value = agg_QS(ens=ens,obs=obs_temp,alpha=alpha)
      )
    df <- rbind(df,df_temp)
  }
  
  # ---------------------------------------------------------------------------- #
  #                                 MULTIVARIATE                                 #
  # ---------------------------------------------------------------------------- #
  
  ### Variogram score
  for (p_vs in p_vs_list){
    df_temp <- data.frame(model = rep(m,n),
                          score_name = rep(paste0("vs_",p_vs),n),
                          score_value=variogram_score(ens=ens, obs=obs_temp, p_vs=p_vs, w_vs=weights_d)
      )
      df <- rbind(df,df_temp)
  }
  
  ### Power-Variation score
  for (p_pvs in p_pvs_list){
    df_temp <- data.frame(model = rep(m,n),
                          score_name = rep(paste0("pvs_",p_pvs),n),
                          score_value = pvar_score_sample(ens=ens ,obs=obs_temp, p=p_pvs,grille=grille)
      )
      df <- rbind(df,df_temp)
  }
  
  ### Patched ES
  for (patch_size in patch_list_es){
    df_temp <- data.frame(model = rep(m,n),
                          score_name = rep(paste0("es_p",patch_size),n),
                          score_value=patched_es_sample(ens=ens, obs=obs_temp,p1=patch_size,p2=patch_size,grille=grille)
    )
    df <- rbind(df,df_temp)
  }
  
  ### Energy Score
  df_temp <- data.frame(model = rep(m,n),
                          score_name = rep("es",n),
                          score_value=energy_score(ens=ens, obs=obs_temp)
    )
  df <- rbind(df,df_temp)
  
  ### Anisotropic score (a1)
  for (h in h_list){
    df_temp <- data.frame(model = rep(m,n),
                          score_name = rep(paste0("as_a1_h",h),n),
                          score_value = SE_Tiso_a1(ens=ens, obs_temp,h=h)
    )
    df <- rbind(df,df_temp)
  }
  
  ### Anisotropic score (a2)
  for (h in h_list){
    df_temp <- data.frame(model = rep(m,n),
                          score_name = rep(paste0("as_a2_h",h),n),
                          score_value = SE_Tiso_a2(ens=ens, obs_temp,h=h)
    )
    df <- rbind(df,df_temp)
  }
  df

  ### Agg. CRPS of spatial mean
  for (patch_size in patch_list_dp){
    df_temp <- data.frame(model = rep(m,n),
                        score_name = rep(paste0("agg_crps_sm_p",patch_size),n),
                        score_value = agg_CRPS_mean(ens=ens,obs=obs_temp,p1=patch_size,p2=patch_size)
    )
    df <- rbind(df,df_temp)
  }
  
  ### Agg. Brier Score of FTE
  for (patch_size in patch_list_dp){
    for (t in t_list_fte){
      df_temp <- data.frame(model = rep(m,n),
                            score_name = rep(paste0("agg_se_fte_t",t,"_p",patch_size),n),
                            score_value = agg_SE_FTE(ens=ens,obs=obs_temp,p1=patch_size,p2=patch_size,grille=grille,t=t)
      )
      df <- rbind(df,df_temp)
    }
  }
    
  df
}
stopCluster(cl)


## Diebold-Mariano tests

### Univariate scoring rules

list_models <- unique(df_exp$model)
list_univ <- c("agg_se","agg_dss","agg_bs_t2.5","agg_bs_t5","agg_bs_t7.5","agg_bs_t10","agg_qs_a0.7","agg_qs_a0.8","agg_qs_a0.9","agg_crps")

df_avg <- df_exp %>% group_by(model,score_name) %>% summarize(avg_score = mean(score_value))
df_dm_univ <- data.frame(score_name=character(),model1=character(),model2=character(),dm_test=character(),p_value=numeric(),z_test=numeric())

for (sr in list_univ){
    for (m_ref in list_models){
      s2 <- subset(df_exp,model==m_ref & score_name==sr)$score_value
      for (m in list_models){
        if (m!=m_ref && ((sr %in% list_univ) && !setequal(c(m_ref,m),c('fc_ecc','fc_ss'))) | !(sr %in% list_univ) ){
          s1 <- subset(df_exp,model==m & score_name==sr)$score_value
          x <- (s1-s2)/sd(s1-s2)
          
          z_test = z.test(x,sigma.x=1)
          if (z_test$estimate>0){
            df_temp <- data.frame(score_name = sr,
                                  model1 = m_ref,
                                  model2 = m,
                                  dm_test = z_test$p.value<.05,
                                  p_value = z_test$p.value,
                                  z_test = z_test$estimate
            )
            df_dm_univ <- rbind(df_dm_univ,df_temp)
          }
          
        }
      }
    }
  # }
}
rownames(df_dm_univ) <- NULL
print(df_dm_univ)

### Multivariate scoring rules

list_models <- unique(df_exp$model)
list_multiv <- c("es","es_p3","vs_0.5","vs_1","vs_2","pvs_0.5","pvs_1","pvs_2","as_a1_h1","as_a1_h2","as_a1_h3","as_a1_h4","as_a1_h5","agg_crps_sm_p3","agg_se_fte_t2.5_p3","agg_se_fte_t5_p3","agg_se_fte_t7.5_p3","agg_se_fte_t10_p3")

df_avg <- df_exp %>% group_by(model,score_name) %>% summarize(avg_score = mean(score_value))
df_dm_multiv <- data.frame(score_name=character(),model1=character(),model2=character(),dm_test=character(),p_value=numeric(),z_test=numeric())

for (sr in list_multiv){
  for (m_ref in list_models){
    s2 <- subset(df_exp,model==m_ref & score_name==sr)$score_value
    for (m in list_models){
      if ((m!=m_ref && !grepl("agg_se_fte",sr)) || (grepl("agg_se_fte",sr) && !setequal(c(m_ref,m),c('fc_ecc','fc_ss')) && m!=m_ref) ){
        s1 <- subset(df_exp,model==m & score_name==sr)$score_value
        x <- (s1-s2)/sd(s1-s2)
        
        z_test = z.test(x,sigma.x=1)
        if (z_test$estimate>0){
          df_temp <- data.frame(score_name = sr,
                                model1 = m_ref,
                                model2 = m,
                                dm_test = z_test$p.value<.05,
                                p_value = z_test$p.value,
                                z_test = z_test$estimate
          )
          df_dm_multiv <- rbind(df_dm_multiv,df_temp)
        }
        
      }
    }
  }
}
rownames(df_dm_multiv) <- NULL
print(df_dm_multiv)


## Expected scoring rules

dff <- df_exp %>%
  group_by(score_name,model) %>%
  dplyr::summarize(avg_score = mean(score_value)) %>% 
  pivot_wider(
    names_from = score_name, 
    values_from = avg_score
  ) %>%
    slice(c(2,1,3))


dff %>%
  select(model,agg_se,agg_dss,agg_bs_t2.5,agg_bs_t5,agg_bs_t7.5,agg_bs_t10,agg_qs_a0.7,agg_qs_a0.8,agg_qs_a0.9,agg_crps)

dff %>%
  select(model,es,es_p3,vs_0.5,vs_1,vs_2,pvs_0.5,pvs_1,pvs_2,as_a1_h1,as_a1_h2,as_a1_h3,agg_crps_sm_p3,agg_se_fte_t2.5_p3,agg_se_fte_t5_p3,agg_se_fte_t7.5_p3,agg_se_fte_t10_p3)