library(tidyverse)
library(MASS)
library(tmvtnorm)
library(truncnorm)
library(parallel)
library(viridis)
library(refund)
library(pracma)
library(caret)
library(LinCDE)
library(randomForest)
source('input_funcs_qstar.R')
source('mid_func.R')


i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

input_grid <- expand.grid('b1' = c(2,5,8), 'b3' = c(2,5,8),'q_sd' = c(10, 12), 'b3_coef' = c(0,3), 'b3_sd' = c(1, 3))


if(input_grid$q_sd[i] == 10){
  delta_vec <- c(0.03, 0.045, 0.06, 0.09, 0.12, 0.16)
} else if(input_grid$q_sd[i] == 12){
  delta_vec <- c(0.02, 0.03, 0.04, 0.06, 0.08, 0.1)
}


func_coef_out <- gen_funcs_and_coefs(52545, 150, 200, 10, input_grid$b1[i], 2, input_grid$b3[i], input_grid$b3_coef[i])


## get true values ###
true_out <- gen_sim_no_func(seed = 5, n=1000000, p=10, b1=input_grid$b1[i], b2=2, b3=input_grid$b3[i],
                   t1=150, t2=200, corr=0, q_sd = input_grid$q_sd[i], b3_sd = input_grid$b3_sd[i],
                   func_coef_out = func_coef_out)
true_mu <- mean(true_out$mu)
lb <- unique(true_out$lb)
ub <- unique(true_out$ub)

q_true<-  mclapply(delta_vec, function(x){
  print(x)
  q_out <- gen_q_no_func(true_out, NA, func_coef_out = func_coef_out, delta = x, cores = 4)
  return(mean(q_out$mu_q))
}, mc.cores = 6)
true_mu_q <- unlist(q_true)
write_rds(c(true_mu, true_mu_q, lb, ub), paste0('true_res_Qstar/', i, '.rds'))

#uncomment after the true results have been saved
# true_mu <- read_rds(paste0('true_res_Qstar/', i, '.rds'))[1]
# true_mu_q <- read_rds(paste0('true_res_Qstar/', i, '.rds'))[2:7]
# lb <- read_rds(paste0('true_res_Qstar/', i, '.rds'))[8]
# ub <- read_rds(paste0('true_res_Qstar/', i, '.rds'))[9]


res <- mclapply(seq(10000001, 10000001+999), gen_sim,
              n=1000, p=10, b1=input_grid$b1[i], b2=2, b3=input_grid$b3[i], t1=150, t2=200,
              corr=0, q_sd = input_grid$q_sd[i], b3_sd = input_grid$b3_sd[i], func_coef_out = func_coef_out,
              lb = lb, ub = ub,
              mc.cores = 25)

est_res <- mclapply(seq(1, length(res)), function(m){
  print(paste0('df # ', m))
  final_df <- lapply(c(0.99), function(k){
    q_sd <- unique(res[[m]][[1]]$q_sd)
    print(k)
    mod_df <- get_coef_fpca_df(res[[m]], k)
    fpca_num <- c(sum(str_detect(colnames(mod_df), 'coef') == TRUE & str_detect(colnames(mod_df), '_') == FALSE),
                  sum(str_detect(colnames(mod_df), 'coef') == TRUE & str_detect(colnames(mod_df), '_2') == TRUE),
                  sum(str_detect(colnames(mod_df), 'coef') == TRUE & str_detect(colnames(mod_df), '_3') == TRUE))
    lincde_fpca <- tryCatch(
      {
        set.seed(85675 + 1000*m)
        vec1 <- sort(sample(seq(1, 1000), 500, replace = FALSE))
        vec2 <- seq(1, 1000)[!(seq(1, 1000) %in% vec1)]
        mod_df1 <- mod_df[vec1,]
        mod_df2 <- mod_df[vec2,]
        
        design_pts <- seq(lb, ub, length.out = 200)
        hyper_grid <- expand.grid('mtry' = seq(max(3, round((ncol(mod_df)-1)/3) -7) , round((ncol(mod_df)-1)/3) + 5, by = 1))
        
        print('lincde fpca')
        dr_lincde1_fpca  <- get_w3_star(n.trees = 1000, depth = 2, centering = F, centeringMethod = "linearRegression",
                                   design_pts, delta_vec, mod_df1 , mod_df2 , 'coef1_2', hyper_grid, m*2000)
        dr_lincde2_fpca  <- get_w3_star(n.trees = 1000, depth = 2,centering = F, centeringMethod = "linearRegression",
                                   design_pts, delta_vec, mod_df2 , mod_df1 , 'coef1_2', hyper_grid, m*4020)
        
        cbind(type = 'est_fpca', get_est(dr_lincde1_fpca, dr_lincde2_fpca, mod_df1$y, mod_df2$y, delta_vec))
      },
      error = function(e) {
        tryCatch({
          print('second')
          set.seed(85675 + 1000*m + 5)
          vec1 <- sort(sample(seq(1, 1000), 500, replace = FALSE))
          vec2 <- seq(1, 1000)[!(seq(1, 1000) %in% vec1)]
          mod_df1 <- mod_df[vec1,]
          mod_df2 <- mod_df[vec2,]
          
          design_pts <- seq(lb, ub, length.out = 200)
          hyper_grid <- expand.grid('mtry' = seq(max(3, round((ncol(mod_df)-1)/3) -7), round((ncol(mod_df)-1)/3) + 7, by = 1))
          
          dr_lincde1_fpca  <- get_w3_star(n.trees = 1000, depth = 2, centering = F, centeringMethod = "linearRegression",
                                     design_pts, delta_vec, mod_df1 , mod_df2 , 'coef1_2', hyper_grid, m*2000)

          dr_lincde2_fpca  <- get_w3_star(n.trees = 1000, depth =  2,centering = F, centeringMethod = "linearRegression",
                                     design_pts, delta_vec, mod_df2 , mod_df1 , 'coef1_2', hyper_grid, m*4020)
          cbind(type = 'est_fpca', get_est(dr_lincde1_fpca, dr_lincde2_fpca, mod_df1$y, mod_df2$y, delta_vec))
        },
        error = function(e){
          fake_df <- as.data.frame(matrix(rep(NA, 50), nrow = 5))
          colnames(fake_df) <- c('type', 'delta', 'out_mod', 'w_norm', 'dr_norm', 'dr_norm2', 'var_est', 'var_est2', 'var_contrast_est', 'var_contrast_est2')
          fake_df$type = 'est_fpca'
          return(fake_df)
        })
      }
    )
  

    
    final <- as.data.frame(rbind(lincde_fpca)) %>%
      mutate(seed = m,
             est_y = mean(c(mod_df1$y, mod_df2$y)),
             true_mu = true_mu,
             true_mu_q = rep(true_mu_q),
             pve = k,
             lb = lb,
             ub = ub,
             n_col = ncol(mod_df),
             n_fpca1 = fpca_num[1],
             n_fpca2 = fpca_num[2],
             n_fpca3 = fpca_num[3])
    return(final)
  })

    return( Reduce(rbind, final_df))
}

, mc.cores = 25)

final_res <- do.call(rbind, est_res)
write_rds(final_res, paste0('res_Qstar/', i, '.rds'))



