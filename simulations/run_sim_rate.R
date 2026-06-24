library(tidyverse)
library(MASS)
library(tmvtnorm)
library(truncnorm)
library(parallel)
library(viridis)
library(refund)
library(pracma)
library(caret)
library(gbm)
library(LinCDE)
library(randomForest)
library(ranger)
source('input_funcs_new.R')
source('mid_func.R')


i <- 1
val <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
n_options <- c(1000,2000,4000,8000,16000,32000)
mult <- trunc((val - 1)/40)
new_val <- val - mult*40
seeds <- seq((new_val-1)*25 +10000001, new_val*25 + 10000000)
n_df <- n_options[mult+1]

input_grid <- expand.grid('b1' = c(2,5,8), 'b3' = c(2,5,8), 'corr' = c(0, 0.2),'q_sd' = c(10, 12))


if(input_grid$q_sd[i] == 10){
  delta_vec <- c(0.03, 0.045, 0.06, 0.09, 0.12, 0.16)
} else if(input_grid$q_sd[i] == 12){
  delta_vec <- c(0.02, 0.03, 0.04, 0.06, 0.08, 0.1)
}


func_coef_out <- gen_funcs_and_coefs(52545, 150, 200, 10, input_grid$b1[i], 2, input_grid$b3[i])



true_mu <- read_rds(paste0('true_res_Q/', i, '.rds'))[1]
true_mu_q <- read_rds(paste0('true_res_Q/', i, '.rds'))[2:7]
lb <- read_rds(paste0('true_res_Q/', i, '.rds'))[8]
ub <- read_rds(paste0('true_res_Q/', i, '.rds'))[9]


res <- mclapply(seeds, gen_sim,
              n=n_df, p=10, b1=input_grid$b1[i], b2=2, b3=input_grid$b3[i], t1=150, t2=200,
              corr=input_grid$corr[i], q_sd = input_grid$q_sd[i], func_coef_out = func_coef_out,
              lb = lb, ub = ub)

m_vals <- seeds - 10000000
est_res <- mclapply(seq(1, length(res)), function(m){
  print(paste0('df # ', m_vals[m]))
  final_df <- lapply(c(0.99), function(k){
    q_sd <- unique(res[[m]][[1]]$q_sd)
    print(k)
    mod_df <- get_coef_fpca_df(res[[m]], k)
    lincde_fpca <- tryCatch(
      {
        set.seed(85675 + 1000*m_vals[m])
        vec1 <- sort(sample(seq(1, n_df), n_df/2, replace = FALSE))
        vec2 <- seq(1, n_df)[!(seq(1, n_df) %in% vec1)]
        mod_df1 <- mod_df[vec1,]
        mod_df2 <- mod_df[vec2,]
        
        design_pts <- seq(lb, ub, length.out = 200)
        hyper_grid <- cbind(expand.grid('mtry' = seq(max(3, round((ncol(mod_df)-1)/3) -7) , round((ncol(mod_df)-1)/3) + 5, by = 1)),
                            splitrule = 'variance', min.node.size = 5)
        
        print('lincde fpca')
        dr_lincde1_fpca  <- get_w3_t2_fast(n.trees = 1000, depth = 2, centering = F, centeringMethod = "linearRegression",
                                   design_pts, delta_vec, mod_df1 , mod_df2 , 'coef1_2', hyper_grid, m_vals[m]*2000, cores = 1)
        dr_lincde2_fpca  <- get_w3_t2_fast(n.trees = 1000, depth = 2,centering = F, centeringMethod = "linearRegression",
                                   design_pts, delta_vec, mod_df2 , mod_df1 , 'coef1_2', hyper_grid, m_vals[m]*4020, cores = 1)
        
        cbind(type = 'est_fpca', get_est(dr_lincde1_fpca, dr_lincde2_fpca, mod_df1$y, mod_df2$y, delta_vec))
      },
      error = function(e) {
        tryCatch({
          print('second')
          set.seed(85675 + 1000*m + 5)
          vec1 <- sort(sample(seq(1, n_df), n_df/2, replace = FALSE))
          vec2 <- seq(1, n_df)[!(seq(1, n_df) %in% vec1)]
          mod_df1 <- mod_df[vec1,]
          mod_df2 <- mod_df[vec2,]
          
          design_pts <- seq(lb, ub, length.out = 200)
          hyper_grid <- expand.grid('mtry' = seq(max(3, round((ncol(mod_df)-1)/3) -7), round((ncol(mod_df)-1)/3) + 7, by = 1))
          
          dr_lincde1_fpca  <- get_w3_t2_fast(n.trees = 1000, depth = 2, centering = F, centeringMethod = "linearRegression",
                                     design_pts, delta_vec, mod_df1 , mod_df2 , 'coef1_2', hyper_grid, m_vals[m]*2000)

          dr_lincde2_fpca  <- get_w3_t2_fast(n.trees = 1000, depth =  2,centering = F, centeringMethod = "linearRegression",
                                     design_pts, delta_vec, mod_df2 , mod_df1 , 'coef1_2', hyper_grid, m_vals[m]*4020)
          cbind(type = 'est_fpca', get_est(dr_lincde1_fpca, dr_lincde2_fpca, mod_df1$y, mod_df2$y, delta_vec))
        },
        error = function(e){
          fake_df <- as.data.frame(matrix(rep(NA, 45), nrow = 5))
          colnames(fake_df) <- c('type', 'delta', 'out_mod', 'w_norm', 'dr_norm', 'dr_norm2', 'var_est', 'var_est2', 'var_contrast_est')
          fake_df$type = 'est_fpca'
          return(fake_df)
        })
      }
    )
    

    
    final <- as.data.frame(rbind(lincde_fpca)) %>%
      mutate(seed = m,
             est_y = mean(c(mod_df1$y, mod_df2$y)),
             true_mu = true_mu,
             true_mu_q = rep(true_mu_q, 1),
             pve = k,
             lb = lb,
             ub = ub,
             n_col = ncol(mod_df))
    return(final)
  })

    return( Reduce(rbind, final_df))
}

, mc.cores = 25)

final_res <- do.call(rbind, est_res)
write_rds(final_res, paste0('res_rate_Q/', val, '.rds'))



