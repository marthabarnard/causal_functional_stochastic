library(tidyverse)
library(MASS)
library(viridis)

#### Aggregate n = 1000 q() results ####

# join saved files
input_grid <- expand.grid('b1' = c(2,5,8), 'b3' = c(2,5,8), 'corr' = c(0, 0.2),'q_sd' = c(10, 12))
current_files <- sort(as.numeric(str_remove(list.files('res_Q/'), '.rds')))
all_df <- do.call(rbind,lapply(current_files, function(x){
  df <- read_rds(paste0('res_Q/', x, '.rds')) %>%
    mutate(b1 = input_grid$b1[x], b3 = input_grid$b3[x],
           corr = input_grid$corr[x], q_sd = input_grid$q_sd[x])
  return(df)
}))

# add additional metrics to get cis
all_df <- all_df %>%
  mutate_at(c('dr_norm', 'var_est', 'var_est2', 'var_contrast_est', 'est_y', 'var_contrast_est2'), as.numeric) %>%
  mutate(lb1 = dr_norm - qnorm(0.975)*sqrt(var_est/1000),
         ub1 = dr_norm + qnorm(0.975)*sqrt(var_est/1000),
         lb3 = ifelse(var_contrast_est2 < 0, 0, dr_norm - est_y - qnorm(0.975)*sqrt(var_contrast_est2/1000)),
         ub3 = ifelse(var_contrast_est2 < 0, 0, dr_norm - est_y + qnorm(0.975)*sqrt(var_contrast_est2/1000)),
         diff_est = dr_norm - est_y,
         diff_true = true_mu_q - true_mu) %>%
  mutate(ci1 = ifelse(true_mu_q >= lb1 & true_mu_q <= ub1,1,0),
         ci3 = ifelse(diff_true >=lb3 & diff_true <= ub3,1,0))

# aggregate across all simulated datasets for mu
agg_df <- all_df %>%
  group_by(delta, pve, b1, b3, corr, q_sd, type) %>%
  summarise(est = mean(as.numeric(dr_norm)), true = unique(true_mu_q), emp_sd = sd(as.numeric(dr_norm)), var_est = mean(as.numeric(var_est)), ci1 = mean(ci1), emp_sd_var_est = var(as.numeric(dr_norm) - true_mu_q)*1000, rmse = sqrt(mean((as.numeric(dr_norm) - true_mu_q)^2))) %>%
  mutate(abs_perc_bias = abs(100*(est - true)/true))
write_rds(agg_df, 'saved_results/n1000_q_mu.rds')

# aggregate across all simulated datasets for tau
agg_df2 <- all_df %>%
  group_by(delta, pve, b1, b3, corr, q_sd, type) %>%
  summarise(est = mean(as.numeric(diff_est)), true = unique(diff_true), emp_sd = sd(as.numeric(diff_est)), var_est = mean(as.numeric(var_contrast_est2)), ci1 = mean(ci3), emp_sd_var_est = var(as.numeric(diff_est) - diff_true)*1000, rmse = sqrt(mean((as.numeric(diff_est) - diff_true)^2))) %>%
  mutate(abs_perc_bias = abs(100*(est - true)/true))
write_rds(agg_df2, 'saved_results/n1000_q_tau.rds')


#### Aggregate n = 1000 q*() results ####
# join saved files
input_grid <- expand.grid('b1' = c(2,5,8), 'b3' = c(2,5,8),'q_sd' = c(10, 12), 'b3_coef' = c(0,3), 'b3_sd' = c(1, 3))

all_df3 <- do.call(rbind,lapply(seq(1, 72), function(x){
  df <- read_rds(paste0('res_Qstar/', x, '.rds')) %>%
    mutate(b1 = input_grid$b1[x], b3 = input_grid$b3[x],
           b3_coef = input_grid$b3_coef[x], q_sd = input_grid$q_sd[x], b3_sd = input_grid$b3_sd[x])
  return(df)
}))

# add additional metrics to get cis
all_df3 <- all_df3 %>%
  mutate_at(c('dr_norm', 'var_est', 'var_est2', 'var_contrast_est', 'est_y', 'var_contrast_est2'), as.numeric) %>%
  mutate(lb1 = ifelse(var_est <0, 0, dr_norm - qnorm(0.975)*sqrt(var_est/1000)),
         ub1 = ifelse(var_est<0, 0, dr_norm + qnorm(0.975)*sqrt(var_est/1000)),
         lb3 = ifelse(var_contrast_est2 < 0, 0, dr_norm - est_y - qnorm(0.975)*sqrt(var_contrast_est2/1000)),
         ub3 = ifelse(var_contrast_est2 < 0, 0, dr_norm - est_y + qnorm(0.975)*sqrt(var_contrast_est2/1000)),
         diff_est = dr_norm - est_y,
         diff_true = true_mu_q - true_mu) %>%
  mutate(ci1 = ifelse(true_mu_q >= lb1 & true_mu_q <= ub1,1,0),
         ci3 = ifelse(diff_true >=lb3 & diff_true <= ub3,1,0))


# aggregate across all simulated datasets for mu
agg_df3 <- all_df3 %>%
  group_by(delta,b1, b3, b3_coef, q_sd, b3_sd, type) %>%
  summarise(est = mean(as.numeric(dr_norm)), true = unique(true_mu_q), emp_sd = sd(as.numeric(dr_norm)), var_est = mean(as.numeric(var_est)), ci1 = mean(ci1), emp_sd_var_est = var(as.numeric(dr_norm) - true_mu_q)*1000) %>%
  mutate(abs_perc_bias = abs(100*(est - true)/true))
write_rds(agg_df3, 'saved_results/n1000_qstar_mu.rds')



#### Aggregate varying n results ####

# join saved files
n_options <- c(1000, 2000,4000,8000,16000,32000)
n_df_list <- lapply(seq(1, 240), function(j){
  mult <- trunc((j - 1)/40)
  n_j <- n_options[mult+1]
  sub_df <- read_rds(paste0('res_rate_Q/', j, '.rds')) %>%
    mutate(n = n_j)
  return(sub_df)
})
n_df_all <- do.call(rbind, n_df_list)


# aggregate across all datasets get rmse and cis
n_agg_df <- n_df_all %>%
  ungroup() %>%
  mutate_at(c('dr_norm', 'var_est', 'var_est2', 'var_contrast_est2', 'est_y'), as.numeric) %>%
  mutate(lb1 = dr_norm - qnorm(0.975)*sqrt(var_est/n),
         ub1 = dr_norm + qnorm(0.975)*sqrt(var_est/n),
         lb3 = ifelse(var_contrast_est2 < 0, 0, dr_norm - est_y - qnorm(0.975)*sqrt(var_contrast_est2/n)),
         ub3 = ifelse(var_contrast_est2 < 0, 0, dr_norm - est_y + qnorm(0.975)*sqrt(var_contrast_est2/n)),
         diff_est = dr_norm - est_y,
         diff_true = true_mu_q - true_mu) %>%
  mutate(ci1 = ifelse(true_mu_q >= lb1 & true_mu_q <= ub1,1,0),
         ci3 = ifelse(diff_true >=lb3 & diff_true <= ub3,1,0)) %>%
  group_by(delta, type, pve, n) %>%
  summarise(rmse = sqrt(mean((as.numeric(dr_norm) - true_mu_q)^2)),
            rmse_contrast = sqrt(mean((as.numeric(diff_est) - diff_true)^2)),
            var_est = mean(as.numeric(var_est)),
            var_contrast_est = mean(as.numeric(var_contrast_est2)),
            ci1 = mean(ci1),
            ci3 = mean(ci3),
            true_var_est = var(as.numeric(dr_norm) - true_mu_q)*unique(n),
            true_var_contrast_est = var(as.numeric(diff_est) - diff_true)*unique(n),
            diff_est = mean(diff_est),
            diff_true = unique(diff_true))
write_rds(n_agg_df, 'saved_results/n_vary.rds')



