library(tidyverse)
library(refund)
library(mgcv)
library(tidyfun)
library(patchwork)
library(broom)
library(viridis)
library(splines)
library(splines2)
library(refundr)
library(pracma)
library(LinCDE)
library(caret)
library(randomForest)
source('real_data_helper_funcs.R')

set.seed(779)
all_df <- read_rds('nhanes_real_data_clean.rds') %>%
  mutate(extreme_vals = if_any(c(3:1443), ~.x >60)) %>%
  filter(extreme_vals == FALSE) %>%
  dplyr::select(-extreme_vals)

gamma_x <- function(x, t1, t2){
  (((x-t1)/(t2 - t1))^1.5)*(((t2-x)/(t2 - t1)))
}
t1 <- 60*17
t2 <- 60*20
fpca_df <- all_df %>%
  dplyr::select(seq(1,1441))

fpca_res <- get_coef_fpca_df(fpca_df, t1, t2)
write_rds(fpca_res, 'model_res2/fpca_res_5_8.rds')


mod_df <- fpca_res[[1]] %>%
  left_join(all_df %>%dplyr::select(c(1, 1442:1467)), by = 'SEQN') %>%
  dplyr::select(-SEQN) %>%
  rename(y = mort5yr)

vec1 <- sort(sample(seq(1, nrow(mod_df)), nrow(mod_df)/2, replace = FALSE))
vec2 <- seq(1, nrow(mod_df))[!(seq(1, nrow(mod_df)) %in% vec1)]
mod_df1 <- mod_df[vec1,]
write_rds(mod_df1, 'model_run_data2/mod_df1_5_8.rds')
mod_df2 <- mod_df[vec2,]
write_rds(mod_df2, 'model_run_data2/mod_df2_5_8.rds')

design_pts <- seq(-300, 600, length.out = 750)
hyper_grid <- expand.grid('mtry' = seq(max(3, round((ncol(mod_df)-1)/3) -7) , round((ncol(mod_df)-1)/3) + 5, by = 1))
delta_vec <- c(seq(0.005, 0.05, by = 0.0025))

print('start model fit 1')
dr_lincde1_fpca  <- get_w3_t2(n.trees = 2000, depth = 2, centering = F, centeringMethod = "linearRegression",
                              design_pts, delta_vec, mod_df1, mod_df2, 'coef1_2', hyper_grid, 2002, '1_5_8')
write_rds(dr_lincde1_fpca, 'model_res2/mods1_5_8.rds')
print('end model fit 1')

print('start model fit 2')
dr_lincde2_fpca  <- get_w3_t2(n.trees = 2000, depth = 2,centering = F, centeringMethod = "linearRegression",
                              design_pts, delta_vec, mod_df2, mod_df1, 'coef1_2', hyper_grid, 4022, '2_5_8')
write_rds(dr_lincde2_fpca, 'model_res2/mods2_5_8.rds')
print('end model fit 2')