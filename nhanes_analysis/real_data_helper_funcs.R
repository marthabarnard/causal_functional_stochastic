get_coef_fpca_df <- function(fpca_df, t1, t2, pve = 0.99){
  n <- nrow(fpca_df)
  int1 <- seq(1, t1 -1)
  int2 <- seq(t1, t2)
  int3 <- seq(t2+1, 1440)
  
  ### [0, t1] ###
  fpca_int1 <- fpca.face(as.matrix(fpca_df[,int1+1]), center = TRUE, argvals = int1+1, pve = pve)
  func_df_int1_coefs_fpca <- as.data.frame(cbind(fpca_df$SEQN, fpca_int1$scores))
  colnames(func_df_int1_coefs_fpca) <- c('SEQN', paste0('coef', seq(1, ncol(fpca_int1$scores))))
  
  
  ### [t1, t2] ###
  fpca_int2 <- fpca.face(as.matrix(fpca_df[,int2+1]), center = TRUE, argvals = int2+1, pve = pve)
  
  
  psi2 <- fpca_int2$efunctions
  
  psi2_new <- cbind(gamma_x(int2, t1, t2), psi2)
  new_basis <- gramSchmidt(psi2_new)[[1]]
  colnames(new_basis) <- paste0('x', seq(1, ncol(new_basis)))
  #new_basis_big <- replicate(n, new_basis, simplify = FALSE)
  # mod_df2_fpca <- fpca_df[,c(1, int2+1)] %>%
  #   pivot_longer(cols = as.character(int2)) %>%
  #   arrange(SEQN, as.numeric(name)) %>%
  #   cbind(Reduce(rbind, new_basis_big)) %>%
  #   dplyr::select(-name)
  
  func_df_int2_fpca <- lapply(seq(1, n), function(j){
    mod_df2_id <- data.frame(cbind('obs' = unlist(fpca_df[j, int2+1]) - fpca_int2$mu, new_basis))
    mod <- lm(obs~0+., data = mod_df2_id)
    # plot(new_basis %*% coef(mod) + fpca_int2$mu)
    # plot(unlist(fpca_df[j, int2+1]))
    return(list(c(fpca_df$SEQN[j], coef(mod)), new_basis %*% coef(mod) + fpca_int2$mu))
  })
  
  func_df_int2_coefs_fpca <- t(sapply(func_df_int2_fpca, function(x){x[[1]]}))
  func_df_int2_pred_fpca <- t(sapply(func_df_int2_fpca, function(x){x[[2]]}))
  
  func_df_int2_coefs_fpca <- as.data.frame(func_df_int2_coefs_fpca)
  colnames(func_df_int2_coefs_fpca) <- c('SEQN', paste0('coef', seq(1,ncol(new_basis)), '_2'))
  
  ### [t2, T] ###
  fpca_int3 <- fpca.face(as.matrix(fpca_df[,int3+1]), center = TRUE, argvals = int3+1, pve = pve)
  
  
  func_df_int3_coefs_fpca <- as.data.frame(cbind('SEQN' = fpca_df$SEQN, fpca_int3$scores))
  colnames(func_df_int3_coefs_fpca) <- c('SEQN', paste0('coef', seq(1, ncol(fpca_int3$scores)), '_3'))

  
  mean_list <- list('int1' = fpca_int1$mu, 'int2' = fpca_int2$mu, 'int3' = fpca_int3$mu)
  func_pred <- cbind(fpca_int1$Yhat, func_df_int2_pred_fpca, fpca_int3$Yhat)
  
  
  mod_df_fpca <- func_df_int1_coefs_fpca %>%
    full_join(func_df_int2_coefs_fpca, by = 'SEQN') %>%
    full_join(func_df_int3_coefs_fpca, by = 'SEQN') 
  
  return(list('mod_df' = mod_df_fpca, 'mean_list' = mean_list, 'func_pred_df' = func_pred))
  
}


get_w3_t2 <- function(n.trees, depth = 2, centering = F, centeringMethod = "linearRegression", design_pts, delta, df1, df2, t_col, hyper_grid, seed, mod_name){
  set.seed(seed)
  rf_model1 <- train(x = df1[, -which(colnames(df1) == 'y')],
                     y = as.factor(df1$y),
                     method = "rf",
                     tuneGrid = hyper_grid,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 2))
  # rf_model1_y2 <- train(x = df1[, -which(colnames(df1) == 'y')],
  #                       y = df1$y^2,
  #                       method = "rf",
  #                       tuneGrid = hyper_grid,
  #                       preProcess = c("center", "scale"),
  #                       trControl = trainControl(method = "cv", number = 2))
  #splineDf = 20, prior = "uniform",basis = "nsTransform", 
  dens_mod <- LinCDE.boost(y = df1[, t_col], X = df1[, -which(colnames(df1) == 'y' | colnames(df1) == t_col)], depth = depth, 
                           n.trees = n.trees, centering = centering, centeringMethod = centeringMethod, numberBin = length(design_pts), minY = min(design_pts), maxY = max(design_pts), verbose = FALSE)
  #write_rds(dens_mod, paste0('model_res/density_mod_res_', mod_name, '.rds'))
  write_rds(dens_mod, paste0('model_res2/density_mod_res_', mod_name, '.rds'))
  pred <- predict(dens_mod, X = df2[, -which(colnames(df2) == 'y' | colnames(df2) == t_col)], densityOnly = FALSE)
  pred_y <- predict(dens_mod, X = df2[, -which(colnames(df2) == 'y' | colnames(df2) == t_col)], y = df2[, t_col], densityOnly = FALSE)
  hat_pi <- pred[[1]]
  
  loglik <- max(pred_y$testLogLikelihoodHistory)
  
  design_pts <-pred[[2]]
  n2 <- nrow(df2)
  
  nu2 <- lapply(delta, function(j){
    nu_vec <- unname(apply(hat_pi, 1, function(x){
      mean(exp(j*design_pts)*x)*(design_pts[length(design_pts)] - design_pts[1])
    }))
    return(nu_vec)
  })
  
  xi2_all <- lapply(seq(1, length(delta)), function(k){
    j <- delta[k]
    hat_mu_2_design <- sapply(design_pts, function(i){
      design_df <- df2 %>%
        dplyr::select(-y) %>%
        mutate(!!sym(t_col) := i)
      return(predict(rf_model1, newdata = design_df,  type = "prob")$`1`)
    })
    
    xi_vec <- sapply(seq(1, n2), function(x){
      mean(hat_mu_2_design[x, ]*hat_pi[x,]*exp(j*design_pts)/nu2[[k]][x])*(design_pts[length(design_pts)] - design_pts[1])
    })
    
    xi2_vec <- sapply(seq(1, n2), function(x){
      mean((hat_mu_2_design[x, ]^2)*hat_pi[x,]*exp(j*design_pts)/nu2[[k]][x])*(design_pts[length(design_pts)] - design_pts[1])
    })
    
    return(cbind(xi_vec, xi2_vec))
  })
  
  xi2 <- lapply(xi2_all, function(x){x[,1]})
  xi2_squared <- lapply(xi2_all, function(x){x[,2]})
  
  
  mu_pred <- predict(rf_model1, newdata = df2[, -which(colnames(df2) == 'y')], type = "prob")$`1`
  
  y2_pred <- mu_pred #predict(rf_model1_y2, newdata = df2[, -which(colnames(df2) == 'y')])
  
  mu_q_random <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    mu_q_vec <- sapply(seq(1, n2), function(x){
      set.seed(j*4 + x*80 +  5556)
      q <- exp(k*design_pts)*hat_pi[x,]/nu2[[j]][x]
      sample(design_pts, size=1, replace=TRUE, prob=q/(sum(q)))
    })
    return(mu_q_vec)
  })
  
  mu_q_mean <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    mu_q_vec <- sapply(seq(1, n2), function(x){
      sum(design_pts*exp(k*design_pts)*hat_pi[x,]/nu2[[j]][x])
    })
    return(mu_q_vec)
  })
  
  
  w <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    return(exp(k*df2[, t_col])/nu2[[j]])
  })
  
  return(list('xi' = xi2, 'w' = w, 'mu_q_random' = mu_q_random,'mu_q_mean' = mu_q_mean, 'mu_pred' = mu_pred, 'loglik' = loglik, 'y2_pred' = y2_pred, 'xi_squared' = xi2_squared))
  
}

