
get_coef_fpca_df <- function(out, pve = 0.99){
  out_df <- out[[1]]
  func_df <- out[[2]] %>%
    mutate(name = as.numeric(name))
  n <- nrow(out_df)
  t1 <- unique(out_df$t1)
  t2 <- unique(out_df$t2)
  int1 <- seq(1, t1 -1)
  int2 <- seq(t1, t2)
  int3 <- seq(t2+1, 300)
  
  
  fpca_df <- func_df %>%
    dplyr::select(id, name, obs) %>%
    pivot_wider(names_from = name, values_from = obs)
  
  
  ### [0, t1] ###
  fpca_int1 <- fpca.face(as.matrix(fpca_df[,int1+1]), center = FALSE, argvals = int1 + 1, pve = pve)
  func_df_int1_coefs_fpca <- as.data.frame(cbind('id' = seq(1, n), fpca_int1$scores))
  colnames(func_df_int1_coefs_fpca) <- c('id', paste0('coef', seq(1, ncol(fpca_int1$scores))))
  
  
  ### [t1, t2] ###
  fpca_int2 <- fpca.face(as.matrix(fpca_df[,int2+1]), center = FALSE, argvals = int2 + 1, pve = pve)
  
  psi2 <- sqrt(n) * fpca_int2$efunctions
  evalues2 <- fpca_int2$evalues / n
  
  psi2_new <- cbind(gamma_x(int2, t1, t2), psi2)
  new_basis <- gramSchmidt(psi2_new)[[1]]
  colnames(new_basis) <- paste0('x', seq(1, ncol(new_basis)))
  new_basis_big <- replicate(n, new_basis, simplify = FALSE)
  mod_df2_fpca <- func_df %>%
    filter(name %in% int2) %>%
    arrange(id, as.numeric(name)) %>%
    cbind(Reduce(rbind, new_basis_big)) %>%
    dplyr::select(-name, -value, -mean_value)

  func_df_int2_coefs_fpca <- t(sapply(seq(1, n), function(j){
    mod_df2_id <- mod_df2_fpca %>%
      filter(id == j) %>%
      dplyr::select(-id)
    mod <- lm(obs~0+., data = mod_df2_id)
    return(c(j, coef(mod)))
  }))
  
  func_df_int2_coefs_fpca <- as.data.frame(func_df_int2_coefs_fpca)
  colnames(func_df_int2_coefs_fpca) <- c('id', paste0('coef', seq(1,ncol(new_basis)), '_2'))
  
  ### [t2, T] ###
  fpca_int3 <- fpca.face(as.matrix(fpca_df[,int3+1]), center = FALSE, argvals = int3 + 1, pve = pve)

  
  func_df_int3_coefs_fpca <- as.data.frame(cbind('id' = seq(1, n), fpca_int3$scores))
  if(is.null(dim(fpca_int3$scores))){
    colnames(func_df_int3_coefs_fpca) <- c('id', paste0('coef', seq(1, 1), '_3'))
  } else {
    colnames(func_df_int3_coefs_fpca) <- c('id', paste0('coef', seq(1, ncol(fpca_int3$scores)), '_3'))
  }

  
  
  out_df_join <- out_df %>%
    mutate(id = seq(1, unique(n))) %>%
    dplyr::select(id,paste0('x', seq(1, unique(out_df$p))), y)
  mod_df_fpca <- out_df_join %>%
    full_join(func_df_int1_coefs_fpca, by = 'id') %>%
    full_join(func_df_int2_coefs_fpca, by = 'id') %>%
    full_join(func_df_int3_coefs_fpca, by = 'id') %>%
    dplyr::select(-id)
  
  return(mod_df_fpca)
  
}


get_mu <- function(df1, df2, hyper_grid, seed){
  set.seed(seed)
  rf_model1 <- train(x = df1[, -which(colnames(df1) == 'y')],
                     y = df1$y,
                     method = "rf",
                     tuneGrid = hyper_grid,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 2))
  
  hat_mu_2 <- predict(rf_model1, newdata = df2[, -which(colnames(df2) == 'y')])
  return(hat_mu_2)
}



get_w1 <- function(q_sd, design_pts, delta, df1, df2, t_col, hyper_grid, seed){
  set.seed(seed)
  rf_model1 <- train(x = df1[, -which(colnames(df1) == 'y')],
                     y = df1$y,
                     method = "rf",
                     tuneGrid = hyper_grid,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 2))
  rf_model1_y2 <- train(x = df1[, -which(colnames(df1) == 'y')],
                        y = df1$y^2,
                        method = "rf",
                        tuneGrid = hyper_grid,
                        preProcess = c("center", "scale"),
                        trControl = trainControl(method = "cv", number = 2))
  
  rf_df1 <- train(x = df1[, -which(colnames(df1) == 'y' | colnames(df1) == t_col)],
                  y = df1[, t_col],
                  method = "rf",
                  tuneGrid = hyper_grid,
                  preProcess = c("center", "scale"),
                  trControl = trainControl(method = "cv", number = 2))
  hat_pi_2 <- predict(rf_df1, newdata = df2[, -which(colnames(df2) == 'y' | colnames(df2) == t_col)])
  n2 <- nrow(df2)
  
  hat_f_2 <- sapply(hat_pi_2, function(x){
    return(q_func(design_pts, mu = x, sigma = q_sd, delta = 0, lb = min(design_pts), ub = max(design_pts)))
  })
  
  nu2 <- lapply(delta, function(j){
    nu_vec <- sapply(seq(1, n2), function(i){
      mean(exp(j*design_pts)*hat_f_2[,i])*(design_pts[length(design_pts)] - design_pts[1])
    })
    return(nu_vec)
  })
  
  xi2_all <- lapply(seq(1, length(delta)), function(k){
    j <- delta[k]
    hat_mu_2_design <- sapply(design_pts, function(i){
      design_df <- df2 %>%
        dplyr::select(-y) %>%
        mutate(!!sym(t_col) := i)
      return(predict(rf_model1, newdata = design_df))
    })
    
    xi_vec <- sapply(seq(1, n2), function(x){
      mean(hat_mu_2_design[x, ]*hat_f_2[,x]*exp(j*design_pts)/nu2[[k]][x])*(design_pts[length(design_pts)] - design_pts[1])
    })
    
    xi2_vec <- sapply(seq(1, n2), function(x){
      mean((hat_mu_2_design[x, ]^2)*hat_f_2[,x]*exp(j*design_pts)/nu2[[k]][x])*(design_pts[length(design_pts)] - design_pts[1])
    })
    
    return(cbind(xi_vec, xi2_vec))
  })
  
  xi2 <- lapply(xi2_all, function(x){x[,1]})
  xi2_squared <- lapply(xi2_all, function(x){x[,2]})
  
  
  mu_pred <- predict(rf_model1, newdata = df2[, -which(colnames(df2) == 'y')])
  
  y2_pred <- predict(rf_model1_y2, newdata = df2[, -which(colnames(df2) == 'y')])
  
  w <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    return(exp(k*df2[, t_col])/nu2[[j]])
  })
  
  return(list('xi' = xi2, 'w' = w, 'mu_pred' = mu_pred, 'y2_pred' = y2_pred, 'xi_squared' = xi2_squared))
}



get_w3_t2 <- function(n.trees, depth = 2, centering = F, centeringMethod = "linearRegression", design_pts, delta, df1, df2, t_col, hyper_grid, seed){
  set.seed(seed)
  rf_model1 <- train(x = df1[, -which(colnames(df1) == 'y')],
                     y = df1$y,
                     method = "rf",
                     tuneGrid = hyper_grid,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 2))
  rf_model1_y2 <- train(x = df1[, -which(colnames(df1) == 'y')],
                     y = df1$y^2,
                     method = "rf",
                     tuneGrid = hyper_grid,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 2))
  dens_mod <- LinCDE.boost(y = df1[, t_col], X = df1[, -which(colnames(df1) == 'y' | colnames(df1) == t_col)], depth = depth, 
                           n.trees = n.trees, centering = centering, centeringMethod = centeringMethod, numberBin = length(design_pts), minY = min(design_pts), maxY = max(design_pts), verbose = FALSE) 
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
      return(predict(rf_model1, newdata = design_df))
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
  
  mu_pred <- predict(rf_model1, newdata = df2[, -which(colnames(df2) == 'y')])
  
  y2_pred <- predict(rf_model1_y2, newdata = df2[, -which(colnames(df2) == 'y')])
  
  mu_q_random <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    mu_q_vec <- sapply(seq(1, n2), function(x){
      set.seed(j*4 + x*80 +  5556)
      q <- exp(k*design_pts)*hat_pi[x,]/nu2[[j]][x]
      sample(design_pts, size=1, replace=TRUE, prob=q/(sum(q)))
    })
    return(mu_q_vec)
  })
  
  
  w <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    return(exp(k*df2[, t_col])/nu2[[j]])
  })
  
  return(list('xi' = xi2, 'w' = w, 'mu_q_random' = mu_q_random, 'mu_pred' = mu_pred, 'loglik' = loglik, 'y2_pred' = y2_pred, 'xi_squared' = xi2_squared))
  
}


get_w3_t2_fast <- function(n.trees, depth = 2, centering = F, centeringMethod = "linearRegression", design_pts, delta, df1, df2, t_col, hyper_grid, seed, cores){
  set.seed(seed)
  train_n <- min(nrow(df1), 2000)
  train_df1 <- sample_n(df1, train_n)
  
  rf_train1 <- train(x = train_df1[, -which(colnames(df1) == 'y')],
                     y = train_df1$y,
                     method = "ranger",
                     tuneGrid = hyper_grid,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 2, returnData = F))
  
  rf_model1 <- ranger(y~., data = df1, mtry = rf_train1$bestTune$mtry, num.threads = cores)
  
  rf_train1_y2 <- train(x = train_df1[, -which(colnames(df1) == 'y')],
                        y = train_df1$y^2,
                        method = "ranger",
                        tuneGrid = hyper_grid,
                        preProcess = c("center", "scale"),
                        trControl = trainControl(method = "cv", number = 2))
  
  rf_model1_y2 <- ranger(y^2~., data = df1, mtry = rf_train1_y2$bestTune$mtry, num.threads = cores)
  dens_mod <- LinCDE.boost(y = df1[, t_col], X = df1[, -which(colnames(df1) == 'y' | colnames(df1) == t_col)], depth = depth, 
                           n.trees = n.trees, centering = centering, centeringMethod = centeringMethod, numberBin = length(design_pts), minY = min(design_pts), maxY = max(design_pts), verbose = FALSE) 
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
      return(predict(rf_model1, data = design_df)$predictions)
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

  
  mu_pred <- predict(rf_model1, data = df2[, -which(colnames(df2) == 'y')])$predictions
  
  y2_pred <- predict(rf_model1_y2, data = df2[, -which(colnames(df2) == 'y')])$predictions
  
  mu_q_random <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    mu_q_vec <- sapply(seq(1, n2), function(x){
      set.seed(j*4 + x*80 +  5556)
      q <- exp(k*design_pts)*hat_pi[x,]/nu2[[j]][x]
      sample(design_pts, size=1, replace=TRUE, prob=q/(sum(q)))
    })
    return(mu_q_vec)
  })
  
  
  w <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    return(exp(k*df2[, t_col])/nu2[[j]])
  })
  
  return(list('xi' = xi2, 'w' = w, 'mu_q_random' = mu_q_random, 'mu_pred' = mu_pred, 'loglik' = loglik, 'y2_pred' = y2_pred, 'xi_squared' = xi2_squared))
  
}




get_w3_star <- function(n.trees, depth = 2, centering = F, centeringMethod = "linearRegression", design_pts, delta, df1, df2, t_col, hyper_grid, seed){
  set.seed(seed)
  rf_model1 <- train(x = df1[, -which(colnames(df1) == 'y')],
                     y = df1$y,
                     method = "rf",
                     tuneGrid = hyper_grid,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 2))
  rf_model1_y2 <- train(x = df1[, -which(colnames(df1) == 'y')],
                        y = df1$y^2,
                        method = "rf",
                        tuneGrid = hyper_grid,
                        preProcess = c("center", "scale"),
                        trControl = trainControl(method = "cv", number = 2))
  
  time3_cols <- colnames(df1)[str_detect(colnames(df1), '_3')]
  f_df1 <- df1 %>%
      dplyr::select(-all_of(time3_cols), - y)
  f_df2 <- df2 %>%
      dplyr::select(-all_of(time3_cols), - y)

  dens_mod <- LinCDE.boost(y = df1[, t_col], X = f_df1[, -which(colnames(f_df1) == t_col)], depth = depth,
                             n.trees = n.trees, centering = centering, centeringMethod = centeringMethod, numberBin = length(design_pts), minY = min(design_pts), maxY = max(design_pts), verbose = FALSE)
  pred <- predict(dens_mod, X = f_df2[, -which(colnames(f_df2) == t_col)], densityOnly = FALSE)
  pred_y <- predict(dens_mod, X = f_df2[, -which(colnames(f_df2) == t_col)], y = f_df2[, t_col], densityOnly = FALSE)
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
      return(predict(rf_model1, newdata = design_df))
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

  mu_pred <- predict(rf_model1, newdata = df2[, -which(colnames(df2) == 'y')])
  
  y2_pred <- predict(rf_model1_y2, newdata = df2[, -which(colnames(df2) == 'y')])
  
  mu_q_random <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    mu_q_vec <- sapply(seq(1, n2), function(x){
      set.seed(j*4 + x*80 +  5556)
      q <- exp(k*design_pts)*hat_pi[x,]/nu2[[j]][x]
      sample(design_pts, size=1, replace=TRUE, prob=q/(sum(q)))
    })
    return(mu_q_vec)
  })
  
  
  w <- lapply(seq(1, length(delta)), function(j){
    k <- delta[j]
    return(exp(k*df2[, t_col])/nu2[[j]])
  })
  
  return(list('xi' = xi2, 'w' = w, 'mu_q_random' = mu_q_random, 'mu_pred' = mu_pred, 'loglik' = loglik, 'y2_pred' = y2_pred, 'xi_squared' = xi2_squared))
  
}



get_est <- function(res1, res2, y1, y2, delta){
  #2 and 1 are switched to match with y (consistent with the cross-fitting)
  w1 <- res2$w
  w2 <- res1$w
  xi1 <- res2$xi
  xi2 <- res1$xi
  mu1 <- res2$mu_pred
  mu2 <- res1$mu_pred
  y2_1 <- res2$y2_pred
  y2_2 <- res1$y2_pred
  xi1_squared <- res2$xi_squared
  xi2_squared <- res1$xi_squared
  
  k <- length(delta)
  out_mod <- sapply(seq(1, k), function(j){
    mean(c(xi1[[j]], xi2[[j]]))
  })
  w_norm <- sapply(seq(1, k), function(j){
    mean(c(w1[[j]]*y1/mean(w1[[j]]), w2[[j]]*y2/mean(w2[[j]])))
  })
  
  
  dr_norm <- sapply(seq(1, k), function(j){
    mean(c(w1[[j]]*(y1-xi1[[j]])/mean(w1[[j]])+ xi1[[j]], w2[[j]]*(y2-xi2[[j]])/mean(w2[[j]]) + xi2[[j]]))
  })
  
  dr_norm2 <- sapply(seq(1, k), function(j){
    mean(c(w1[[j]]*(y1-mu1)/mean(w1[[j]])+ xi1[[j]], w2[[j]]*(y2-mu2)/mean(w2[[j]]) + xi2[[j]]))
  })
  
  var_est <- sapply(seq(1, k), function(j){
    mean(c(((w1[[j]]/mean(w1[[j]]))^2)*((y2_1 - mu1^2) + (mu1 - xi1[[j]])^2), ((w2[[j]]/mean(w2[[j]]))^2)*((y2_2 - mu2^2) + (mu2 - xi2[[j]])^2))) + var(c(xi1[[j]], xi2[[j]]))
  })
  var_est2 <- sapply(seq(1, k), function(j){
    mean(c(((w1[[j]]/mean(w1[[j]]))^2)*((y1 - mu1)^2 + (mu1 - xi1[[j]])^2), ((w2[[j]]/mean(w2[[j]]))^2)*((y2 - mu2)^2 + (mu2 - xi2[[j]])^2))) + var(c(xi1[[j]], xi2[[j]]))
  })
  var_contrast_est <- sapply(seq(1, k), function(j){
    mean(c(((w1[[j]]/mean(w1[[j]]))^2 - 1)*(y2_1 - mu1^2) + ((w1[[j]]/mean(w1[[j]]))^2)*((mu1 - xi1[[j]])^2), ((w2[[j]]/mean(w2[[j]]))^2 - 1)*(y2_2 - mu2^2) + ((w2[[j]]/mean(w2[[j]]))^2)*((mu2 - xi2[[j]])^2))) + var(c(mu1, mu2) - c(xi1[[j]], xi2[[j]])) - 2*(var(c(mu1, mu2))-var(c(xi1[[j]], xi2[[j]])))
  })
  
  var_contrast_est2 <- sapply(seq(1, k), function(j){
    mean(c(((w1[[j]]/mean(w1[[j]]))^2 - 1)*(y2_1 - mu1^2) + ((w1[[j]]/mean(w1[[j]]))^2)*((mu1 - xi1[[j]])^2), ((w2[[j]]/mean(w2[[j]]))^2 - 1)*(y2_2 - mu2^2) + ((w2[[j]]/mean(w2[[j]]))^2)*((mu2 - xi2[[j]])^2))) + var(c(mu1, mu2) - c(xi1[[j]], xi2[[j]])) - 2*mean(c(xi1_squared[[j]], xi2_squared[[j]])-c(xi1[[j]], xi2[[j]])^2)
  })
  
  return(cbind('delta' = delta, 'out_mod' = out_mod,  'w_norm' = w_norm,'dr_norm' = dr_norm, 
               'dr_norm2' = dr_norm2, 'var_est' = var_est, 'var_est2' = var_est2, 'var_contrast_est' = var_contrast_est, 'var_contrast_est2' = var_contrast_est2))
}



get_est_test <- function(res1, res2, y1, y2, delta){
  w1 <- res2$w
  w2 <- res1$w
  xi1 <- res2$xi
  xi2 <- res1$xi
  mu1 <- res2$mu_pred
  mu2 <- res1$mu_pred
  y2_1 <- res2$y2_pred
  y2_2 <- res1$y2_pred
  
  k <- length(delta)
  out_mod <- sapply(seq(1, k), function(j){
    mean(c(xi1[[j]], xi2[[j]]))
  })

  w_norm <- sapply(seq(1, k), function(j){
    mean(c(w1[[j]]*y1/mean(w1[[j]]), w2[[j]]*y2/mean(w2[[j]])))
  })

  
  dr_norm <- sapply(seq(1, k), function(j){
    mean(c(w1[[j]]*(y1-xi1[[j]])/mean(w1[[j]])+ xi1[[j]], w2[[j]]*(y2-xi2[[j]])/mean(w2[[j]]) + xi2[[j]]))
  })
  
  dr_norm2 <- sapply(seq(1, k), function(j){
    mean(c(w1[[j]]*(y1-mu1)/mean(w1[[j]])+ xi1[[j]], w2[[j]]*(y2-mu2)/mean(w2[[j]]) + xi2[[j]]))
  })
  
  var_est <- sapply(seq(1, k), function(j){
    mean(c(((w1[[j]]/mean(w1[[j]]))^2)*((y2_1 - mu1^2) + (mu1 - xi1[[j]])^2), ((w2[[j]]/mean(w2[[j]]))^2)*((y2_2 - mu2^2) + (mu2 - xi2[[j]])^2))) + var(c(xi1[[j]], xi2[[j]]))
  })
  
  xi_mean <- sapply(seq(1, k), function(j){
    mean(c(xi1[[j]], xi2[[j]]))
  })
  
  mu_mean <- rep(mean(c(mu1, mu2)), k)

  return(cbind('delta' = delta, 'out_mod' = out_mod,  'w_norm' = w_norm,'dr_norm' = dr_norm, 
               'dr_norm2' = dr_norm2, 'var_est' = var_est, 'muXA_mean' = mu_mean, 'EQmean' = xi_mean))
}


### variance estimation functions!! ####

get_var_est <- function(res1, res2, y1, y2, delta){
  w1 <- res2$w
  w2 <- res1$w
  xi1 <- res2$xi
  xi2 <- res1$xi
  mu1 <- res2$mu_pred
  mu2 <- res1$mu_pred
  y2_1 <- res2$y2_pred
  y2_2 <- res1$y2_pred
  
  k <- length(delta)
  
  var_est <- sapply(seq(1, k), function(j){
    mean(c(((w1[[j]]/mean(w1[[j]]))^2)*((y2_1 - mu1^2) + (mu1 - xi1[[j]])^2), ((w2[[j]]/mean(w2[[j]]))^2)*((y2_2 - mu2^2) + (mu2 - xi2[[j]])^2))) + var(c(xi1[[j]], xi2[[j]]))
  })
  
  return(cbind('delta' = delta, 'var_est' = var_est))
}







