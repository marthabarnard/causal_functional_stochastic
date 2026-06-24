q_func <- function(a, mu, sigma, delta, lb, ub){
  dtruncnorm(a, lb, ub, mu, sigma)*exp(a*delta)
}

gamma_x <- function(x, t1, t2){
  50*(((x-t1)/(t2 - t1))^4)*(((t2-x)/(t2 - t1))^4)
}


gamma_x2 <- function(x, t2){
  5*(((x-t2)/(300 - t2))^4)*(((300-x)/(300 - t2))^1)
}

gamma_x1 <- function(x, t1){
  (((x)/(t1))^3)*(((t1-x)/(t1))^1)
}

algorithm_gen_coefs <- function(start_val, p, shuffle = T){
  res_vec <- vector(mode ='numeric', length = p)
  res_vec[1] <- start_val
  if(p == 1){
    return(res_vec)
  } else{
    for(j in seq(2, p)){
      val <- 2*(j/p)*start_val -3*(j/p)*start_val^2
      if(j %% 3 == 1){
        res_vec[j] <- -val
      } else{
        res_vec[j] <- val
      }
    }
    if(shuffle == F){
      return(res_vec*(10/max(10, p)))
    }
    if(shuffle == T){
      vals1 <- sort(seq(1, p, by = 2), decreasing = T)
      vals2 <- seq(2, p, by = 2)
      return((10/max(10, p))*res_vec[c(vals1, vals2)])
    }}
  
}


gen_funcs_and_coefs <- function(seed, t1, t2, p, b1, b2, b3){
  set.seed(seed)
  t1 <- 150
  t2 <- 200
  int1 <- seq(1, t1 -1)
  int2 <- seq(t1, t2)
  int3 <- seq(t2+1, 300)
 
  func1 <- cbind(sapply(seq(0, b1-1), function(x){gamma_x1(int1, t1-1)*int1^x}))
  func1 <- gramSchmidt(func1)$Q
  func2 <- cbind(sapply(seq(0, b2-1), function(x){gamma_x(int2, t1, t2)*int2^x}))
  func2 <- gramSchmidt(func2)$Q
  func3 <- cbind(sapply(seq(0, b3-1), function(x){gamma_x2(int3,t2+1)*int3^x}))
  func3 <- gramSchmidt(func3)$Q
  
  vcov1 <- apply(diag(rep(1,p)),2, function(x){ifelse(x == 0, 0, x)})
  vcov2 <- apply(diag(rep(1,b1)),2, function(x){ifelse(x == 0, 0, x)})
  vcov3 <- apply(diag(rep(1,b3)),2, function(x){ifelse(x == 0, 0, x)})
  #each row of below is the coeficients for the basis function (i.e., row 1 is the lin comb coefs for basis function 1)
  func1_lincomb <- rtmvnorm(b1, rep(0, p), vcov1, lower = rep(-2.5, p), upper = rep(2.5, p)) 
  
  func3_lincomb <- rtmvnorm(b3, rep(0,b1), vcov2, lower = rep(-2.5,  b1), upper = rep(2.5, b1))
  
  func2_1_lincomb <- rtruncnorm(b3, -3, 3, 0, 2)
  func2_lincomb <- rbind(func2_1_lincomb,
                         rtmvnorm(b2-1, rep(0, b3), vcov3, lower = rep(-1, b3), upper = rep(1, b3))) 
  rownames(func2_lincomb) <- NULL
  return(list('func1' = func1, 'func2' = func2, 'func3' = func3,
              'func1_lincomb' = func1_lincomb, 'func3_lincomb' = func3_lincomb, 'func2_lincomb' = func2_lincomb))
}



gen_sim <- function(seed,
                    n,
                    p,
                    b1,
                    b2,
                    b3,
                    t1,
                    t2,
                    corr,
                    q_sd,
                    func_coef_out,
                    lb = NULL,
                    ub = NULL){
  
  ### pull results from func_coef_out ###
  func1 <- func_coef_out[['func1']]
  func3 <- func_coef_out[['func3']]
  func2 <- func_coef_out[['func2']]
  func1_lincomb <- func_coef_out[['func1_lincomb']]
  func3_lincomb <- func_coef_out[['func3_lincomb']]
  func2_lincomb <- func_coef_out[['func2_lincomb']]
  
  
  ### generate covariates ###
  set.seed(seed)
  x0 <- rep(1, n) #intercept
  vcov <- apply(diag(rep(1,p)),2, function(x){ifelse(x == 0, corr, x)})
  v1_vp <- rtmvnorm(n, rep(0, p), vcov, lower = rep(-2.5, p), upper = rep(2.5, p)) 
  
  x1_x2 <- v1_vp[, 1:round(p/2)]
  x3_x4 <- apply(v1_vp[,(round(p/2) +1):p], 2, function(x){ifelse(x < 0, 1, 0)})
  x <- cbind(x0, x1_x2, x3_x4)
  colnames(x) <- paste0('x', seq(0, p))
  
  ### set up time and functional basis ###
  int1 <- seq(1, t1 -1)
  int2 <- seq(t1, t2)
  int3 <- seq(t2+1, 300)
  
  
  ps_coef <- seq(0.5, 3, length.out =p)
  b1_coef <- seq(0.5, 3, length.out =b1)/20
  b3_coef <- seq(0.5, 3, length.out =b3)
  
  func_int1_coefs <- sapply(seq(1, b1), function(j){
    coefs_j <- ps_coef*func1_lincomb[j, ]
    return(x%*%c(0, coefs_j))
  })
  
  #compute functions and add error to coefficients
  funcs_int1_all <- lapply(seq(1, nrow(func_int1_coefs)), function(j){
    i <- func_int1_coefs[j, ]
    error1 <- rtruncnorm(b1, -1, 1, 0, 1)
    coefs <- i + error1
    f <- rowSums(sapply(seq(1, b1), function(x){
      return(coefs[x]*func1[int1,x])
    }))
    return(list(f, coefs))
  })
  
  funcs_int1 <- t(sapply(funcs_int1_all, function(x){x[[1]]}))
  int1_coefs <- t(sapply(funcs_int1_all, function(x){x[[2]]}))
  
  ### set up [t2, T] ###
  func_int3_coefs <- sapply(seq(1, b3), function(j){
    coefs_j <- c(b1_coef)*func3_lincomb[j, ]
    return(cbind(int1_coefs)%*%c(coefs_j))
  })
  
  funcs_int3_all <- lapply(seq(1, nrow(func_int3_coefs)), function(j){
    i <- func_int3_coefs[j, ]
    error1 <- rtruncnorm(b3, -1, 1, 0, 1)
    coefs <- i + error1
    f <- rowSums(sapply(seq(1, b3), function(x){
      return(coefs[x]*func3[,x])
    }))
    return(list(f, coefs))
  })
  
  
  funcs_int3 <- t(sapply(funcs_int3_all, function(x){x[[1]]}))
  int3_coefs <- t(sapply(funcs_int3_all, function(x){x[[2]]}))
  
  ### set up [t1, t2] function - generate basis coefficients given X, [0, t1], and [t2, T]###
  func_int2_coefs <- cbind(sapply(seq(1, b2), function(j){
      coefs_j <- c(b3_coef)*func2_lincomb[j, ]
    return(cbind(int3_coefs)%*%c(coefs_j))
  }))
  
  
  if(is.null(lb) == T){
    lb <- floor(range(func_int2_coefs[, 1])[1]) - 40
  }
  if(is.null(ub) == T){
    ub <- ceiling(range(func_int2_coefs[, 1])[2]) + 40
  }
  
  funcs_int2_all <- lapply(seq(1, nrow(func_int2_coefs)), function(j){
    i <- func_int2_coefs[j, ]
    error2 <- c(0, rtruncnorm(b2-1, -1, 1, 0, 1))
    coefs <- i + error2
    coefs[1] <- rtruncnorm(1, lb, ub, coefs[1], q_sd)
    f <- rowSums(sapply(seq(1,b2), function(x){
      return(coefs[x]*func2[,x])
    }))
    return(list(f, coefs))
  })
  
  funcs_int2 <- t(sapply(funcs_int2_all, function(x){x[[1]]}))
  int2_coefs <- t(sapply(funcs_int2_all, function(x){x[[2]]}))
  
  
  func_df <- data.frame(cbind(funcs_int1, funcs_int2, funcs_int3))
  colnames(func_df) <- c(int1, int2, int3)
  func_df <- func_df %>%
    mutate(id = seq(1, n)) %>%
    pivot_longer(cols = (c(int1, int2, int3)))
  
  funcs_int2_mean <- sapply(seq(1, nrow(func_int2_coefs)), function(j){
    i <- func_int2_coefs[j, ]
    coefs <- i
    f <- rowSums(sapply(seq(1, b2), function(x){
      return(coefs[x]*func2[,x])
    }))
    return(f)
  })
  
  
  
  func_df$mean_value <- as.vector(rbind(t(funcs_int1), funcs_int2_mean, t(funcs_int3)))
  

  norm_error <- rnorm(nrow(func_df), 0, 0.25)
  
  func_df$obs <- func_df$value + norm_error
  
  
  
  ### generate outcomes ###
  x_coefs <- c(4, algorithm_gen_coefs(0.03, p))
  t1_coefs <- algorithm_gen_coefs(0.04, b1, shuffle = F)
  t2_coef1 <-  algorithm_gen_coefs(0.03, b2-1, shuffle = T)
  t2_coef2 <-  algorithm_gen_coefs(0.05, b2-1, shuffle = T)
  t3_coefs <- algorithm_gen_coefs(0.02, b3, shuffle = F)
  q_coef1 <- 2*abs(max(c(t1_coefs, t2_coef1, t3_coefs)))
  q_coef2 <- 2*abs(max(c(t1_coefs, t2_coef2, t3_coefs)))
  t2_coefs1 <- c(q_coef1, t2_coef1)
  t2_coefs2 <- c(q_coef2, t2_coef2)
  
  
  mu_pt1 <- x %*% x_coefs + int1_coefs %*% t1_coefs + int3_coefs%*% t3_coefs
  mu <- ifelse(mu_pt1 < 0, mu_pt1 + int2_coefs %*% t2_coefs1,  mu_pt1+int2_coefs %*% t2_coefs2)
  y <- mu + rnorm(n, 0, 0.5)
  
  all_coefs <- cbind(int1_coefs, int2_coefs, int3_coefs)
  colnames(all_coefs) <- c(paste0('coef', seq(1, b1)), paste0('coef', seq(1, b2), '_2'), paste0('coef', seq(1, b3), '_3'))
  
  only_t2_coefs1 <- matrix(rep(c(q_coef1, t2_coef1), n), nrow = n, byrow = T)
  colnames(only_t2_coefs1) <- paste0('t2_mu1_', seq(1, b2))
  
  only_t2_coefs2 <- matrix(rep(c(q_coef2, t2_coef2), n), nrow = n, byrow = T)
  colnames(only_t2_coefs2) <- paste0('t2_mu2_', seq(1, b2))
  
  out_df <- data.frame('seed' = seed, 'n' = n, 'p' =p , 'b1' = b1, 'b2' = b2, 'b3' = b3,  'corr' = corr,
                       'q_sd' = q_sd, 't1' = t1, 't2' = t2, 'lb' = lb, 'ub' = ub, x,
                       all_coefs, 'mean_coef1_2' = func_int2_coefs[,1],
                       only_t2_coefs1, only_t2_coefs2, 'mu_pt1' = mu_pt1, 'mu' = mu, 'y' = y)
  
  
  return(list('out' = out_df, 'func_df' = func_df))
}



gen_sim_no_func <- function(seed,
                    n,
                    p,
                    b1,
                    b2,
                    b3,
                    t1,
                    t2,
                    corr,
                    q_sd,
                    func_coef_out,
                    lb = NULL,
                    ub = NULL){
  
  ### pull results from func_coef_out ###
  func1 <- func_coef_out[['func1']]
  func3 <- func_coef_out[['func3']]
  func2 <- func_coef_out[['func2']]
  func1_lincomb <- func_coef_out[['func1_lincomb']]
  func3_lincomb <- func_coef_out[['func3_lincomb']]
  func2_lincomb <- func_coef_out[['func2_lincomb']]
  
  
  ### generate covariates ###
  set.seed(seed)
  x0 <- rep(1, n) #intercept
  vcov <- apply(diag(rep(1,p)),2, function(x){ifelse(x == 0, corr, x)})
  v1_vp <- rtmvnorm(n, rep(0, p), vcov, lower = rep(-2.5, p), upper = rep(2.5, p)) 
  
  x1_x2 <- v1_vp[, 1:round(p/2)]
  x3_x4 <- apply(v1_vp[,(round(p/2) +1):p], 2, function(x){ifelse(x < 0, 1, 0)})
  x <- cbind(x0, x1_x2, x3_x4)
  colnames(x) <- paste0('x', seq(0, p))
  
  ### set up time and functional basis ###
  int1 <- seq(1, t1 -1)
  int2 <- seq(t1, t2)
  int3 <- seq(t2+1, 300)
  
  
  ps_coef <- seq(0.5, 3, length.out =p)
  b1_coef <- seq(0.5, 3, length.out =b1)/20
  b3_coef <- seq(0.5, 3, length.out =b3)
  
  func_int1_coefs <- sapply(seq(1, b1), function(j){
    coefs_j <- ps_coef*func1_lincomb[j, ]
    return(x%*%c(0, coefs_j))
  })
  
  #compute functions and add error to coefficients
  funcs_int1_all <- lapply(seq(1, nrow(func_int1_coefs)), function(j){
    i <- func_int1_coefs[j, ]
    error1 <- rtruncnorm(b1, -1, 1, 0, 1)
    coefs <- i + error1
    f <- rowSums(sapply(seq(1, b1), function(x){
      return(coefs[x]*func1[int1,x])
    }))
    return(list(f, coefs))
  })
  
  funcs_int1 <- t(sapply(funcs_int1_all, function(x){x[[1]]}))
  int1_coefs <- t(sapply(funcs_int1_all, function(x){x[[2]]}))
  
  ### set up [t2, T] ###
  func_int3_coefs <- sapply(seq(1, b3), function(j){
    coefs_j <- c(b1_coef)*func3_lincomb[j, ]
    return(cbind(int1_coefs)%*%c(coefs_j))
  })
  
  funcs_int3_all <- lapply(seq(1, nrow(func_int3_coefs)), function(j){
    i <- func_int3_coefs[j, ]
    error1 <- rtruncnorm(b3, -1, 1, 0, 1)
    coefs <- i + error1
    f <- rowSums(sapply(seq(1, b3), function(x){
      return(coefs[x]*func3[,x])
    }))
    return(list(f, coefs))
  })
  
  
  funcs_int3 <- t(sapply(funcs_int3_all, function(x){x[[1]]}))
  int3_coefs <- t(sapply(funcs_int3_all, function(x){x[[2]]}))
  
  ### set up [t1, t2] function - generate basis coefficients given X, [0, t1], and [t2, T]###
  func_int2_coefs <- cbind(sapply(seq(1, b2), function(j){
    coefs_j <- c(b3_coef)*func2_lincomb[j, ]
    return(cbind(int3_coefs)%*%c(coefs_j))
  }))
  
  
  if(is.null(lb) == T){
    lb <- floor(range(func_int2_coefs[, 1])[1]) - 40
  }
  if(is.null(ub) == T){
    ub <- ceiling(range(func_int2_coefs[, 1])[2]) + 40
  }
  
  funcs_int2_all <- lapply(seq(1, nrow(func_int2_coefs)), function(j){
    i <- func_int2_coefs[j, ]
    error2 <- c(0, rtruncnorm(b2-1, -1, 1, 0, 1))
    coefs <- i + error2
    coefs[1] <- rtruncnorm(1, lb, ub, coefs[1], q_sd)
    f <- rowSums(sapply(seq(1,b2), function(x){
      return(coefs[x]*func2[,x])
    }))
    return(list(f, coefs))
  })
  
  funcs_int2 <- t(sapply(funcs_int2_all, function(x){x[[1]]}))
  int2_coefs <- t(sapply(funcs_int2_all, function(x){x[[2]]}))
  

  
  funcs_int2_mean <- sapply(seq(1, nrow(func_int2_coefs)), function(j){
    i <- func_int2_coefs[j, ]
    coefs <- i
    f <- rowSums(sapply(seq(1, b2), function(x){
      return(coefs[x]*func2[,x])
    }))
    return(f)
  })
  
  
  ### generate outcomes ###
  x_coefs <- c(4, algorithm_gen_coefs(0.03, p))
  t1_coefs <- algorithm_gen_coefs(0.04, b1, shuffle = F)
  t2_coef1 <-  algorithm_gen_coefs(0.03, b2-1, shuffle = T)
  t2_coef2 <-  algorithm_gen_coefs(0.05, b2-1, shuffle = T)
  t3_coefs <- algorithm_gen_coefs(0.02, b3, shuffle = F)
  q_coef1 <- 2*abs(max(c(t1_coefs, t2_coef1, t3_coefs)))
  q_coef2 <- 2*abs(max(c(t1_coefs, t2_coef2, t3_coefs)))
  t2_coefs1 <- c(q_coef1, t2_coef1)
  t2_coefs2 <- c(q_coef2, t2_coef2)
  
  
  mu_pt1 <- x %*% x_coefs + int1_coefs %*% t1_coefs + int3_coefs%*% t3_coefs
  mu <- ifelse(mu_pt1 < 0, mu_pt1 + int2_coefs %*% t2_coefs1,  mu_pt1+int2_coefs %*% t2_coefs2)
  y <- mu + rnorm(n, 0, 0.5)
  
  all_coefs <- cbind(int1_coefs, int2_coefs, int3_coefs)
  colnames(all_coefs) <- c(paste0('coef', seq(1, b1)), paste0('coef', seq(1, b2), '_2'), paste0('coef', seq(1, b3), '_3'))
  
  only_t2_coefs1 <- matrix(rep(c(q_coef1, t2_coef1), n), nrow = n, byrow = T)
  colnames(only_t2_coefs1) <- paste0('t2_mu1_', seq(1, b2))
  
  only_t2_coefs2 <- matrix(rep(c(q_coef2, t2_coef2), n), nrow = n, byrow = T)
  colnames(only_t2_coefs2) <- paste0('t2_mu2_', seq(1, b2))
  
  out_df <- data.frame('seed' = seed, 'n' = n, 'p' =p , 'b1' = b1, 'b2' = b2, 'b3' = b3,  'corr' = corr,
                       'q_sd' = q_sd, 't1' = t1, 't2' = t2, 'lb' = lb, 'ub' = ub, x,
                       all_coefs, 'mean_coef1_2' = func_int2_coefs[,1],
                       only_t2_coefs1, only_t2_coefs2, 'mu_pt1' = mu_pt1, 'mu' = mu, 'y' = y)
  
  
  return(out_df)
}



gen_q <- function(out_df, func_df, func_coef_out, delta, cores){
  n <- nrow(out_df)
  t1 <- unique(out_df$t1)
  t2 <- unique(out_df$t2)
  int1 <- seq(1, t1 -1)
  int2 <- seq(t1, t2)
  int3 <- seq(t2+1, 300)
  
  lb <- unique(out_df$lb)
  ub <- unique(out_df$ub)
  sigma <- unique(out_df$q_sd)
  
  q_all <- mclapply(seq(1, n), function(j){
    mu1 <- out_df$mean_coef1_2[j]
    actual_coef <- out_df$coef1_2[j]
    f <- sapply(seq(lb, ub, by = 0.1), q_func, mu = mu1, sigma = sigma, delta = 0, lb = lb, ub = ub)
    q <- sapply(seq(lb, ub, by = 0.1), q_func, mu = mu1, sigma = sigma, delta = delta, lb = lb, ub = ub)
    
    f_val <- q_func(actual_coef, mu = mu1, sigma = sigma, delta = 0, lb = lb, ub = ub) /(sum(f))
    q_val <- q_func(actual_coef, mu = mu1, sigma = sigma, delta = delta, lb = lb, ub = ub)/(sum(q))
    q_mean <- sum(seq(lb, ub, by = 0.1)*(q/(sum(q))))
    q_random <- sample(seq(lb,ub, by = 0.1), size=1, replace=TRUE, prob=q/(sum(q)))
    
    return(c('q_mean' = q_mean,'q_random' = q_random, 'd_ratio' = q_val/f_val))
  }, mc.cores = cores)
  
  
  q_mean_coef <- sapply(q_all, function(j){
    return(j[1])
  })
  q_random_coef <- sapply(q_all, function(j){
    return(j[2])
  })
  q_d_ratio <- sapply(q_all, function(j){
    return(j[3])
  })
  
  new_coefs <- as.matrix(cbind(q_random_coef, out_df[paste0('coef', seq(2, unique(out_df$b2)), '_2')]))
  
  t2_coefs1 <- out_df %>%
    dplyr::select(paste0('t2_mu1_', seq(1, unique(out_df$b2)))) %>%
    unique() %>%
    unlist() %>%
    unname()
  t2_coefs2 <- out_df %>%
    dplyr::select(paste0('t2_mu2_', seq(1, unique(out_df$b2)))) %>%
    unique() %>%
    unlist() %>%
    unname()
  
  mu_q <-  ifelse(out_df$mu_pt1 < 2, out_df$mu_pt1 + new_coefs %*% t2_coefs1,  
                  out_df$mu_pt1 + new_coefs %*% t2_coefs2)
  
  func2 <- func_coef_out$func2
  gamma1 <- func2[, 1]
  funcs_int2_new <- t(sapply(seq(1, nrow(out_df)), function(j){
    og <- out_df$mean_coef1_2[j]
    q <- q_mean_coef[j]
    #error2 <- rnorm(1, 0, 5)
    error2 <- 0 #have to do this now that I am doing variations on truncated normal
    f <- q*gamma1- og*gamma1
    return(f)
  }))
  
  
  func_new_df <- data.frame(funcs_int2_new)
  colnames(func_new_df) <- c(int2)
  
  func_new_df <- func_new_df %>%
    mutate(id = seq(1, n)) %>%
    pivot_longer(cols = as.character(c(int2)), values_to = 'q_value') %>%
    right_join(func_df, by = c('id', 'name')) %>%
    mutate(mean_value = ifelse(is.na(q_value), NA, mean_value + q_value)) %>%
    dplyr::select(-q_value, - value, - obs) %>%
    mutate(delta = delta) %>%
    arrange(id, as.numeric(name))
  
  
  out_q_df <- data.frame('delta' = delta, 'q_mean_coef' = q_mean_coef, 'q_random_coef' = q_random_coef, 'q_d_ratio' = q_d_ratio, 'mu_q' = mu_q)
  
  return(list('out_q' = out_q_df, 'func_q_df' = func_new_df))
  
}


gen_q_no_func <- function(out_df, func_df, func_coef_out, delta, cores){
  n <- nrow(out_df)
  t1 <- unique(out_df$t1)
  t2 <- unique(out_df$t2)
  int1 <- seq(1, t1 -1)
  int2 <- seq(t1, t2)
  int3 <- seq(t2+1, 300)
  
  lb <- unique(out_df$lb)
  ub <- unique(out_df$ub)
  sigma <- unique(out_df$q_sd)
  
  q_all <- mclapply(seq(1, n), function(j){
    mu1 <- out_df$mean_coef1_2[j]
    actual_coef <- out_df$coef1_2[j]
    f <- sapply(seq(lb, ub, by = 0.1), q_func, mu = mu1, sigma = sigma, delta = 0, lb = lb, ub = ub)
    q <- sapply(seq(lb, ub, by = 0.1), q_func, mu = mu1, sigma = sigma, delta = delta, lb = lb, ub = ub)
    
    f_val <- q_func(actual_coef, mu = mu1, sigma = sigma, delta = 0, lb = lb, ub = ub) /(sum(f))
    q_val <- q_func(actual_coef, mu = mu1, sigma = sigma, delta = delta, lb = lb, ub = ub)/(sum(q))
    q_mean <- sum(seq(lb, ub, by = 0.1)*(q/(sum(q))))
    q_random <- sample(seq(lb,ub, by = 0.1), size=1, replace=TRUE, prob=q/(sum(q)))
    
    return(c('q_mean' = q_mean,'q_random' = q_random, 'd_ratio' = q_val/f_val))
  }, mc.cores = cores)
  
  
  q_mean_coef <- sapply(q_all, function(j){
    return(j[1])
  })
  q_random_coef <- sapply(q_all, function(j){
    return(j[2])
  })
  q_d_ratio <- sapply(q_all, function(j){
    return(j[3])
  })
  
  new_coefs <- as.matrix(cbind(q_random_coef, out_df[paste0('coef', seq(2, unique(out_df$b2)), '_2')]))
  
  t2_coefs1 <- unname(unlist(unique(out_df[, paste0('t2_mu1_', seq(1, unique(out_df$b2)))])))
  t2_coefs2 <- unname(unlist(unique(out_df[, paste0('t2_mu2_', seq(1, unique(out_df$b2)))])))
  
  
  mu_q <-  ifelse(out_df$mu_pt1 < 2, out_df$mu_pt1 + new_coefs %*% t2_coefs1,  
                  out_df$mu_pt1 + new_coefs %*% t2_coefs2)
  
  
  out_q_df <- data.frame('delta' = delta, 'q_mean_coef' = q_mean_coef, 'q_random_coef' = q_random_coef, 'q_d_ratio' = q_d_ratio, 'mu_q' = mu_q)
  
  return(out_q_df)
  
}


plot_q_dist <- function(delta_vec, mu, sigma, lb, ub){
  q_dists <- lapply(delta_vec, function(j){
    first_q <- sapply(seq(lb,ub, by = 0.1), q_func, mu = mu, sigma = sigma, delta = j, lb = lb, ub = ub)
    j_df <- data.frame('delta' = j, 'vals' = seq(lb,ub, by = 0.1), 'pdf' = first_q/(sum(first_q)))
    return(j_df)
  })
  delta_plot <- do.call(rbind, q_dists)
  ggplot(delta_plot, aes(x = vals, y = pdf, color = as.factor(delta))) +
    geom_line() +
    theme_bw() +
    scale_color_viridis(option = 'A', discrete = T, begin = 0.2, end = 0.9)
}
