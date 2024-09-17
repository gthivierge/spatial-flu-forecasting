# quantile regression

options(scipen = 999)

# library(glmnet, quietly = TRUE)
# library(hqreg, quietly = TRUE)


get_lvcf_models <- function(cutoff_dt, ar_data_full_r, model_df_i, in_season, as_of, outcome, ar_lst) {

  
  if(is.null(ar_lst)){
    ar_lst <- arrange_data(ar_data_full_r, nbr_lst = list(), cutoff_dt, lags = c(), in_season, as_of, outcome, lag_diff = FALSE, predictor = outcome, "lvcf", indic = FALSE)
  } else{
    ar_lst <- ar_lst
  }
  
  ar_X <- as.data.frame(ar_lst[3]) %>%
    dplyr::select(-c(issue, {outcome}))
  ar_Y <- as.data.frame(ar_lst[3]) %>%
    dplyr::select(c(date, {outcome}))
  
  X_i <- ar_X[which(ar_X$date <= cutoff_dt),] %>%
    dplyr::select(-date)
  
  Y_i <- ar_Y[which(ar_Y$date <= cutoff_dt),] %>%
    dplyr::select(-date)
  
  # shift the Y values so we are predicting 2 weeks ahead
  X_i <- X_i[1:(nrow(X_i)-2),]
  y_i <- as.data.frame(Y_i[3:nrow(Y_i),])
  names(y_i) <- {outcome}
  
  Yhat <- ar_Y[which(ar_Y$date <= cutoff_dt),] %>%
    dplyr::select({outcome}) %>%
    lag(n=2) # the outcome week is t+2, and the "last value" will be just the current value
  
  yhat <- c(Yhat[3:nrow(Yhat),])
  
  full_xy <- bind_cols(X_i, y_i) %>%
    as.data.frame() %>%
    bind_cols(prediction = yhat)
  
  model_df_i$beta_0 <- NA_real_
  model_df_i$prediction_date <- cutoff_dt %m-% weeks(2)
  model_df_i$target_date <- cutoff_dt
  model_df_i$mse <- sum((yhat - y_i)^2)/nrow(y_i)
  model_df_i$mad <- sum(abs(yhat - y_i))/nrow(y_i)
  model_df_i$prediction <- full_xy %>%
    tail(n=1) %>%
    dplyr::select(prediction) %>% pull()
   
  
  return(model_df_i)
  
}


full_qreg <- function(data, tau, neg_pred, outcome, ar_lst){
  
  newrow <- data[nrow(data),]
  data2 <- data[-nrow(data),]



    
    N_row <- nrow(data)
    
    newcols <- colnames(data2)[-which(colnames(data2) %in% c({outcome}))]
    
    newform <- reformulate(termlabels = c(sprintf("`%s`", newcols)), response = {outcome})
   
    z <- rq(formula = newform, tau = tau, data = data2) 
    
    Y <- data2[, which(names(data2) == {outcome})]
    
    if(neg_pred == TRUE){
      mse <- mean(z$residuals^2)    
      mad <- mean(abs(z$residuals))
    } else{
      fitted_2 <- pmax(z$fitted.values, 0)
      mse <- sum((Y - fitted_2)^2)/length(fitted_2)
      mad <- sum(abs(Y - fitted_2))/length(fitted_2)
    }
    
    pred_0 <- predict(z, newdata = newrow)
    pred <- as.numeric(if_else(neg_pred == FALSE, pmax(pred_0, 0), pred_0))
    
    m <- list("coeff" = z$coefficients, "MSE" = mse, "MAD" = mad, "pred" = pred)
  
  return(m)
}


get_qx_models <- function(cutoff_dt, ar_data_full_r, nbr_lst, quantiles, model_df_i, 
                          neg_pred, outcome, as_of, lags, in_season, lag_diff, ar_lst = NULL, predictor,
                           prim_reg){
  
  if(is.null(ar_lst)){
    ar_lst <- arrange_data(ar_data_full_r, nbr_lst, cutoff_dt, lags, in_season, as_of, outcome, lag_diff, predictor, "quantile", natl)
  } else{
    ar_lst <- ar_lst
  }
  
  # for pooling, reorder to predict on correct state
  ar_ids <- as.data.frame(ar_lst[3]) %>%
    mutate(prim = if_else(region == prim_reg, 1, 0)) %>%
    dplyr::select(prim)
  
  ar_X <- as.data.frame(ar_lst[1]) %>%
    dplyr::select(-issue) %>%
    bind_cols(ar_ids)
  ar_Y <- as.data.frame(ar_lst[2]) %>%
    bind_cols(ar_ids)
  
  # subset data using dates of interest
  
  X_i <- ar_X[which(ar_X$date <= cutoff_dt),] %>%
    arrange(prim, date) %>%
    filter(prim == 1 | date < cutoff_dt) %>%
    dplyr::select(-c(date, prim))
  
  Y_i <- ar_Y[which(ar_Y$date <= cutoff_dt),] %>%
    arrange(prim, date) %>%
    filter(prim == 1 | date < cutoff_dt) %>%
    dplyr::select(-c(date, prim))
  
  # shift the Y values so we are predicting 2 weeks ahead
  X_i <- X_i[1:(nrow(X_i)-2),]
  y_i <- as.data.frame(Y_i[3:nrow(Y_i),])
  names(y_i) <- {outcome}
  
  full_xy <- bind_cols(X_i, y_i) %>%
    as.data.frame() %>%
    dplyr::select(-int)
  
  
  # quantiles <- as.numeric(quantiles$quantile)
  
  reg <- lapply(quantiles, 
                full_qreg, 
                data = full_xy, 
                neg_pred = neg_pred,
                outcome = outcome)
  
  reg_2 <- as.data.frame(do.call(rbind, reg))
  
  coeff_df <- as.data.frame(do.call(rbind, reg_2$coeff))

  mse_df <- as.data.frame(do.call(rbind, reg_2$MSE))
  names(mse_df)<-"mse"
  mad_df <- as.data.frame(do.call(rbind, reg_2$MAD))
  names(mad_df) <- "mad"
  # 
  model_df_i$quantile <- quantiles
    model_df_i$prediction_date <- rep(cutoff_dt  %m-% weeks(2), length(quantiles))
  model_df_i$target_date <- rep(cutoff_dt, length(quantiles))

  model_df_i$mse <- mse_df$mse
  model_df_i$mad <- mad_df$mad
  model_df_i[,6:ncol(model_df_i)] <- coeff_df
  model_df_i$prediction <- as.numeric(unlist(reg_2$pred))
  # out <- list("model_df_i" = model_df_i, "predictions" = fitted_vals)
  return(model_df_i)
}


######### linear regression ##########

full_linreg <- function(data, neg_pred, outcome, quantiles){
  
  newrow <- data[nrow(data),]
  data2 <- data[-nrow(data),]
  

        
    # fm <- paste({outcome}, "~ .", sep = " ")
    
    # find self vs neighbor lags, and remove the ili column
    
   st_cols <- str_extract(names(data2), "(?<=lag_).+(?=_[0-2]wk)")[-ncol(data2)]
    
   # only want to constrain the neighbor lags
   # coefficients must be positive

   
   A <- data2 %>%
     dplyr::select(-c({outcome})) %>% 
     mutate(int = 1) %>%
     relocate(int) %>%
     as.matrix()
   
   b <- data2[, which(names(data2) == {outcome})]

   # intercept and neighbor lags positive
   
   
   # if(constrain == TRUE){
   #   
   #   # bootstrap <- TRUE
   #   
   #   lower <- c(0, if_else(is.na(st_cols), -10000000000, 0))
   #   lin_out <- lsei(A, b, lower = lower, type = 2)
   # } else{
   #   if(bootstrap == FALSE){
   #     # fm <- paste({outcome}, "~ .", sep = " ")
   #     # z <- lm(formula = fm, data = data2)
   #     # sum_z <- summary(z)
   #   } else{
   #     lin_out <- lsei(A, b, type = 2)
   #   }
   # }
   

  lin_out <- lsei(A, b, type = 2)

    if(neg_pred == TRUE){
     #  mse <- mean(z$residuals^2)
     #  rse <- sum_z$sigma
     #  mad <- mean(abs(z$residuals))
     #  pred <- as.numeric(predict(z, newdata = newrow))
     # 
     # betas <- z$coefficients
      
      betas <- lin_out$X
      
      train_pred <- A %*% betas
      
      train_resid <- train_pred - b
      
      test_x <- newrow %>%
        dplyr::select(-c({outcome})) %>% 
        mutate(int = 1) %>%
        relocate(int) %>%
        as.numeric()
      
      mse <- mean(train_resid^2)
      mad <- mean(abs(train_resid))
      pred <- sum(test_x * betas)
      
    } else{ # compute statistics with positive predictions only

      betas <- lin_out$X
      
      train_pred <- pmax(A %*% betas, 0)
      
      train_resid <- train_pred - b
      
      test_x <- newrow %>%
        dplyr::select(-c({outcome})) %>% 
        mutate(int = 1) %>%
        relocate(int) %>%
        as.numeric()
      
      mse <- mean(train_resid^2)
      mad <- mean(abs(train_resid))
      
      k <- length(betas) - 1
      sse <- sum(train_resid^2)
      n <- length(train_resid)
      rse <- sqrt(sse/(n-(1+k)))
      pred_int <- sqrt(rse^2 * (1 + t(test_x) %*% solve(t(A) %*% A) %*% test_x))


      se_df <- n - (1 + k)
    
    se <- pred_int
    pred <- sum(test_x * betas)
  }
  
  m <- list("coeff" = betas, 
            "MSE" = mse, 
            "se" = se, 
            "se_df" = se_df, 
            "prediction" = pred,
            "MAD" = mad)  
  
  
  # m <- list("coeff" = z$coefficients, 
  #           "MSE" = mse, 
  #           "se" = se, 
  #           "se_df" = se_df, 
  #           "prediction" = pred,
  #           "MAD" = mad)
  return(m)
}

get_lin_models <- function(cutoff_dt, ar_data_full_r, nbr_lst, quantiles, model_df_i, neg_pred, 
                            outcome, as_of, lags, in_season, lag_diff, ar_lst = NULL, predictor,
                          conformal, prim_reg, indic){

  if(is.null(ar_lst)){
    ar_lst <- arrange_data(ar_data_full_r, nbr_lst, cutoff_dt, lags, in_season, as_of, outcome, lag_diff, predictor, "linear", indic)
  } else{
    ar_lst <- ar_lst
  }
  
  
  # for pooling, reorder to predict on correct state
  ar_ids <- as.data.frame(ar_lst[3]) %>%
    mutate(prim = if_else(region == prim_reg, 1, 0)) %>%
    dplyr::select(prim)
  
  ar_X <- as.data.frame(ar_lst[1]) %>%
    dplyr::select(-issue) %>%
    bind_cols(ar_ids)
  ar_Y <- as.data.frame(ar_lst[2]) %>%
    bind_cols(ar_ids)
  
  # subset data using dates of interest
  
  X_i <- ar_X[which(ar_X$date <= cutoff_dt),] %>%
    arrange(prim, date) %>%
    filter(prim == 1 | date < cutoff_dt) %>%
    dplyr::select(-c(date, prim))
  
  Y_i <- ar_Y[which(ar_Y$date <= cutoff_dt),] %>%
    arrange(prim, date) %>%
    filter(prim == 1 | date < cutoff_dt) %>%
    dplyr::select(-c(date, prim))
  
  # shift the Y values so we are predicting 2 weeks ahead
  X_i <- X_i[1:(nrow(X_i)-2),]
  y_i <- as.data.frame(Y_i[3:nrow(Y_i),])
  names(y_i) <- {outcome}
  
  full_xy <- bind_cols(X_i, y_i) %>%
    as.data.frame() %>%
    dplyr::select(-int)
  
  # quantiles <- as.numeric(quantiles$quantile)
  
  reg_2 <- full_linreg(full_xy, 
                       neg_pred = neg_pred,
                       outcome = outcome,
                       quantiles = quantiles)
  

# df <- as.data.frame(reg_2)

    # df <- as.data.frame(reg_2)
    coeff_df <- reg_2$coeff
    model_df_i$prediction_date <- cutoff_dt %m-% weeks(2)
    model_df_i$target_date <- cutoff_dt
    model_df_i$mse <- reg_2$MSE
    model_df_i$mad <- reg_2$MAD
    model_df_i$se <- reg_2$se
    model_df_i$se_df <- reg_2$se_df
    
    
    if(neg_pred == TRUE){
      model_df_i$prediction <- reg_2$prediction    
    } else{
      model_df_i$prediction <- pmax(reg_2$prediction, 0)
    }
    
    last_beta_col <- ncol(model_df_i) - 3
    model_df_i[,6:last_beta_col] <- coeff_df
    
  return(model_df_i)
}

######## poisson regression ##########
full_poisson <- function(data){

  # log of 1/2 the smallest number

  # add to all the numbers (pseudocount)
  
  # data[data == 0] <- .Machine$double.eps
  newrow <- data[nrow(data),]
  data2 <- data[-nrow(data),]
      
 newcols <- colnames(data2)[-which(colnames(data2) %in% c("num_patients", "num_ili"))]
      
 newform <- reformulate(termlabels = c(sprintf("`%s`", newcols), "offset(num_patients)"), response = "num_ili")
      
      
    z <- glm(newform, family = "poisson",
              data = data2)
    
    mse <- mean(z$residuals^2)
    mad <- mean(abs(z$residuals))
    

    
    pred <- predict.glm(z, 
                    family = "poisson",
                    newdata = newrow, 
                    type = "response", 
                    se.fit = TRUE)
    
    
    newrow2 <- newrow %>%
      dplyr::select(-c(num_ili, num_patients)) %>%
      mutate(int = 1) %>%
      relocate(int)
    
    sf<-summary(z)
    Xo <- as.matrix(newrow2)
    # Xo<-as.matrix(newrow2)
    lam0<-pred$fit
    #note: trace(V*xx') is equivalent to x'Vx
    term1<-apply(Xo,1,function(x){x%*% crossprod(sf$cov.scaled, x)})
    #term1<-diag(Xo%*% sf$cov.scaled %*% t(Xo)) #more expensive way to do same thing
    Vo<- 1+(1)^2*lam0*term1
    sqlamV<-sqrt(lam0*Vo)
    
    m <- list("coeff" = z$coefficients, 
              "MSE" = mse, 
              "se" = sqlamV/(exp(newrow$num_patients)) * 100,
              "se_df" = round(as.numeric(pred$fit)), 
              "prediction" = as.numeric(pred$fit)/(exp(newrow$num_patients)) * 100,
              "MAD" = mad)

  return(m)
}

get_poisson_models <- function(cutoff_dt, ar_data_full_r, nbr_lst, quantiles, model_df_i, 
                               neg_pred, outcome, as_of, lags, in_season, lag_diff, ar_lst = NULL, predictor, prim_reg){
 
  if(is.null(ar_lst)){
    ar_lst <- arrange_data(ar_data_full_r, nbr_lst, cutoff_dt, lags, in_season, as_of, outcome, lag_diff, predictor, "poisson")
  } else{
    ar_lst <- ar_lst
  }
  c1 <- 0.002207895 # pseudocount - half of smallest ILI 
  
  ar_ids <- as.data.frame(ar_lst[3]) %>%
    mutate(prim = if_else(region == prim_reg, 1, 0)) %>%
    dplyr::select(prim)
  
  
  if(pooled == TRUE){
    ar_X <- as.data.frame(ar_lst[1]) %>%
      left_join((ar_data_full_r %>%
                   dplyr::select(c(date, issue, region, num_patients, num_providers))), by = c("date", "issue", "region")) %>%
      dplyr::select(-c(issue, int)) %>%
      rowwise() %>%
      mutate(across(.cols = c(num_patients, num_providers), .fns = ~logxplusc(.x, c1))) %>%
      ungroup() %>%
      dplyr::select(-region) %>%
      bind_cols(ar_ids)
  } else{
    ar_X <- as.data.frame(ar_lst[1]) %>%
      left_join((ar_data_full_r %>%
                   dplyr::select(c(date, issue, num_patients, num_providers))), by = c("date", "issue")) %>%
      dplyr::select(-c(issue, int)) %>%
      rowwise() %>%
      mutate(across(.cols = c(num_patients, num_providers), .fns = ~logxplusc(.x, c1))) %>%
      ungroup() %>%
      dplyr::select(-region) %>%
      bind_cols(ar_ids)
  }
  
  # for pooling, reorder to predict on correct state

  ar_Y <- as.data.frame(ar_lst[2]) %>%
    bind_cols(ar_ids)
  
  # subset data using dates of interest
  
  X_i <- ar_X[which(ar_X$date <= cutoff_dt),] %>%
    arrange(prim, date) %>%
    filter(prim == 1 | date < cutoff_dt) %>%
    dplyr::select(-c(date, prim))
  
  Y_i <- ar_Y[which(ar_Y$date <= cutoff_dt),] %>%
    arrange(prim, date) %>%
    filter(prim == 1 | date < cutoff_dt) %>%
    dplyr::select(-c(date, prim))
  

  X_i <- X_i[1:(nrow(X_i)-2),]
  y_i <- as.data.frame(Y_i[3:nrow(Y_i),])

  names(y_i) <- "num_ili"
  
  full_xy <- bind_cols(X_i, y_i) %>%
    as.data.frame()
  
  reg_2 <- full_poisson(full_xy)
  
  coeff_df <- reg_2$coeff
  
  
  model_df_i$prediction_date <- cutoff_dt %m-% weeks(2)
  model_df_i$target_date <- cutoff_dt
  model_df_i$mse <- reg_2$MSE
  model_df_i$se <- reg_2$se
  model_df_i$se_df <- reg_2$se_df
  model_df_i$mad <- reg_2$MAD
  model_df_i$prediction <- reg_2$prediction  
  
  last_beta_col <- ncol(model_df_i) - 2
  model_df_i[,6:last_beta_col] <- coeff_df
  # out <- list("model_df_i" = model_df_i, "predictions" = fitted_vals)
  return(model_df_i)
}
