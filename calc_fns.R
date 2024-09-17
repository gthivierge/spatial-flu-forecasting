# various little helper functions
# do calculations, etc
# excluding interval scores and weighted interval scores


make_indic <- function(x, tol){
  
  if(!is.na(x)){
    if(x < 0 & abs(x) > tol){
      y <- -1
    } else if(x > tol){
      y <- 1
    } else{
      y <- 0
    }
  } else{
    y <- 0
  }
  
  return(y)
}



calculate_lags <- function(df, lagvar, lags, lag_diff, self, drop_cols = TRUE, indic, tol){
  
  new_diff <- function(x, lag){
    x1 <- rep(NA, lag)
    x2 <- diff(x, lag = lag)
    
    y <- c(x1, x2)
    return(y)
  }
  
  if(self == TRUE){
    lag_diff <- FALSE # don't do differencing for state's lagged values, only do neighbors
  }
  
  if(lag_diff == TRUE){
    lags <- lags[which(lags > 0)]
    map_lag_diff <- lags %>% map(~purrr::partial(new_diff, lag = .x))
    if(drop_cols == TRUE){
      df2 <- df %>%
        mutate(across(.cols = {{lagvar}}, .fns = map_lag_diff, .names = "lag_{.col}_{lags}wk")) %>%
        rename_with(.cols = {{lagvar}}, ~paste0("lag_", .x, "_0wk")) %>%
        dplyr::select(contains("wk"))
    } else{
      df2 <- df %>%
        mutate(across(.cols = {{lagvar}}, .fns = map_lag_diff, .names = "lag_{.col}_{lags}wk"))
    }
  } else{ # don't do differencing
    map_lag_nodiff <- lags %>% purrr::map(~partial(dplyr::lag, n = .x))
    if(drop_cols == TRUE){
      df2 <- df %>%
        mutate(across({{lagvar}}, .fns = map_lag_nodiff, .names = "lag_{.col}_{lags}wk")) %>%
        dplyr::select(contains("wk"))
      
    } else{
      df2 <- df %>%
        mutate(across(.cols = {{lagvar}}, .fns = map_lag_nodiff, .names = "lag_{.col}_{lags}wk"))
    }
    # map_lag <- lags %>% map(~partial(dplyr::lag, n = .x))
  }
  
  # only do the indicators if lags are already true
  # and if self is true, lags = false
  # so no indicators on own state
  if(lag_diff == TRUE & indic == TRUE){
    df2 <- df2 %>%
      rowwise() %>%
      dplyr::select(!contains("0")) %>%
      mutate(across(starts_with("lag"), .fns = ~make_indic(.x, tol)))
  }
  
  return(df2)
}


logxplusc <- function(x,c){
  return(log(x+c))
}

# if computing multiple quantiles with quantile regression, sort them so they don't cross
sort_quantiles <- function(out_qr){
  out_qr_sorted <- out_qr %>%
    group_by(prediction_date) %>%
    arrange(prediction) %>%
    mutate(quantile = quantiles) %>%
    ungroup()
  return(out_qr_sorted)
}


mk_long <- function(full_df, model_type){
  if(class(full_df) == "data.frame"){
  
    if(model_type == "quantile" | model_type == "linear"){
      long_df <- full_df %>%
        pivot_longer(starts_with("lag"), names_to = "lag_col", values_to = "lag_value") %>%
        pivot_longer(c("mse","WIS"), names_to = "error_type", values_to = "error_value") %>%
        pivot_longer(starts_with("beta"), names_to = "beta_col", values_to = "beta_value") %>%
        # dplyr::select(-int) %>%
        as.data.frame()
      
    }  else if(model_type == "poisson"){
      long_df <- full_df %>%
        pivot_longer(starts_with("lag"), names_to = "lag_col", values_to = "lag_value") %>%
        pivot_longer(c("mse","se","se_df"), names_to = "error_type", values_to = "error_value") %>%
        pivot_longer(starts_with("beta"), names_to = "beta_col", values_to = "beta_value") %>%
        # dplyr::select(-int) %>%
        mutate(quantile = NA)   
    } else if(model_type == "lvcf"){
      long_df <- full_df %>%
        pivot_longer(c("mse", "mad"), names_to = "error_type", values_to = "error_value") %>%
        pivot_longer(starts_with("beta"), names_to = "beta_col", values_to = "beta_value") 
    }
    
    return(long_df)
  } else
    return(NA)
}


make_big_df <- function(flu_data, prim_reg, lags = c(0,1,2), 
                               model_type, inc_nbrs = TRUE, neg_pred,
                               pred_start_dt, in_season = TRUE,
                        quantiles, as_of, outcome, lag_diff, predictor,
                        indic, tol, natl, pooled, cores) {
  
  if(length(prim_reg)==1){
    out_lst <- fit_models(flu_data = flu_data, 
                        lags = lags,
                        prim_reg = prim_reg,
                        model_type = model_type,
                        inc_nbrs = inc_nbrs, 
                        neg_pred = neg_pred,
                        pred_start_dt = pred_start_dt, 
                        in_season = in_season,
                        quantiles = quantiles,
                        as_of = as_of,
                        outcome = outcome,
                        lag_diff = lag_diff,
                        predictor = predictor,
                        indic = indic, 
                        tol = tol,
                        natl = natl,
                        pooled = pooled)
    
    
    long_lst <- mk_long(out_lst, model_type = model_type)
    
    df_out <- as.data.frame(long_lst)
    
  } else{
    
    # prim_reg2 <- prim_reg[7]
    out_lst <- mclapply(X = prim_reg, FUN = fit_models, mc.cores = cores, mc.silent = TRUE, flu_data = flu_data, 
                        lags = lags,
                        model_type = model_type,
                        inc_nbrs = inc_nbrs, 
                        neg_pred = neg_pred,
                        pred_start_dt = pred_start_dt, 
                        in_season = in_season,
                        quantiles = quantiles,
                        as_of = as_of,
                        outcome = outcome,
                        lag_diff = lag_diff,
                        predictor = predictor,
                        indic = indic, 
                        tol = tol,
                        natl = natl,
                        pooled = pooled)
    
    names(out_lst) <- prim_reg
    long_lst <- lapply(out_lst, mk_long, model_type = model_type)
    
    
    df_lst <- rbindlist(long_lst[!is.na(long_lst)])
    df_out <- as.data.frame(df_lst)
  }
  
  return(df_out)
}


make_small_df <- function(big_df, model_type){
  if(length(names(big_df)) > 0){
    
  if(model_type == "quantile"){
    small_df <- kit::funique(big_df[,which(names(big_df) %in% c("target_date", "prediction", "region", "epiweek", "nbrs", "model_type", "error_type", "error_value", "quantile", "WIS", "mad"))]) %>%
      filter(!is.na(region)) %>%
      filter(!is.na(prediction)) %>%
      pivot_wider(names_from = error_type, values_from = error_value) %>%
      mutate(se = NA, se_df = NA) %>%
      dplyr::select(target_date, epiweek, region, prediction, quantile, WIS, mse, mad, se, se_df, nbrs, model_type)
  } else if(model_type == "linear"){
    small_df <- kit::funique(big_df[, which(names(big_df) %in% c("target_date", "prediction", "region", "epiweek", "nbrs", "model_type", "error_type", "error_value", "WIS", "mad"))]) %>%
      filter(!is.na(region)) %>%
      mutate(se = NA, quantile = NA, se_df = NA) %>%
      pivot_wider(names_from = error_type, values_from = error_value) %>%
      dplyr::select(target_date, epiweek, region, prediction, quantile, WIS, mse, mad, se, se_df, nbrs, model_type)
   } else if(model_type == "poisson"){
    small_df <- kit::funique(big_df[,which(names(big_df) %in% c("target_date", "prediction", "region", "epiweek", "nbrs", "model_type", "error_type", "error_value", "quantile", "WIS", "mad"))]) %>%
      filter(!is.na(region)) %>%
      filter(!is.na(prediction)) %>%
      pivot_wider(names_from = error_type, values_from = error_value) %>%
      dplyr::select(target_date, epiweek, region, prediction, quantile, WIS, mse, mad, se, se_df, nbrs, model_type)
  } else if(model_type == "lvcf"){
    small_df <- kit::funique(big_df[,which(names(big_df) %in% c("target_date", "prediction", "region", "epiweek", "nbrs", "model_type", "error_type", "error_value", "quantile", "WIS", "mad"))]) %>%
      filter(!is.na(region)) %>%
      filter(!is.na(prediction)) %>%
      pivot_wider(names_from = error_type, values_from = error_value) %>%
      mutate(se = NA, se_df = NA) %>%
      dplyr::select(target_date, epiweek, region, prediction, quantile, WIS, mse, mad, se, se_df, nbrs, model_type)
  }
  
  small_df$target_date <- as.Date(small_df$target_date)
  return(small_df)
  }
}

conformal_wis <- function(out_matrix, quantiles, df){
  

  out <- as.data.frame(out_matrix)
  names(out) <- paste0("quantile_", tau)
  # out$quantile_0.5 <- df$prediction
  out$target_date <- df$target_date
 
  out$Y <- df$ili

  out <- na.omit(out)
  
  # quantiles <- pull(quantiles)
  
  q_lower <- quantiles[which(quantiles < 0.5)]
  q_med <- 0.5
  alpha_k <- 2 * q_lower  
  width <- 1 - 2 * q_lower
  
  IS_all <- as.data.frame(cbind(out$Y, out$quantile_0.5))
  names(IS_all) <- c("ili", "median")
  
  cover
  
  for(k in 1:length(q_lower)){
    q_col_lower <- out[,which(names(out) == paste0("quantile_", q_lower[k]))] 
    q_col_upper <- out[,which(names(out) == paste0("quantile_", 1-q_lower[k]))]
    # IS_nm <- paste0("IS_", width)
    IS <- scoringutils::interval_score(out$Y, q_col_lower, q_col_upper, width[k] * 100, weigh = FALSE)
    # IS_by_hand(.05, q_col_upper, q_col_lower, Y)
    IS_all <- cbind(IS_all, IS)
    names(IS_all)[length(IS_all)] <- paste0("IS_", width[k])
  }
  
  # names(IS_all)[1:2] <- c(outcome, "median")
  
  WIS <- weighted_IS(q_lower, IS_all, "ili")
  
  newdf <- bind_cols("target_date" = out$target_date, "WIS.conf" = WIS)
  # names(newdf) <- c("target_date", "WIS.conf")

  df2 <- left_join(df, newdf, by = "target_date")
    
  return(df2)
}
