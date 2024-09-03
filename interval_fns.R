# functions related to interval scores and helpers to add WIS to data frames

## interval score
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008618


weighted_IS <- function(quantiles, IS_all, outcome) {
  q_lower <- quantiles[which(quantiles < 0.5)]
  q_med <- 0.5
  alpha_k <- 2 * q_lower  
  width <- 1 - 2 * q_lower
  
  
  WIS <- rep(NA_real_, nrow(IS_all))
  K <- length(q_lower)
  w_0 <- 0.5
  Y <- IS_all[,which(names(IS_all) == outcome)]
  med <- IS_all[,which(names(IS_all) == "median")]
  
    w_k <- alpha_k/2
    wt_IS <- IS_all %>%
      dplyr::select(starts_with("IS_")) %>%
      adply(., 1, function(x) sum(x * w_k)) %>%
      pull(V1)
    
    WIS <- (1/(K + 0.5)) * (w_0 * abs(Y - med) + wt_IS)
  
  
  return(WIS)
}


add_wis <- function(out_reg, quantiles, outcome, model_type){
  
  # quantiles <- unique(out_qr_sorted$quantile)
  # quantiles <- c(0.025, 0.05, 0.2, 0.25, 0.5, 0.75, 0.8, 0.95, 0.975)
  
  if(model_type == "quantile"){
    reg_temp <- out_reg %>% 
      dplyr::select(target_date, {outcome}, quantile, prediction) %>%
      pivot_wider(names_from = "quantile", values_from = "prediction", names_prefix = "quantile_") 
  } else if(model_type == "linear"){
    quantiles <- unique(out_reg$quantile)
    
    reg_temp <- out_reg %>% 
      dplyr::select(target_date, {outcome}, quantile, prediction, se) %>%
      mutate(bound = if_else(quantile == 0.5, prediction, prediction + se*qnorm(quantile))) %>%
      pivot_wider(names_from = "quantile", values_from = "bound", names_prefix = "quantile_") 
    
  } else if(model_type == "poisson"){
    quantiles <- unique(out_reg$quantile)

    reg_temp <- out_reg %>% 
      mutate(num_ili_pred = prediction/100 * num_patients) %>%
      dplyr::select(target_date, ili, quantile, prediction, num_ili_pred, num_patients, num_ili, se) %>%
      mutate(num_se = se/100 * num_patients) %>% # put back on count scale
      mutate(bound = if_else(quantile == 0.5, num_ili_pred, num_ili_pred + qnorm(quantile) * num_se)) %>%
      mutate(bound = bound / num_patients * 100) %>% # put back on ILI scale
      pivot_wider(names_from = "quantile", values_from = "bound", names_prefix = "quantile_") 
  }
  
  if(model_type == "poisson"){
    outcome <- "ili"
  }
  
  Y <- reg_temp[,which(names(reg_temp) == outcome)] %>% pull()
  
  q_lower <- quantiles[which(quantiles < 0.5)]
  q_med <- 0.5
  alpha_k <- 2 * q_lower  
  width <- 1 - 2 * q_lower
  
  IS_all <- as.data.frame(cbind(Y, reg_temp$quantile_0.5))
  names(IS_all) <- c(outcome, "median")
  
  for(k in 1:length(q_lower)){
    q_col_lower <- reg_temp[,which(names(reg_temp) == paste0("quantile_", q_lower[k]))] %>% pull()
    q_col_upper <- reg_temp[,which(names(reg_temp) == paste0("quantile_", 1-q_lower[k]))] %>% pull()
    # IS_nm <- paste0("IS_", width)
    IS <- scoringutils::interval_score(Y, q_col_lower, q_col_upper, width[k] * 100, weigh = FALSE)
    # IS_by_hand(.05, q_col_upper, q_col_lower, Y)
    IS_all <- cbind(IS_all, IS)
    names(IS_all)[length(IS_all)] <- paste0("IS_", width[k])
  }
  
  # names(IS_all)[1:2] <- c(outcome, "median")

    WIS <- weighted_IS(q_lower, IS_all, outcome)
  
  
  reg_temp$WIS <- WIS
  
  
  reg_temp_2 <- reg_temp %>%
    dplyr::select(c(target_date, WIS))
  
  if(model_type == "quantile"){
    out_reg_WIS <- out_reg %>%
      left_join(reg_temp_2, by = "target_date") %>%
      arrange(target_date, quantile) %>%
      as.data.frame()
  } else { # same format for poisson and quantile
    out_reg_WIS <- out_reg %>%
      left_join(reg_temp_2, by =c("target_date")) %>%
      arrange(target_date, quantile) %>%
      dplyr::select(-quantile) %>%
      unique() %>%
      as.data.frame()
  }
  
  return(out_reg_WIS)
}


add_wis_linear <- function(out_lin, outcome){
  
  quantiles <- unique(out_lin$quantile)
  
  lin_temp <- out_lin %>% 
    dplyr::select(target_date, {outcome}, quantile, prediction, se) %>%
    mutate(bound = if_else(quantile == 0.5, prediction, prediction + se*qnorm(quantile))) %>%
    pivot_wider(names_from = "quantile", values_from = "bound", names_prefix = "quantile_") 
  
  q_lower <- quantiles[which(quantiles < 0.5)]
  q_med <- 0.5
  alpha_k <- 2 * q_lower  
  width <- 1 - 2 * q_lower
  
  Y <- lin_temp[,which(names(lin_temp) == outcome)]
  
  IS_all <- cbind(Y, lin_temp$quantile_0.5)
  
  names(IS_all) <- c({outcome}, "median")
  
  
  for(k in 1:length(q_lower)){
    q_col_lower <- lin_temp[,which(names(lin_temp) == paste0("quantile_", q_lower[k]))] %>% pull()
    q_col_upper <- lin_temp[,which(names(lin_temp) == paste0("quantile_", 1-q_lower[k]))] %>% pull()
    IS <- scoringutils::interval_score(Y, q_col_lower, q_col_upper, width[k] * 100, weigh = FALSE)
    IS_all <- cbind(IS_all, IS)
    names(IS_all)[length(IS_all)] <- paste0("IS_", width[k])
  }
  
  names(IS_all)[1:2] <- c({outcome}, "median")
  
  WIS <- weighted_IS(q_lower, IS_all, outcome)
  
  lin_temp$WIS <- WIS
  
  lin_temp_2 <- lin_temp %>%
    dplyr::select(c(target_date, WIS))
  
  out_lin_WIS <- out_lin %>%
    left_join(lin_temp_2, by =c("target_date")) %>%
    arrange(target_date, quantile) %>%
    dplyr::select(-quantile) %>%
    unique() %>%
    as.data.frame()
  
  return(out_lin_WIS)
}


add_wis_poisson <- function(out_poi){
  
  quantiles <- unique(out_poi$quantile)
  
  poi_temp <- out_poi %>% 
    mutate(num_ili_pred = prediction * num_patients/100) %>%
    dplyr::select(target_date, ili, quantile, prediction, num_ili_pred, num_patients, se) %>%
    mutate(bound = if_else(quantile == 0.5, num_ili_pred, num_ili_pred + qnorm(quantile) * se)) %>%
    mutate(bound = bound / num_patients * 100) %>%
    pivot_wider(names_from = "quantile", values_from = "bound", names_prefix = "quantile_") 
  
  q_lower <- quantiles[which(quantiles < 0.5)]
  q_med <- 0.5
  alpha_k <- 2 * q_lower  
  width <- 1 - 2 * q_lower
  
  IS_all <- cbind(poi_temp$ili, poi_temp$quantile_0.5)
  names(IS_all) <- c("num_ili", "median")
  
  
  WIS <- weighted_IS(q_lower, IS_all, "num_ili")
  
  poi_temp$WIS <- WIS
  
  poi_temp_2 <- poi_temp %>%
    dplyr::select(c(target_date, WIS))
  
  out_poi_WIS <- out_poi %>%
    left_join(poi_temp_2, by =c("target_date")) %>%
    arrange(target_date, quantile) %>%
    dplyr::select(-quantile) %>%
    unique() %>%
    as.data.frame()
  
  return(out_poi_WIS)
}


# 
# 
# library(epipredict)
# 
# wis <- function(x, actual,...) {
#   UseMethod("wis")
# }
# 
# wis.distribution <- function(x, actual, ...) {
#   rlang::check_dots_empty()
#   epipredict:::map2_dbl(vctrs::vec_data(x), actual, wis_dist_quantile)
# }
# 
# wis_dist_quantile <- function(x, actual) {
#   q <- vctrs::field(x, "q")
#   if (all(is.na(q))) return(NA)
#   if (is.na(actual)) return(NA)
#   tau <- vctrs::field(x, "tau")
#   2 * mean(pmax(
#     tau * (actual - q),
#     (1 - tau) * (q - actual), na.rm = TRUE
#   ))
# }
# 
