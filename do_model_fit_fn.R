# function to manipulate data and call the model of choice (linear, quantile, poisson)

flu_data <- read.csv("flu_data_all.csv", header = TRUE)

# filter to in-season data
mk_in_season <- function(ar_data_full_r, nbr_lags){
  
  if(!is.null(nbr_lags)){
    nbr_lags <- nbr_lags %>%
      filter(WEEK <= 18 | WEEK >= 40)
  }
  
  
  ar_data_full_r <- ar_data_full_r %>%
    filter(WEEK <= 18 | WEEK >= 40)
  
  min_yr <- min(ar_data_full_r$YEAR)
  max_yr <- max(ar_data_full_r$YEAR)
  
  yrs <- seq(min_yr, max_yr,  by= 1)
  season_name <- c("2010-11",
                   "2011-12",
                   "2012-13",
                   "2013-14",
                   "2014-15",
                   "2015-16",
                   "2016-17",
                   "2017-18",
                   "2018-19")
  
  season_n <- seq(1,9)
  
  tbl <- as.data.frame(cbind(season_n, season_name)) %>%
    mutate(season_n = as.numeric(season_n))
  
  
  ar_data_full_r <- ar_data_full_r %>%
    mutate(season = YEAR - min_yr + if_else(WEEK <= 18, 0, 1)) %>%
    left_join(tbl, join_by(season == season_n))
  
  return(list(ar_data_full_r, nbr_lags))
}

# get list of neighbors of target state according to model variant
get_nbrs <- function(prim_reg, inc_nbrs, natl){
  st_mat <- read.csv( "states_mat_adj.csv", header=TRUE, row.names=1)
    if(inc_nbrs == TRUE){
      nbr_r <- st_mat[which(rownames(st_mat)==prim_reg),]
      nbr_lst <- colnames(nbr_r)[which(nbr_r > 0)]
      nbr_lst <- unique(nbr_lst)
    } else{
      nbr_lst <- list()
    }
  
  if(natl == TRUE){
    nbr_lst <- append(nbr_lst, "US")
  }
  
  return(nbr_lst)
}

# binarize the lagged week to week change in ILI in the neighbor states
do_nbr_lags <- function(nbr_lst, predictor, lag_diff, model_type, as_of, cutoff_dt, indic, tol){
  st_mat <- read.csv("states_mat_adj.csv", header=TRUE, row.names=1)
  
  ar_data_nbr <- flu_data %>%
    mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
    filter(region %in% nbr_lst) %>%
    dplyr::select(c(region, date, epiweek, YEAR, WEEK, {predictor}, issue)) %>%
    mutate(issue = as.Date(issue, c("%m/%d/%y")))
  
  if(as_of == FALSE){ # filter down the data to just most recent issue
    ar_nbr_asof <- ar_data_nbr %>%
      filter(date <= cutoff_dt) %>%
      group_by(region, epiweek) %>%
      slice_max(issue) %>%
      ungroup()
  } else{
    ar_nbr_asof <- ar_data_nbr %>%
      filter(date <= cutoff_dt) %>%
      filter(issue <= cutoff_dt) %>%
      group_by(region, date) %>%
      slice_max(issue) %>%
      ungroup()
  }  
  
  # take log value for poisson
  if(model_type == "poisson"){
    c1 <- min(flu_data$ili[flu_data$ili > 0])/2
    ar_nbr_asof <- ar_nbr_asof %>%
      rowwise() %>%
      mutate(across(.cols = c({predictor}), .fns = ~logxplusc(.x, c1))) %>%
      ungroup()
  }
  
  ar_x_nbr <- ar_nbr_asof %>%
    arrange(date) %>%
    dplyr::select(region, date, {predictor}, WEEK) %>%
    pivot_wider(names_from = region, values_from = {predictor})
  
  
  # calculate lagged values for neighbors
  nbr_lags_0  <- lapply(nbr_lst, calculate_lags, df = ar_x_nbr, lags = lags, lag_diff = lag_diff, self = FALSE, indic = indic, tol = tol)
  nbr_lags <- do.call(cbind, nbr_lags_0)
  nbr_lags$WEEK <- ar_x_nbr$WEEK
  nbr_lags$date <- ar_x_nbr$date
  
  return(nbr_lags)
}

# filter to final data release for each week
filter_issues <- function(ar_data_full_r, cutoff_dt, model_type){
  if(as_of == FALSE){ # filter down the data to just most recent issue
    ar_data_asof <- ar_data_full_r %>%
      filter(date <= cutoff_dt) %>%
      group_by(region, epiweek) %>%
      slice_max(issue) %>%
      ungroup()
  } else{
    ar_data_asof <- ar_data_full_r %>%
      filter(date <= cutoff_dt) %>%
      filter(issue <= cutoff_dt) %>%
      group_by(date) %>%
      slice_max(issue) %>%
      ungroup()
  }
  return(ar_data_asof)
}

# helper function to prepare data
arrange_data <- function(ar_data_full_r, nbr_lst, cutoff_dt, lags, 
                         in_season,
                         as_of, outcome, lag_diff, predictor, model_type, indic, natl){
  
  if(as_of == FALSE){ # filter down the data to just most recent issue
    ar_data_asof <- ar_data_full_r %>%
      filter(date <= cutoff_dt) %>%
      group_by(region, epiweek) %>%
      slice_max(issue) %>%
      ungroup() %>%
      arrange(date)
  } else{
    ar_data_asof <- ar_data_full_r %>%
      filter(date <= cutoff_dt) %>%
      filter(issue <= cutoff_dt) %>%
      group_by(date) %>%
      slice_max(issue) %>%
      ungroup()
  }
  
  if(model_type == "lvcf"){
    if(in_season == TRUE){
      
      in_season_out <- mk_in_season(ar_data_asof, NULL)
      ar_data_asof <- as.data.frame(in_season_out[1])
    }
    
    ar_X <- NULL
    ar_Y <- NULL
    
  } else{
  
  # for poisson, take log first, before doing diffs/lags
  if(model_type == "poisson"){
    c1 <- min(flu_data$ili[flu_data$ili > 0])/2
    ar_data_asof <- ar_data_asof %>%
      rowwise() %>%
      mutate(across(.cols = c({predictor}), ~logxplusc(.x, c1))) %>%
      ungroup()
    
  }
  # calculate lagged values for Primary Region
  ar_data_asof <- calculate_lags(ar_data_asof, {predictor}, lags, lag_diff = FALSE, self = TRUE, drop_cols=FALSE, indic = FALSE, tol = tol)
  
  names(ar_data_asof) <- gsub(paste0({predictor}, "_"), "", names(ar_data_asof))
  
  if(length(nbr_lst) > 0){

    nbr_lags <- do_nbr_lags(nbr_lst, {predictor}, lag_diff, model_type = model_type, as_of = as_of, cutoff_dt = cutoff_dt, indic = indic, tol = tol)
    
    if(in_season == TRUE){
      
      in_season_out <- mk_in_season(ar_data_asof, nbr_lags)
      ar_data_asof <- as.data.frame(in_season_out[1])
      nbr_lags <- as.data.frame(in_season_out[2])
      
    }
    
    if(model_type == "poisson"){
      lag_df <- ar_data_asof %>%
        dplyr::select(date, contains("lag"), issue, region) %>%
        left_join(nbr_lags %>%
                    dplyr::select(-WEEK), by = "date") %>%
        mutate(int = 1)
    }
    else{
      lag_df <- ar_data_asof %>%
      dplyr::select(date, contains("lag"), issue) %>%
      left_join(nbr_lags %>%
                  dplyr::select(-WEEK), by = "date") %>%
      mutate(int = 1)
      }
    
    ar_data_asof <- left_join(ar_data_asof, nbr_lags, by = c("date", "WEEK"))
    
  } else{ # if state has no neighbors / isolated model
    
    if(in_season == TRUE){
      
      in_season_out <- mk_in_season(ar_data_asof, NULL)
      ar_data_asof <- as.data.frame(in_season_out[1])
    }
    if(model_type == "poisson"){
    lag_df <- ar_data_asof %>%
      dplyr::select(date, contains("lag"), issue, region) %>%
      mutate(int = 1)
    } else{
      lag_df <- ar_data_asof %>%
        dplyr::select(date, contains("lag"), issue) %>%
        mutate(int = 1)
    }
  }

  
if(model_type == "poisson"){    
    
    ar_X <- lag_df %>% replace(is.na(.), 0) %>%
      dplyr::select(int, dplyr::starts_with("lag"), date, issue, region)

  } else{
    ar_X <- lag_df %>% replace(is.na(.), 0) %>%
      dplyr::select(int, dplyr::starts_with("lag"), date, issue)
  }
  
    
    
    ar_Y <- ar_data_asof%>%
      dplyr::select(c(date, {outcome}))

  }
  return(list(ar_X, ar_Y, ar_data_asof))
}

# runs for one state at a time

fit_models <- function(flu_data, lags = c(0,1,2), prim_reg = "PA",
                       model_type = "quantile", inc_nbrs = TRUE, neg_pred,
                       pred_start_dt = NA, in_season = TRUE, quantiles,
                      as_of, outcome = "ili", lag_diff, predictor = NULL,
                    indic = FALSE, tol = 0.05,
                       natl, pooled){
  

  # st_mat <- read.csv(paste0(dir_nm, "/", "states_mat_adj.csv"), header=TRUE, row.names=1)
  st_mat <- read.csv( "states_mat_adj.csv", header=TRUE, row.names=1)
  
  if(model_type == "lvcf"){
    
    
    ar_data_full_r <- flu_data %>%
      mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
      mutate(epiweek = date) %>%
      filter(region == prim_reg) %>%
      dplyr::select(c(region, date, epiweek, YEAR, WEEK, {outcome}, issue)) %>%
      mutate(issue = as.Date(issue, c("%m/%d/%y")))
  
    if(nrow(ar_data_full_r) == 0){
      return(NA)
    }
    
    
    # filter dates to be in season, but just to pick which dates to use in lapply
    # this is not the data fed to the regression models -- have to filter issues first
    
    if(in_season == TRUE){
      in_season_out <- mk_in_season(ar_data_full_r, NULL)
      ar_data_full_X <- as.data.frame(in_season_out[1])
    }
    
    ar_X <- ar_data_full_X %>%
      dplyr::select(date) %>%
      distinct()
    
    
    # set date to start doing predictions
    # do this after filtering to be in season

    if(is.na(pred_start_dt)){
      init_dt_n <- length(lags) + 12 * length(lags) + 1  + 3 # more dates than parameters
    } else{
      init_dt_n <- pmax(which.min(abs(as.Date(pred_start_dt) - ar_X$date)), length(lags) + 12 * length(lags) + 1 + 3)
    }
    
    data_dts <- ar_X$date[order(ar_X$date)]
    
    n_beta <- 1
    
    num_cols <- 5 + n_beta
    
    model_df <- matrix(NA_real_, nrow = 1, ncol=num_cols) %>%
      as.data.frame() %>%
      dplyr::rename("prediction_date" = V1,
                    "target_date" = V2,
                    "quantile" = V3,
                    "mse" = V4,
                    "mad" = V5,
                    "beta_0" = V6)
    
    model_df <- model_df %>%
      mutate(prediction = NA_real_)
    
    nbr_lst <- list()
    
    if(as_of == FALSE){
      cutoff_last <- max(ar_data_full_r$date)
      ar_lst <- arrange_data(ar_data_full_r, nbr_lst, cutoff_dt = cutoff_last,
                             lags, in_season, as_of, outcome, lag_diff, predictor,
                             model_type, indic, natl)
    } else{
      ar_lst <- NULL
    }
    
    
    r <- lapply(data_dts[init_dt_n:length(data_dts)], get_lvcf_models, 
                ar_data_full_r,
                in_season = in_season,
                as_of = as_of,
                model_df_i = model_df,
                outcome = {outcome},
                ar_lst = ar_lst)
    
    r2 <- as.data.frame(do.call(rbind, r))
    
    ar_data_full_r_out <- arrange_data(ar_data_full_r, nbr_lst, max(data_dts), lags, 
                                       in_season, as_of, outcome, lag_diff, predictor,
                                       model_type, indic = indic, natl) 
    
    ar_data_full_r_filtered <- as.data.frame(ar_data_full_r_out[3])
    
    full_df <- r2 %>%
      left_join(ar_data_full_r_filtered, join_by("target_date" == "date")) 
    
    full_df$nbrs <- NA
    full_df$model_type <- model_type
    full_df$WIS <- NA
 
  } # end lvcf
  
  # do the data process but only to get the necessary dates, then arrange data inside regression function
  if(model_type == "quantile" | model_type == "linear"){
  
    if(pooled == TRUE){ # use all states
      inc_nbrs <- FALSE
      if(is.null(predictor) | predictor == outcome){
        ar_data_full_r <- flu_data %>%
          mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
          mutate(epiweek = date) %>%
          dplyr::select(c(region, date, epiweek, YEAR, WEEK, {outcome}, issue)) %>%
          mutate(issue = as.Date(issue, c("%m/%d/%y"))) 
      } else{
        ar_data_full_r <- flu_data %>%
          mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
          mutate(epiweek = date) %>%
          dplyr::select(c(region, date, epiweek, YEAR, WEEK, {outcome}, issue, {predictor})) %>%
          mutate(issue = as.Date(issue, c("%m/%d/%y"))) 
      }
    } else{
      if(is.null(predictor) | predictor == outcome){
        ar_data_full_r <- flu_data %>%
          mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
          mutate(epiweek = date) %>%
          filter(region == prim_reg) %>%
          dplyr::select(c(region, date, epiweek, YEAR, WEEK, {outcome}, issue)) %>%
          mutate(issue = as.Date(issue, c("%m/%d/%y")))
      } else{
        ar_data_full_r <- flu_data %>%
          mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
          mutate(epiweek = date) %>%
          filter(region == prim_reg) %>%
          dplyr::select(c(region, date, epiweek, YEAR, WEEK, {outcome}, issue, {predictor})) %>%
          mutate(issue = as.Date(issue, c("%m/%d/%y")))
      }
    }
 
    if(nrow(ar_data_full_r) == 0){
      return(NA)
    }
    
    nbr_lst <- get_nbrs(prim_reg, inc_nbrs, natl)
    
    # filter dates to be in season, but just to pick which dates to use in lapply
    # this is not the data fed to the regression models -- have to filter issues first
      
      if(in_season == TRUE){
        in_season_out <- mk_in_season(ar_data_full_r, NULL)
        ar_data_full_X <- as.data.frame(in_season_out[1])
      }
      
    ar_X <- ar_data_full_X %>%
      dplyr::select(date) %>%
      distinct()

    
    # set date to start doing predictions
    # do this after filtering to be in season
    
    # no state has more than 8 neighbors... plus 4 influencers... 12?
    if(is.na(pred_start_dt)){
      init_dt_n <- length(lags) + 12 * length(lags) + 1 + 3 # more dates than parameters
    } else{
      init_dt_n <- pmax(which.min(abs(as.Date(pred_start_dt) - ar_X$date)), length(lags) + 12 * length(lags) + 1 + 3)
    }
    
    data_dts <- ar_X$date[order(ar_X$date)]
  
  if(model_type == "quantile") {
     # initialize df to store predictions and loss
    
    n_beta <- if_else(length(nbr_lst) > 0, length(nbr_lst) * length(lags), 0) + length(lags) + 1 
    
    if(lag_diff == TRUE & indic == TRUE & sum(lags == 0) > 0){
      n_beta <- n_beta - length(nbr_lst) # remove 0wk cols
    } 
    
    num_cols <- 5 + n_beta
    
    
    model_df <- matrix(NA_real_, nrow = length(quantiles), ncol=num_cols) %>%
      as.data.frame() %>%
      dplyr::rename("prediction_date" = V1,
             "target_date" = V2,
             "quantile" = V3,
             "mse" = V4,
             "mad" = V5,
             "beta_0" = V6) %>%
      mutate(prediction_date = as.Date(prediction_date),
             target_date = as.Date(target_date))
    
    self_cols <- paste0("lag_", lags, "wk")
    
    if(length(nbr_lst) > 0){
      if(lag_diff == TRUE & indic == TRUE){ # if using indicators, don't use current week for neighbors
        lag_pos <- lags[which(lags > 0)]
        nbr_cols <- paste0("lag_", rep(nbr_lst,each=length(lag_pos)), "_", lag_pos, "wk")
      } else{
        nbr_cols <- paste0("lag_", rep(nbr_lst,each=length(lags)), "_", lags, "wk")
      }
    } else{
      nbr_cols <- NULL
    }
    
    v <- c(self_cols, nbr_cols)
    
    nms2 <- gsub("lag", "beta", v)
    colnames(model_df)[7:(length(colnames(model_df)))]<-nms2
    # colnames(model_df)[length(colnames(model_df))]<-"prediction"
    
    if(as_of == FALSE){
      cutoff_last <- max(ar_data_full_r$date)
      ar_lst <- arrange_data(ar_data_full_r, nbr_lst, cutoff_dt = cutoff_last,
                             lags, in_season, as_of, outcome, lag_diff, predictor,
                             model_type, indic, natl)
    } else{
      ar_lst <- NULL
    }
    
    # fit quantile regression models at each quantile
    if(model_type == "quantile"){
    r <- lapply(data_dts[init_dt_n:length(data_dts)], get_qx_models, 
                ar_data_full_r,
                nbr_lst,
                quantiles=quantiles, 
                model_df_i = model_df, 
                neg_pred = neg_pred, 
                outcome = outcome,
                as_of = as_of,
                lags = lags,
                in_season = in_season,
                lag_diff = lag_diff,
                ar_lst = ar_lst,
                predictor = predictor,
                prim_reg = prim_reg)
    } 
    r2 <- as.data.frame(do.call(rbind, r))
    
    ar_data_full_r_out <- arrange_data(ar_data_full_r, nbr_lst, max(data_dts), lags, 
                                            in_season, as_of, outcome, lag_diff, predictor,
                                       model_type, indic = indic, natl) 
    
    ar_data_full_r_filtered <- as.data.frame(ar_data_full_r_out[3])
    
    ar_X <- as.data.frame(ar_data_full_r_out[1])
 
    if(inc_nbrs == TRUE & length(nbr_lst) > 0){
      full_data_tbl <- ar_data_full_r_filtered %>%
        filter(region == prim_reg) %>%
        dplyr::select(-starts_with("lag")) %>%
        left_join(ar_X, by = "date")  
    } else {
      full_data_tbl <- ar_data_full_r_filtered %>%
        filter(region == prim_reg)
    }
    
    if(pooled == TRUE){
      
      r3 <- sort_quantiles(r2)
      
      full_df <- r3 %>%
        left_join(full_data_tbl, join_by("target_date" == "date")) 
      
      full_df$nbrs <- inc_nbrs
      full_df$model_type <- model_type
      
      full_df$WIS <- NA
      
      full_df <- as.data.frame(full_df)
    } else{
      full_df <- r2 %>%
        left_join(full_data_tbl, join_by("target_date" == "date")) 
      
      full_df$nbrs <- inc_nbrs
      full_df$model_type <- model_type
      
      full_df <- sort_quantiles(full_df) 
      full_df <- add_wis(full_df, quantiles, outcome, "quantile")
    }
   
    
  
    } else if(model_type == "linear") {
    
    # initialize df to store predictions and loss
    
    n_beta <- if_else(length(nbr_lst) > 0, length(nbr_lst) * length(lags), 0) + length(lags) + 1 
    
    num_cols <- 5 + n_beta
    
    model_df <- matrix(NA_real_, nrow = 1, ncol=num_cols) %>%
      as.data.frame() %>%
      dplyr::rename("prediction_date" = V1,
             "target_date" = V2,
             "se" = V3,
             "mse" = V4,
             "mad" = V5,
             "beta_0" = V6) %>%
      mutate(prediction_date = as.Date(prediction_date),
             target_date = as.Date(target_date))
   
     self_cols <- paste0("lag_", lags, "wk")
    
    if(length(nbr_lst) > 0){
      nbr_cols <- paste0("lag_", rep(nbr_lst,each=length(lags)), "_", lags, "wk")
    } else{
      nbr_cols <- NULL
    }
    
    v <- c(self_cols, nbr_cols)
    
    nms2 <- gsub("lag", "beta", v)
    colnames(model_df)[7:(length(colnames(model_df)))]<-nms2
    
    model_df <- model_df %>%
      mutate(prediction = NA_real_,
             se_df = NA_real_)
    
    if(as_of == FALSE){
      cutoff_last <- max(ar_data_full_r$date)
      ar_lst <- arrange_data(ar_data_full_r, nbr_lst, cutoff_dt = cutoff_last, lags,
                             in_season, as_of, outcome, lag_diff, predictor, model_type,
                             indic = indic, natl)
      } else{
      ar_lst <- NULL
    }
    
    # fit linear regression models
    r <- lapply(data_dts[init_dt_n:length(data_dts)], get_lin_models, 
                 ar_data_full_r,
                nbr_lst,
                quantiles=quantiles, 
                model_df_i = model_df, 
                neg_pred = neg_pred, 
                outcome = {outcome},
                as_of = as_of,
                lags = lags,
                in_season = in_season,
                lag_diff = lag_diff,
                ar_lst = ar_lst,
                predictor = predictor,
                prim_reg = prim_reg)
    
    r2 <- as.data.frame(do.call(rbind, r))
    
    ar_data_full_r_out <- arrange_data(ar_data_full_r, nbr_lst, max(data_dts), lags, 
                                       in_season,
                                       as_of, outcome, lag_diff, predictor, model_type, indic = indic,
                                       natl) 
    
    ar_data_full_r_filtered <- as.data.frame(ar_data_full_r_out[3])
    
    ar_X <- as.data.frame(ar_data_full_r_out[1])
    
    
    if(inc_nbrs == TRUE & length(nbr_lst) > 0){
      full_data_tbl <- ar_data_full_r_filtered %>%
        dplyr::select(-starts_with("lag")) %>%
        left_join(ar_X, by = "date")  
    } else {
      full_data_tbl <- ar_data_full_r_filtered
    }

    full_df <- r2 %>%
      left_join(full_data_tbl, join_by("target_date" == "date")) 
    
    full_df$nbrs <- inc_nbrs
    full_df$model_type <- model_type
    
    full_df$int <- 1
    
     quantiles <- as.data.frame(quantiles)
    names(quantiles) = "quantile"
    
    lin_q <- full_df %>%
      cross_join(quantiles)
    
    full_df <- add_wis(lin_q, quantiles, outcome, model_type)  
    
  ################################# Poisson ##################################
    
  }
  }  else if(model_type == "poisson"){
    
    if(pooled == TRUE){
      inc_nbrs <- FALSE
      
      if(is.null(predictor) | predictor == outcome){
        ar_data_full_r <- flu_data %>%
          mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
          mutate(epiweek = date) %>%
          dplyr::select(c(region, date, epiweek, YEAR, WEEK, {outcome}, issue, num_patients, num_providers, ili)) %>%
          mutate(issue = as.Date(issue, c("%m/%d/%y")))
      } else{
        ar_data_full_r <- flu_data %>%
          mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
          mutate(epiweek = date) %>%
          dplyr::select(c(region, date, epiweek, YEAR, WEEK, {outcome}, issue, {predictor}, num_patients, num_providers)) %>%
          mutate(issue = as.Date(issue, c("%m/%d/%y")))
      }
    } else{
      
      if(is.null(predictor) | predictor == outcome){
        ar_data_full_r <- flu_data %>%
          mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
          mutate(epiweek = date) %>%
          filter(region == prim_reg) %>%
          dplyr::select(c(region, date, epiweek, YEAR, WEEK, {outcome}, issue, num_patients, num_providers, ili)) %>%
          mutate(issue = as.Date(issue, c("%m/%d/%y")))
      } else{
        ar_data_full_r <- flu_data %>%
          mutate(date = as.Date(epiweek, c("%m/%d/%y"))) %>%
          mutate(epiweek = date) %>%
          filter(region == prim_reg) %>%
          dplyr::select(c(region, date, epiweek, YEAR, WEEK, {outcome}, issue, {predictor}, num_patients, num_providers)) %>%
          mutate(issue = as.Date(issue, c("%m/%d/%y")))
      }
    }
    
    
    
    ar_data_filled <- ar_data_full_r %>%
      arrange(region, date) %>%
      mutate(across(.cols = c({predictor}, {outcome}, num_patients, num_providers), .fns = ~if_else(num_patients == 0, dplyr::lag(.x, 1), .x)))

    if(nrow(ar_data_filled) == 0){
      return(NA)
     }
     # convert ili to log ili
    c1 <- min(flu_data$ili[flu_data$ili > 0])/2

    nbr_lst <- get_nbrs(prim_reg, inc_nbrs, natl)
    
    if(in_season == TRUE){
      in_season_out <- mk_in_season(ar_data_filled, NULL)
      ar_data_full_X <- as.data.frame(in_season_out[1])
    }
    
    ar_X <- ar_data_full_X %>%
      dplyr::select(date) %>%
      distinct()
    
    if(is.na(pred_start_dt)){
      init_dt_n <- length(lags) + 12 * length(lags) + 1 + 3 # more dates than parameters
    } else{
      init_dt_n <- pmax(which.min(abs(as.Date(pred_start_dt) - ar_X$date)), length(lags) + 12 * length(lags) + 1 + 3)
    }
    
    data_dts <- ar_X$date[order(ar_X$date)]
    
    n_beta <- if_else(length(nbr_lst) > 0, length(nbr_lst) * length(lags), 0) + length(lags) + 1
    
    if(lag_diff == TRUE & indic == TRUE & sum(lags == 0) > 0){
      n_beta <- n_beta - length(nbr_lst) # remove 0wk cols
    } 
    num_cols <- 6 + n_beta
    
    model_df <- matrix(NA_real_, nrow = 1, ncol=num_cols) %>%
      as.data.frame() %>%
      dplyr::rename("prediction_date" = V1,
                    "target_date" = V2,
                    "se" = V3,
                    "mse" = V4,
                    "mad" = V5,
                    "beta_0" = V6) %>%
      mutate(prediction_date = as.Date(prediction_date),
             target_date = as.Date(target_date))
    
    self_cols <- paste0("lag_", lags, "wk")
    
    if(length(nbr_lst) > 0){
      if(lag_diff == TRUE & indic == TRUE){ # if using indicators, don't use current week for neighbors
        lag_pos <- lags[which(lags > 0)]
        nbr_cols <- paste0("lag_", rep(nbr_lst,each=length(lag_pos)), "_", lag_pos, "wk")
      } else{
        nbr_cols <- paste0("lag_", rep(nbr_lst,each=length(lags)), "_", lags, "wk")
      }
    } else{
      nbr_cols <- NULL
    }
    
    v <- c(self_cols, nbr_cols)
    
    nms2 <- gsub("lag", "beta", v)
    nms3 <- c(nms2, "beta_num_providers")
    colnames(model_df)[7:(length(colnames(model_df)))]<-nms3
    
    model_df <- model_df %>%
      mutate(prediction = NA_real_,
             se_df = NA_real_)
    
    if(as_of == FALSE){
      cutoff_last <- max(ar_data_full_r$date)
      ar_lst <- arrange_data(ar_data_filled, nbr_lst, cutoff_dt = cutoff_last, lags,
                             in_season, as_of, outcome, lag_diff, predictor, model_type,
                             indic = indic, natl)
    } else{
      ar_lst <- NULL
    }

    # fit poisson regression models
    
    # r <- lapply(data_dts[init_dt_n:length(data_dts)], get_poisson_models, 
    #             ar_X=ar_X, 
    #             ar_Y=ar_Y, 
    #             model_df_i = model_df, 
    #             neg_pred = neg_pred,
    #             regularization = regularization)
    
    r <- lapply(data_dts[init_dt_n:length(data_dts)], get_poisson_models, 
                ar_data_filled,
                nbr_lst,
                model_df_i = model_df, 
                neg_pred = neg_pred, 
                outcome = {outcome},
                as_of = as_of,
                lags = lags,
                in_season = in_season,
                lag_diff = lag_diff,
                ar_lst = ar_lst,
                predictor = {predictor},
                prim_reg = prim_reg)
    
    r2 <- as.data.frame(do.call(rbind, r))
    
    # if(length(nbr_lst) > 0){
    #   full_data_tbl <- ar_data_filled %>%
    #     dplyr::select(-starts_with("lag")) %>%
    #     left_join(
    #       (ar_X %>% dplyr::select(-c(num_patients, num_providers))), 
    #        by = "date")  
    # } else {
    #   full_data_tbl <- ar_data_filled
    # }
    
    ar_data_full_r_out <- arrange_data(ar_data_filled, nbr_lst, max(data_dts), lags, 
                                       in_season,
                                       as_of, outcome, lag_diff, predictor, model_type, indic, natl) 
    
    ar_data_full_r_filtered <- as.data.frame(ar_data_full_r_out[3])
    
    ar_X <- as.data.frame(ar_data_full_r_out[1]) %>%
      dplyr::select(-issue)
    
    if(inc_nbrs == TRUE & length(nbr_lst) > 0){
      full_data_tbl <- ar_data_full_r_filtered %>%
        filter(region == prim_reg) %>%
        dplyr::select(-starts_with("lag")) %>%
        left_join(ar_X, by = c("date", "region"))  
    } else {
      full_data_tbl <- ar_data_full_r_filtered %>%
        filter(region == prim_reg)
    }
    
    full_data_tbl <- full_data_tbl %>%
      mutate(across(.cols = {predictor}, .fns = ~exp(.x + c1)))
    
    
    if(pooled == TRUE){
    
      full_df <- r2 %>%
        left_join(full_data_tbl, join_by("target_date" == "date")) 
      
      full_df$nbrs <- inc_nbrs
      full_df$model_type <- model_type
      
      # full_df$prediction <- predict_poisson(full_df, inc_hum = FALSE)
      
      quantiles <- as.data.frame(quantiles)
      names(quantiles) = "quantile"
      
      poi_q <- full_df %>%
        cross_join(quantiles)
      
      full_df$WIS <- NA
      full_df <- as.data.frame(full_df)
      
    } else{
      full_df <- r2 %>%
        left_join(full_data_tbl, join_by("target_date" == "date")) 
      
      full_df$nbrs <- inc_nbrs
      full_df$model_type <- model_type
      
      # full_df$prediction <- predict_poisson(full_df, inc_hum = FALSE)
      
      quantiles <- as.data.frame(quantiles)
      names(quantiles) = "quantile"
      
      
      poi_q <- full_df %>%
        cross_join(quantiles)
      
      full_df <- add_wis(poi_q, quantiles, "ili", "poisson")
    }
    
  }
  
  return(full_df)
}
