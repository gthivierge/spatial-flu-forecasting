#!/usr/bin/env

# read in csv file
# call function to fit models
# make the model summaries
# store the model summaries

#' @param file File with combinations of model specs 
#' @export

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
cores <- args[2]

if(!is.na(args[3])){
  start_row <- args[3]
} else{
  start_row <- 1
}

if(!is.na(args[4])){
  end_row <- args[4]
} else{
  end_row <- Inf
}

check_pkg <- function(pkg){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
  } 
}

check_pkg("import")
check_pkg("plyr")
check_pkg("dplyr")
check_pkg("tidyr")
check_pkg("ggplot2")
check_pkg("lubridate")
check_pkg("readr")
check_pkg("purrr")
check_pkg("stringr")
check_pkg("quantreg")
check_pkg("kit")
check_pkg("data.table")
check_pkg("parallel")
check_pkg("epidatr")
check_pkg("glmnet")
check_pkg("hqreg")
check_pkg("scoringutils")
check_pkg("limSolve")
# check_pkg("MMWRweek")

# library(plyr, quietly = TRUE)

import::from(plyr, adply)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(lubridate, quietly = TRUE)
library(readr, quietly = TRUE)
library(purrr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(quantreg, quietly = TRUE)
library(kit, quietly = TRUE)
library(data.table, quietly = TRUE)
library(parallel, quietly = TRUE)
library(epidatr, quietly = TRUE)
library(glmnet, quietly = TRUE)
library(hqreg, quietly = TRUE)
library(scoringutils, quietly = TRUE)
library(limSolve, quietly = TRUE)


source("reg_fns.R")
source("calc_fns.R")
source("interval_fns.R")
source("do_model_fit_fn.R")
source("quantile_tracker.R")

# file <- "pooled_v2.csv"
# cores <- 6
# start_row <- 17
 
newfolder <- paste0(gsub(".csv", "", file), "_", today())
 
dir.create(file.path(getwd(), newfolder), showWarnings = FALSE) 
 
model_input <- read.csv(file, header=TRUE)


flu_data <- read.csv("flu_data_all.csv", header = TRUE)

ili_dat <- flu_data %>%
  group_by(region, epiweek) %>%
  slice_max(issue) %>%
  ungroup() %>%
  dplyr::select(epiweek, region, ili, wili) %>%
  mutate(epiweek = as.Date(epiweek, c("%m/%d/%y"))) %>%
  arrange(region, epiweek) %>% 
  mutate(ili = if_else(region == "US", wili, ili)) %>%
  dplyr::select(-wili)

st_mat <- read.csv("states_mat_adj.csv", header=TRUE, row.names=1)

# quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
               # 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99) # from cmu flusights submission

# quantiles <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
quantiles <- c(0.025, 0.5, 0.975)


model_summary <- as.data.frame(matrix(nrow = 0, ncol = 37))
names(model_summary) <- c("target_date",    "epiweek",        "region",         "prediction" ,    "quantile",       "WIS",           
                          "mse",            "mad",            "se",             "se_df",          "nbrs",           "model_type", 
                          "ili", "WIS.conf", "abserr", "sqerr", "cmad", "cmse",
                          "lags",          "neg_pred",       "in_season",      
                          "lag_diff",       "outcome",       "as_of", "predictor",
                           "indic", "tol", "natl", "pooled")         

model_summary$target_date <- as.Date(model_summary$target_date)
model_summary[,which(names(model_summary) %in% c("region", "model_type"))] <- lapply(model_summary[,which(names(model_summary) %in% c("region", "model_type"))], as.character)
model_summary[,which(!names(model_summary) %in% c("region", "model_type", "target_date", "nbrs"))] <- lapply(model_summary[,which(!names(model_summary) %in% c("region", "model_type", "target_date", "nbrs"))], as.numeric)

# make file for storage
new_sum_file <- paste0(newfolder, "/", substr(file, 1, nchar(file) - 4), "_summary", today(),"_", format(Sys.time(), "%H.%M.%S"), ".csv")

fwrite(model_summary, new_sum_file)

# print(paste("starting model fitting", now()))

# end_row = 16
# start_row = 9
for(i in start_row:pmin(end_row, nrow(model_input))){
  
  print(paste("starting row", i, now()))

  
  if(model_input$prim_reg[i] == "ALL"){
    prim_reg <- as.list(names(st_mat))
  } else{
    prim_reg <- list(model_input$prim_reg[i])
  }
  
   # prim_reg <- "NE"

  lags <- c(0:model_input$lags[i])
  model_type <- model_input$model_type[i]
  inc_nbrs <- model_input$inc_nbrs[i]
  neg_pred <- model_input$neg_pred[i]
  pred_start_dt <- as.character(as.Date(model_input$pred_start_dt[i], format = "%m/%d/%y"))
  in_season <- model_input$in_season[i]
  lag_diff <- model_input$lag_diff[i] # for doing differencing
  outcome <- model_input$outcome[i] # ILI or ILI+
  as_of <- model_input$as_of[i]
  predictor <- model_input$predictor[i]
  indic <- model_input$indic[i]
  tol <- model_input$tol[i]
  natl <- model_input$natl[i]
  pooled <- model_input$pooled[i]
  
  
  big_df <- make_big_df(flu_data, prim_reg, lags, model_type, inc_nbrs,
  neg_pred, pred_start_dt, in_season, quantiles, as_of, 
  outcome, lag_diff, predictor, indic, tol, natl,
  pooled, cores)

  
    print("writing big file...")
    
    new_file <- paste0(newfolder, "/", substr(file, 1, nchar(file) - 4), "_row", i, "_", model_type, "_detail.csv.gz")
    
    fwrite(big_df, new_file)
 
  # save summary statistics (small df) separately from all model parameters (big df) 
  small_df_0 <- make_small_df(big_df, model_type)
  
  
  
  # add conformal WIS automatically
  
  
  conf_intervals <-c(0.5, 0.8, 0.95)

  
  df_conf <- small_df_0 %>%
    ungroup() %>%
    left_join(ili_dat, by = join_by(epiweek, region)) %>%
    group_by(region) %>%
    filter(quantile == 0.5 | is.na(quantile))
  
  conf_lst <- group_split(df_conf)
  conf_out <- mclapply(conf_lst, quantile_tracker, conf_intervals, eta = 0.01, mc.cores = 6)
  
  # test_df <- conf_lst[[1]]
  
  
  small_df <- do.call(rbind, conf_out)
  
  # y <- quantile_tracker(test_df, conf_intervals, 0.01)
  
  small_df$lags <- max(lags)
  small_df$model_type <- model_type
  small_df$neg_pred <- neg_pred
  small_df$in_season <- in_season
  small_df$gamma_t <- gamma_t
  small_df$lag_diff <- lag_diff # for doing differencing
  small_df$outcome <- outcome # ILI or ILI+
  small_df$as_of <- as_of
  small_df$predictor <- predictor
  small_df$indic <- indic
  small_df$tol <- tol
  small_df$natl <- natl
  small_df$pooled <- pooled
  
  
  # model_summary <- bind_rows(model_summary, small_df)
  fwrite(small_df, new_sum_file, append = TRUE)
  
  print(paste("completed row", i, now()))
}

# print(paste("writing model summary file...", now()))

