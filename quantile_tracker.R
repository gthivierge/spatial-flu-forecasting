
# Quantile tracker
# See Anastasios N. Angelopoulos, Emmanuel J. Candès, and Ryan J. Tibshirani. Conformal PID Control for Time
# Series Prediction. Advances in Neural Information Processing Systems, 36, 7 2023.


quantile_update <- function(st, qt, tau, eta){
  alpha<- 1-tau
  err <- (st>qt)
  qt+eta*(err-alpha)
}


quantile_tracker <- function(dfr, intervals, eta){

tau <- intervals
  
N <- nrow(dfr)
# eta <- 0.01

dfr$score <- abs(dfr$ili - dfr$prediction)

IS_all <- as.data.frame(cbind(dfr$ili, dfr$prediction))
names(IS_all) <- c("ili", "median")

for(i in 1:length(tau)){
  q_est<-matrix(NA_real_,nrow=N,ncol=1)
  q_est[1,]<-dfr$prediction[1]
  for(t in 1:(N-1)){
    q_est[t+1,]<-quantile_update(dfr$score[t],q_est[t,],tau[i],eta)
  }
  
  interval_low<- dfr$prediction - q_est 
  interval_hi<- dfr$prediction + q_est
  
  
  IS <- scoringutils::interval_score(dfr$ili, interval_low, interval_hi, tau[i] * 100, weigh = FALSE)
  IS_all <- cbind(IS_all, IS)
  names(IS_all)[length(IS_all)] <- paste0("IS_", tau[i])
}

q_lower <- (1-tau)/2

WIS <- weighted_IS(q_lower, IS_all, "ili")

newdf <- bind_cols(dfr, "WIS.conf" = WIS) %>%
  mutate(abserr = abs(ili - prediction),
sqerr = (ili - prediction)^2) %>%
  mutate(cmad = cummean(abserr),
         cmse = cummean(sqerr)) %>%
  dplyr::select(-score)

return(newdf)
}