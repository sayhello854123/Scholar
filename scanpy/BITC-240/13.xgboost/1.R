gc()
library(survival)
library(survminer)
library(survXgboost)
library(xgboost)

head(lung)
lung <- lung[complete.cases(lung), ] # doesn't handle missing values at the moment
lung$status <- lung$status - 1 # format status variable correctly such that 1 is event/death and 0 is censored/alive
label <- ifelse(lung$status == 1, lung$time, -lung$time)

val_ind <- sample.int(nrow(lung), 0.1 * nrow(lung))
x_train <- as.matrix(lung[-val_ind, !names(lung) %in% c("time", "status")])
x_label <- label[-val_ind]
x_val <- xgb.DMatrix(as.matrix(lung[val_ind, !names(lung) %in% c("time", "status")]),
                     label = label[val_ind])
a <- lung[val_ind, !names(lung) %in% c("time", "status")]
a <- as.numeric(a)
A1 <- Matrix(a,sparse = T)
surv_xgboost_model <- xgb.train.surv(
  params = list(
    objective = "survival:cox",
    eval_metric = "cox-nloglik",
    eta = 0.05 # larger eta leads to algorithm not converging, resulting in NaN predictions
  ), data = x_train, label = x_label,
  watchlist = list(val2 = x_val),
  nrounds = 1000, early_stopping_rounds = 30
)
importance <- xgb.importance(colnames(x_train), model = surv_xgboost_model)  

times <- seq(10, 1000, 50)
survival_curves <- predict(object = surv_xgboost_model, newdata = x_train, type = "surv", times = times)
dim(survival_curves)
matplot(times, t(survival_curves[1:5, ]), type = "l")
risk_scores <- predict(object = surv_xgboost_model, newdata = x_train, type = "risk")
hist(risk_scores)
importance <- xgb.importance(colnames(a), model = surv_xgboost_model)  
head(importance)