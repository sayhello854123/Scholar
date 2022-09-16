library('survival')
library('xgboost')
library('gbm')
# set random state
set.seed(0)

data(veteran, package = "randomForestSRC")
cat("Number of samples:", nrow(veteran), "\n")
cat("Columns of dataset:", colnames(veteran), "\n")
veteran[c(1:5), ]

x_cols <- c('trt', 'celltype', 'karno', 'diagtime', 'age', 'prior')
y_cols <- c('T')
train <- sample(1:nrow(veteran), round(nrow(veteran) * 0.80))

data_train <- veteran[train, ]
data_test <- veteran[-train, ]


data_train$T <- data_train$time
data_train$T[data_train$status==0] = -data_train$T[data_train$status==0]
data_test$T <- data_test$time
data_test$T[data_test$status==0] = -data_test$T[data_test$status==0]

Convert to xgb.DMatrix
dtrain <- xgb.DMatrix(as.matrix(data_train[x_cols]), label=data_train[, y_cols])
dtest  <- xgb.DMatrix(as.matrix(data_test[x_cols]), label=data_test[, y_cols])
a <- Matrix(data_train[x_cols],sparse = T)
set.seed(2891)
param <- list(max_depth=3, eta = 0.06, silent = 1, objective = "survival:cox", eval_metric = "cox-nloglik")
model <- xgb.train(param, dtrain, nrounds=100)
importance <- xgb.importance(colnames(data_train[x_cols]), model = model)  
head(importance)

pred.train <- log(predict(model, dtrain))
pred.test  <- log(predict(model, dtest))

Hmisc::rcorr.cens(-pred.train, Surv(data_train$time, data_train$status))
importance <- xgb.importance(colnames(), model = xgb)  