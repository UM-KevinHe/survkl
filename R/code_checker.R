#--------------------------- step by step check ------------------------------#
source("R/utils.R")

source("R/coxkl.R")
source("R/cv.coxkl.R")

source("R/coxkl_ridge.R")
source("R/cv.coxkl_ridge.R")

source("R/coxkl_highdim.R")
source("R/cv.coxkl_highdim.R")

Rcpp::sourceCpp("src/utils.cpp")
Rcpp::sourceCpp("src/KLCox.cpp")
Rcpp::sourceCpp("src/KLCox_highdim.cpp")


#-------------------- 1. check: compare with coxph under eta = 0 --------------------------# 
load("data/ExampleData.RData")
dat <- data.frame(time = ExampleData$time,
                  status = ExampleData$status,
                  stratum = ExampleData$stratum,
                  ExampleData$z)

library(survival)
rhs <- paste(colnames(ExampleData$z), collapse = " + ")
coxph_res <- coxph(
  as.formula(paste("Surv(time, status) ~", rhs, "+ strata(stratum)")),
  data = dat,
  control = coxph.control(eps = 1e-7, iter.max = 100)
)
coxph_res$coefficients
#        Z1         Z2         Z3         Z4         Z5         Z6 
# 0.1586184 -0.2285272  0.4641955 -0.3864544  0.3095823 -0.2555455 
coxph_res$loglik[2]
# [1] -1202.085



result1 <- coxkl(z = ExampleData$z,
                 delta = ExampleData$status,
                 time = ExampleData$time,
                 stratum = ExampleData$stratum,
                 RS = NULL,
                 beta = ExampleData$beta_external,
                 etas = 0)
result1$beta
#               0
# [1,]  0.1586184
# [2,] -0.2285272
# [3,]  0.4641955
# [4,] -0.3864544
# [5,]  0.3095823
# [6,] -0.2555455
result1$likelihood
#         0 
# -1202.085 


#-------------------- 2. check: compare with Lingfeng's original codes --------------------------# 
load("data/ExampleDataHighDim.RData")
result_new <- coxkl(z = ExampleData$z,
                    delta = ExampleData$status,
                    time = ExampleData$time,
                    RS = NULL,
                    beta = ExampleData$beta_external,
                    etas = seq(0, 5, 1))
head(result_new$beta)  # the result is the same with that form Lingfeng's original code
#                0           1           2           3           4           5
# [1,]  0.33158984  0.30507909  0.29725498  0.29348346  0.29126065  0.28979449
# [2,] -0.45556536 -0.38588179 -0.36236166 -0.35044965 -0.34324199 -0.33840898
# [3,] -0.07571637 -0.04594781 -0.03470755 -0.02884734 -0.02525266 -0.02282281
# [4,] -0.07332686  0.00575570  0.03151854  0.04437497  0.05208957  0.05723476
# [5,]  0.62690796  0.54722031  0.52055560  0.50720223  0.49918453  0.49383736
# [6,]  0.19092080  0.09606307  0.06455724  0.04881102  0.03936478  0.03306787




dat <- data.frame(time = ExampleData$time,
                  status = ExampleData$status,
                  stratum = ExampleData$stratum,
                  ExampleData$z)

rhs <- paste(colnames(ExampleData$z), collapse = " + ")
coxph_res <- coxph(
  as.formula(paste("Surv(time, status) ~", rhs)),
  data = dat,
  control = coxph.control(eps = 1e-7, iter.max = 100)
)
head(coxph_res$coefficients)
#         Z1          Z2          Z3          Z4          Z5          Z6 
# 0.33158984 -0.45556536 -0.07571637 -0.07332686  0.62690796  0.19092080 

#-------------------- 3. check: cv.coxkl --------------------------# 
load("data/ExampleData.RData")

etas <- generate_eta(method = "exponential", n = 10, max_eta = 5)
etas <- sample(etas) ## suffle eta orders
cv_res1 <- cv.coxkl(z = ExampleData$z,
                    delta = ExampleData$status,
                    time = ExampleData$time,
                    stratum = NULL,
                    RS = NULL,
                    beta = ExampleData$beta_external,
                    etas = etas,
                    nfolds = 5,
                    criteria = c("LinPred"),   #"V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"
                    message = T)
cv_res1
# $internal_stat
#           eta LinPred_Loss
# 1  0.00000000     3534.119
# 2  0.03374245     3533.104
# 3  0.09002825     3531.598
# 4  0.18391863     3529.504
# 5  0.34053721     3526.852
# 6  0.60179276     3523.880
# 7  1.03759328     3520.997
# 8  1.76455236     3518.588
# 9  2.97719318     3516.823
# 10 5.00000000     3515.652
# 
# $external_stat
# [1] 3513.818

plot(cv_res1$internal_stat[,1], cv_res1$internal_stat[,2], type = "b", xlab = "eta")



#---------------------------4. ridge (coxkl_ridge)------------------------------#
load("data/ExampleDataHighDim.RData")
result3 <- coxkl_ridge(z = ExampleData$z,
                       delta = ExampleData$status,
                       time = ExampleData$time,
                       stratum = NULL,
                       RS = NULL,
                       beta = ExampleData$beta_external,
                       message = T)
#   Cross-validation over lambda sequence:
#   |==============================| 100%

head(result3$beta[, 1:5])
#         138.1936    125.9168    114.7307    104.5384     95.2515
# [1,]  0.13888135  0.14639228  0.15397647  0.16160406  0.16924477
# [2,] -0.17547848 -0.18474565 -0.19417029 -0.20372421 -0.21337737
# [3,] -0.04192810 -0.04351869 -0.04505306 -0.04652448 -0.04792778
# [4,]  0.03213009  0.03161102  0.03088162  0.02993718  0.02877540
# [5,]  0.23628338  0.24810751  0.26010622  0.27225002  0.28450815
# [6,]  0.09007171  0.09407818  0.09810041  0.10212425  0.10613519

#---------------------------5. cv.coxkl_ridge ------------------------------#
load("data/ExampleDataHighDim.RData")
set.seed(1)
etas <- generate_eta(method = "exponential", n = 50, max_eta = 10)
etas <- sample(etas) ## suffle eta orders
cv.result3 <- cv.coxkl_ridge(z = ExampleData$z,
                             delta = ExampleData$status,
                             time = ExampleData$time,
                             stratum = NULL,
                             RS = NULL,
                             beta = ExampleData$beta_external,
                             etas = etas,
                             nfolds = 5, 
                             cv.criteria = "CIndex_pooled",
                             message = T)
head(cv.result3$integrated_stat.best_per_eta)
#           eta    lambda CIndex_pooled
# 1 0.000000000 104.53838     0.7086646
# 2 0.009953651  95.29411     0.7090724
# 3 0.020888146  95.34000     0.7098879
# 4 0.032900138  95.38928     0.7101427
# 5 0.046095806  95.44212     0.7108053
# 6 0.060591790  87.01481     0.7113660

cv.result3$external_stat
# 0.7388379

plot(cv.result3$integrated_stat.best_per_eta[,1], cv.result3$integrated_stat.best_per_eta[,3], type = "b", xlab = "eta")
abline(h = cv.result3$external_stat, col = "red", lty = 2)
#---------------------------6.coxkl_highdim------------------------------#
load("data/ExampleDataHighDim.RData")
result4 <- coxkl_highdim(z = ExampleData$z,
                         delta = ExampleData$status,
                         time = ExampleData$time,
                         stratum = NULL,
                         RS = NULL,
                         beta = ExampleData$beta_external,
                         eta = 0,
                         alpha = 1.0,
                         message = T)

library(grpreg)
X <- as.matrix(ExampleData$z)
y <- survival::Surv(ExampleData$time, ExampleData$status)
group <- 1:ncol(X)
fit_lasso <- grpreg::grpsurv(X, y, group = group, penalty = "grLasso", alpha = 1.0)

max(abs(fit_lasso$lambda - result4$lambda))
# [1] 1e-09

#---------------------------7. cv.coxkl_highdim ------------------------------#
load("data/ExampleDataHighDim.RData")
set.seed(1)
etas <- generate_eta(method = "exponential", n = 50, max_eta = 10)
etas <- sample(etas) 
cv.result5 <- cv.coxkl_highdim(z = ExampleData$z,
                               delta = ExampleData$status,
                               time = ExampleData$time,
                               stratum = NULL,
                               RS = NULL,
                               beta = ExampleData$beta_external,
                               etas = etas,
                               alpha = 1.0,
                               nfolds = 5, 
                               cv.criteria = "CIndex_pooled",
                               message = T)
head(cv.result5$integrated_stat.best_per_eta)
#           eta     lambda CIndex_pooled
# 1 0.000000000 0.03048043     0.7136595
# 2 0.009953651 0.03049407     0.7139653
# 3 0.020888146 0.02845257     0.7145770
# 4 0.032900138 0.02846728     0.7148318
# 5 0.046095806 0.02656338     0.7158002
# 6 0.060591790 0.02849992     0.7162080

cv.result5$external_stat
# 0.7388379


plot(cv.result5$integrated_stat.best_per_eta[,1], cv.result5$integrated_stat.best_per_eta[,3], type = "b", xlab = "eta")
abline(h = cv.result5$external_stat, col = "red", lty = 2)









