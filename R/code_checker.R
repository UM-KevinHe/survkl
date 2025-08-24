#--------------------------- step by step check ------------------------------#
load("data/ExampleData.RData")
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
dat <- data.frame(time = ExampleData$time,
                  status = ExampleData$status,
                  stratum = ExampleData$stratum,
                  ExampleData$z)

library(survival)
coxph_res <- coxph(
  Surv(time, status) ~ Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + strata(stratum),
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
                 etas = 0,
                 Mstop = 50,
                 backtrack = F,
                 message = F)
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



#-------------------- 3. check: cv.coxkl --------------------------# 
load("data/ExampleData.RData")

etas <- generate_eta(method = "exponential", n = 10, max_eta = 5)
cv_res1 <- cv.coxkl(z = ExampleData$z,
                    delta = ExampleData$status,
                    time = ExampleData$time,
                    stratum = NULL,
                    RS = NULL,
                    beta = ExampleData$beta_external,
                    etas = etas,
                    nfolds = 5,
                    criteria = c("V&VH"),
                    # criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
                    message = F)
cv_res1
#         eta  VVH_Loss
# 1  0.00000000 -2067.738
# 2  0.03374245 -2066.098
# 3  0.09002825 -2073.831
# 4  0.18391863 -2070.157
# 5  0.34053721 -2071.9793we2
# 6  0.60179276 -2066.785
# 7  1.03759328 -2065.022
# 8  1.76455236 -2064.752
# 9  2.97719318 -2064.582
# 10 5.00000000 -2066.030

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
#         138.1936    128.8798    120.1937    112.0931    104.5384
# [1,]  0.13888135  0.14450653  0.15017707  0.15588048  0.16160406
# [2,] -0.17547848 -0.18241304 -0.18944000 -0.19654779 -0.20372421
# [3,] -0.04192810 -0.04312603 -0.04429336 -0.04542704 -0.04652448
# [4,]  0.03213009  0.03176027  0.03127296  0.03066580  0.02993718
# [5,]  0.23628338  0.24513399  0.25408686  0.26312975  0.27225002
# [6,]  0.09007171  0.09307455  0.09608820  0.09910679  0.10212425


#---------------------------5. cv.coxkl_ridge ------------------------------#
cv.result3 <- cv.coxkl_ridge(z = ExampleData$z,
                             delta = ExampleData$status,
                             time = ExampleData$time,
                             stratum = NULL,
                             RS = NULL,
                             beta = ExampleData$beta_external,
                             etas = seq(0, 5, 1),
                             nfolds = 5, 
                             cv.eta.criteria = "V&VH", 
                             cv.lambda.criteria = "V&VH",
                             message = T)
cv.result3
#   eta  VVH_Loss
# 1   0 -1181.542
# 2   1 -1167.634
# 3   2 -1169.899
# 4   3 -1166.570
# 5   4 -1162.638
# 6   5 -1163.420

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

#---------------------------7. cv.coxkl_highdim ------------------------------#
result5 <- cv.coxkl_highdim(z = ExampleData$z,
                            delta = ExampleData$status,
                            time = ExampleData$time,
                            stratum = NULL,
                            RS = NULL,
                            beta = ExampleData$beta_external,
                            etas = seq(0, 5, 1),
                            alpha = 1.0,
                            nfolds = 5, 
                            cv.eta.criteria = "V&VH",
                            cv.lambda.criteria = "V&VH",
                            message = T)
# > result5
#   eta  VVH_Loss
# 1   0 -1172.989
# 2   1 -1167.563
# 3   2 -1164.809
# 4   3 -1164.275
# 5   4 -1164.608
# 6   5 -1164.935






