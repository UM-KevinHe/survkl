#' Example high-dimensional survival data
#'
#' A simulated survival dataset in a high-dimensional linear setting
#' with 50 covariates (6 signals + 44 AR(1) noise), Weibull baseline
#' hazard, and controlled censoring. Includes internal train/test sets,
#' and an external-data–estimated coefficient vector.
#'
#' @name ExampleData_highdim
#' @docType data
#' @usage data(ExampleData_highdim)
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{train}{A list with components:
#'     \describe{
#'       \item{z}{Data frame of size \eqn{n_\mathrm{train}\times 50} with covariates \code{Z1}–\code{Z50}.}
#'       \item{status}{Vector of event indicators (\code{1}=event, \code{0}=censored).}
#'       \item{time}{Numeric vector of observed times \eqn{\min(T, C)}.}
#'       \item{stratum}{Vector of stratum labels (here all \code{1}).}
#'     }
#'   }
#'   \item{test}{A list with the same structure as \code{train}, with size \eqn{n_\mathrm{test}\times 50} for \code{z}.}
#'   \item{beta_external}{Numeric vector (length 50, named \code{Z1}–\code{Z50}) of Cox coefficients
#'     estimated on an external dataset using only \code{Z1}–\code{Z6} and expanded to length 50
#'     (zeros for \code{Z7}–\code{Z50}).}
#' }
#'
#' @details Data-generating mechanism:
#' \itemize{
#'   \item Covariates: 50 variables with signals \code{Z1}–\code{Z6} and noise \code{Z7}–\code{Z50}.
#'     \itemize{
#'       \item \code{Z1}, \code{Z2} ~ bivariate normal with AR(1) correlation \eqn{\rho=0.5}.
#'       \item \code{Z3}, \code{Z4} ~ independent Bernoulli(0.5).
#'       \item \code{Z5} ~ \eqn{N(2,1)}, \code{Z6} ~ \eqn{N(-2,1)} (group indicator fixed at 1).
#'       \item \code{Z7}–\code{Z50} ~ multivariate normal with AR(1) correlation \eqn{\rho=0.5}.
#'     }
#'   \item True coefficients: \eqn{\beta = (0.3,-0.3,0.3,-0.3,0.3,-0.3,0,\ldots,0)} (length 50).
#'   \item Event times: Weibull baseline hazard
#'     \eqn{h_0(t)=\lambda\nu\, t^{\nu-1}} with \eqn{\lambda=1}, \eqn{\nu=2}.
#'     Given linear predictor \eqn{\eta = Z^\top \beta}, draw \eqn{U\sim\mathrm{Unif}(0,1)} and set
#'     \deqn{T = \left(\frac{-\log U}{\lambda\, e^{\eta}}\right)^{1/\nu}.}
#'   \item Censoring: \eqn{C\sim \mathrm{Unif}(0,\text{ub})} with \code{ub} tuned iteratively to
#'     achieve the target censoring rate (internal: \code{0.70}; external: \code{0.50}).
#'     Observed time is \eqn{\min(T,C)}, status is \eqn{\mathbf{1}\{T \le C\}}.
#'   \item External coefficients: Fit a Cox model
#'     \code{Surv(time, status) ~ Z1 + ... + Z6} on the external data (Breslow ties),
#'     then place the estimated coefficients into a length-50 vector (zeros elsewhere).
#' }
#'
#' @examples
#' data(ExampleData_highdim)
#' 
#' head(ExampleData_highdim$train$z)
#' table(ExampleData_highdim$train$status)
#' summary(ExampleData_highdim$train$time)
#'
#' head(ExampleData_highdim$test$z)
#' table(ExampleData_highdim$test$status)
#' summary(ExampleData_highdim$test$time)
#' 
"ExampleData_highdim"


#' Example low-dimensional survival data
#'
#' A simulated survival dataset in a low-dimensional linear setting
#' with 6 covariates (2 correlated continuous, 2 binary, 2 mean-shifted normals),
#' Weibull baseline hazard, and controlled censoring. Includes internal train/test sets,
#' and three external-quality coefficient vectors.
#'
#' @name ExampleData_lowdim
#' @docType data
#' @usage data(ExampleData_lowdim)
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{train}{A list with components:
#'     \describe{
#'       \item{z}{Data frame of size \eqn{n_\mathrm{train}\times 6} with covariates \code{Z1}–\code{Z6}.}
#'       \item{status}{Vector of event indicators (\code{1}=event, \code{0}=censored).}
#'       \item{time}{Numeric vector of observed times \eqn{\min(T, C)}.}
#'       \item{stratum}{Vector of stratum labels (here all \code{1}).}
#'     }
#'   }
#'   \item{test}{A list with the same structure as \code{train}, with size \eqn{n_\mathrm{test}\times 6} for \code{z}.}
#'   \item{beta_external_good}{Numeric vector (length 6; named \code{Z1}–\code{Z6}) of Cox coefficients estimated on a
#'     "Good" external dataset using all \code{Z1}–\code{Z6}.}
#'   \item{beta_external_fair}{Numeric vector (length 6; names \code{Z1}–\code{Z6}) of Cox coefficients estimated on a
#'     "Fair" external dataset using a reduced subset \code{Z1}, \code{Z3}, \code{Z5}, \code{Z6};
#'     coefficients for variables not used are \code{0}.}
#'   \item{beta_external_poor}{Numeric vector (length 6; names \code{Z1}–\code{Z6}) of Cox coefficients estimated on a
#'     "Poor" external dataset using \code{Z1} and \code{Z5} only; remaining entries are \code{0}.}
#' }
#'
#' @details Data-generating mechanism:
#' \itemize{
#'   \item Covariates: 6 variables \code{Z1}–\code{Z6}.
#'     \itemize{
#'       \item \code{Z1}, \code{Z2} ~ bivariate normal with AR(1) correlation \eqn{\rho=0.5}.
#'       \item \code{Z3}, \code{Z4} ~ independent Bernoulli(0.5).
#'       \item \code{Z5} ~ \eqn{N(2,1)}, \code{Z6} ~ \eqn{N(-2,1)} (group indicator fixed at 1 for internal train/test).
#'     }
#'   \item True coefficients: \eqn{\beta = (0.3,-0.3,0.3,-0.3,0.3,-0.3)} (length 6).
#'   \item Event times: Weibull baseline hazard
#'     \eqn{h_0(t)=\lambda\nu \, t^{\nu-1}} with \eqn{\lambda=1}, \eqn{\nu=2}.
#'     Given linear predictor \eqn{\eta = Z^\top \beta}, draw \eqn{U\sim\mathrm{Unif}(0,1)} and set
#'     \deqn{T = \left(\frac{-\log U}{\lambda \, e^{\eta}}\right)^{1/\nu}.}
#'   \item Censoring: \eqn{C\sim \mathrm{Unif}(0,\text{ub})} with \code{ub} tuned iteratively to
#'     achieve the target censoring rate (internal: \code{0.70}; external: \code{0.50}).
#'     Observed time is \eqn{\min(T,C)}, status is \eqn{\mathbf{1}\{T \le C\}}.
#'   \item External coefficients: For each quality level ("Good", "Fair", "Poor"), fit a Cox model
#'     \code{Surv(time, status) ~ Z1 + ...} on the corresponding external data (Breslow ties)
#'     using the specified covariate subset; place estimates into a length-6 vector named \code{Z1}–\code{Z6}
#'     with zeros for variables not included.
#' }
#'
#' @examples
#' data(ExampleData_lowdim)
#' 
#' head(ExampleData_lowdim$train$z)
#' table(ExampleData_lowdim$train$status)
#' summary(ExampleData_lowdim$train$time)
#'
#' head(ExampleData_lowdim$test$z)
#' table(ExampleData_lowdim$test$status)
#' summary(ExampleData_lowdim$test$time)
#' 
"ExampleData_lowdim"



#' Study to Understand Prognoses Preferences Outcomes and Risks of Treatment
#' @name support
#' @docType data
#' @usage data(support)
#' 
#' @description The \code{support} dataset tracks five response variables: hospital
#'   death, severe functional disability, hospital costs, and time until death
#'   and death itself. The patients are followed for up to 5.56 years. See Bhatnagar et al. (2020) for details.
#'
#' @details Some of the original data was missing. Before imputation, there were
#'   a total of 9,104 individuals and 47 variables. Following Bhatnagar et al. (2020), a few variables 
#'   were removed. Three response variables were removed:
#'   hospital charges, patient ratio of costs to charges and patient
#'   micro-costs. Hospital death was also removed as it was directly informative
#'   of the event of interest, namely death. Additionally, functional disability and
#'   income were removed as they are ordinal covariates. Finally, 8
#'   covariates were removed related to the results of previous findings: SUPPORT
#'   day 3 physiology score (\code{sps}), APACHE III day 3 physiology score
#'   (\code{aps}), SUPPORT model 2-month survival estimate, SUPPORT model
#'   6-month survival estimate, Physician's 2-month survival estimate for pt.,
#'   Physician's 6-month survival estimate for pt., Patient had Do Not
#'   Resuscitate (DNR) order, and Day of DNR order (<0 if before study). Of
#'   these, \code{sps} and \code{aps} were added on after imputation, as they
#'   were missing only 1 observation. First the imputation is done manually using the normal
#'   values for physiological measures recommended by Knaus et al. (1995). Next,
#'   a single dataset was imputed using \pkg{mice} with default settings. After
#'   imputation, the covariate for surrogate activities of daily
#'   living was not imputed. This is due to collinearity between the other two
#'   covariates for activities of daily living. Therefore, surrogate activities
#'   of daily living were removed. See details in the R package (casebase) by Bhatnagar et al. (2020).
#'
#' @format A data frame with 9,104 observations and 34 variables after imputation
#'   and the removal of response variables like hospital charges, patient ratio
#'   of costs to charges and micro-costs following Bhatnagar et al. (2020). 
#'   Ordinal variables, namely functional disability and income, were also removed. 
#'   Finally, Surrogate activities of daily living were removed due to sparsity. 
#'   There were 6 other model scores in the data-set and they were removed; only aps and sps were kept.
#'   \describe{ 
#'   \item{age}{ stores a double representing age. } 
#'   \item{death}{
#'   death at any time up to NDI (National Death Index) date: 12/31/1994. } 
#'   \item{sex}{ 0=female, 1=male. } 
#'   \item{slos}{ days from study entry to discharge. } 
#'   \item{d.time}{ days of
#'   follow-up. } 
#'   \item{dzgroup}{ each level of dzgroup: ARF/MOSF w/Sepsis,
#'   COPD, CHF, Cirrhosis, Coma, Colon Cancer, Lung Cancer, MOSF with
#'   malignancy. } 
#'   \item{dzclass}{ ARF/MOSF, COPD/CHF/Cirrhosis, Coma and cancer disease classes. } 
#'   \item{num.co}{ the number of comorbidities. }
#'   \item{edu}{ years of education of patients. } 
#'   \item{scoma}{ the SUPPORT coma score based on Glasgow D3. } 
#'   \item{avtisst}{ average TISS, days 3-25. }
#'   \item{race}{ indicates race: White, Black, Asian, Hispanic or other. }
#'   \item{hday}{ day in Hospital at Study Admit.} 
#'   \item{diabetes}{diabetes (Com27-28, Dx 73).} 
#'   \item{dementia}{dementia (Comorbidity 6).} 
#'   \item{ca}{cancer state. } 
#'   \item{meanbp}{ mean arterial blood pressure day 3. } 
#'   \item{wblc}{ white blood cell count on day 3. } 
#'   \item{hrt}{ heart rate day 3. }
#'   \item{resp}{ respiration rate day 3. } 
#'   \item{temp}{ temperature, in Celsius, on day 3. } 
#'   \item{pafi}{ PaO2/(0.01*FiO2) day 3. } 
#'   \item{alb}{serum albumin day 3. } 
#'   \item{bili}{ bilirubin day 3. } 
#'   \item{crea}{ serum creatinine day 3. } 
#'   \item{sod}{ serum sodium day 3. } 
#'   \item{ph}{ serum pH (in arteries) day 3. } 
#'   \item{glucose}{ serum glucose day 3. } 
#'   \item{bun}{ bun day 3. } 
#'   \item{urine}{ urine output day 3. } 
#'   \item{adlp}{ adl patient day 3. }  
#'   \item{adlsc}{ imputed adl calibrated to surrogate, if a surrogate was used for a follow up.} 
#'   \item{sps}{SUPPORT physiology score.}
#'   \item{aps}{apache III physiology score.} }
#'   
#' @source Available at the following website:
#'   \url{https://biostat.app.vumc.org/wiki/Main/SupportDesc}.
#' 
#' @references 
#' 
#' Bhatnagar, S., Turgeon, M., Islam, J., Hanley, J. A., and Saarela, O. (2020) casebase: Fitting Flexible Smooth-in-Time
#' Hazards and Risk Functions via Logistic and Multinomial Regression. 
#' \emph{R package version 0.9.0},
#' <https://CRAN.R-project.org/package=casebase>.
#' 
#' Knaus, W. A., Harrell, F. E., Lynn, J., Goldman, L., Phillips, R. S., Connors, A. F., et al. (1995) 
#' The SUPPORT prognostic model: Objective estimates of survival for seriously ill hospitalized adults. 
#' \emph{Annals of Internal Medicine}, \strong{122(3)}: 191-203.
#' \cr
#' 
#' 
#' @examples
#' \dontrun{
#' rm(list = ls())
#' library(devtools)
#' 
#' setwd("~/University of Michigan Dropbox/Lingfeng Luo/Lingfeng Research/R Tutorial Package/coxkl/coxkll/")
#' 
#' load_all()
#' 
#' data(support)
#' set.seed(123)
#' 
#' support <- support[support$ca %in% c("metastatic"),]
#' time <- support$d.time
#' death <- support$death
#' diabetes <-  model.matrix(~factor(support$diabetes))[,-1]
#' #sex: female as the reference group
#' sex <- model.matrix(~support$sex)[,-1]
#' #age: continuous variable
#' age <-support$age
#' age[support$age<=50] <- "<50"
#' age[support$age>50 & support$age<=60] <- "50-59"
#' age[support$age>60 & support$age<70] <- "60-69"
#' age[support$age>=70] <- "70+"
#' age <- factor(age, levels = c("60-69", "<50", "50-59", "70+"))
#' z_age <- model.matrix(~age)[,-1]
#' z <- data.frame(z_age, sex, diabetes)
#' colnames(z) <- c("age_50", "age_50_59", "age_70", "diabetes", "male")
#' data <- data.frame(time, death, z)
#' 
#' 
#' n <- nrow(data)
#' n_ext  <- floor(0.87 * n)
#' n_int  <- floor(0.03 * n)
#' n_test <- n - n_ext - n_int
#' idx     <- sample(seq_len(n))
#' idx_ext <- idx[       1:n_ext]
#' idx_int <- idx[(n_ext + 1):(n_ext + n_int)]
#' idx_test<- idx[(n_ext + n_int + 1):n]
#' 
#' external_data <- data[idx_ext, ]
#' internal_data <- data[idx_int, ]
#' test_data     <- data[idx_test, ]
#' 
#' library(survival)
#' ext_cox <- coxph(
#'   Surv(time, death) ~ age_50 + age_50_59 + age_70 + diabetes + male,
#'   data = external_data
#' )
#' beta_external <- coef(ext_cox)
#' 
#' result1 <- cv.coxkl(
#'   z        = internal_data[, c("age_50", "age_50_59", "age_70", "diabetes", "male")],
#'   delta    = internal_data$death,
#'   time     = internal_data$time,
#'   beta     = beta_external,
#'   eta_list = seq(0, 5, by = 1)
#' )
#' plot(result1$result)
#' #' do not run:
#' 
#' 
#' library(devtools)
#' rstudioapi::getActiveDocumentContext()$path
#' 
#' setwd("~/University of Michigan Dropbox/Lingfeng Luo/Lingfeng Research/R Tutorial Package/coxkl/coxkll/")
#' 
#' load_all()
#' 
#' data(support)
#' set.seed(123)
#' 
#' support <- support[support$ca %in% c("metastatic"),]
#' time <- support$d.time
#' death <- support$death
#' diabetes <-  model.matrix(~factor(support$diabetes))[,-1]
#' #sex: female as the reference group
#' sex <- model.matrix(~support$sex)[,-1]
#' #age: continuous variable
#' age <-support$age
#' age[support$age<=50] <- "<50"
#' age[support$age>50 & support$age<=60] <- "50-59"
#' age[support$age>60 & support$age<70] <- "60-69"
#' age[support$age>=70] <- "70+"
#' age <- factor(age, levels = c("60-69", "<50", "50-59", "70+"))
#' z_age <- model.matrix(~age)[,-1]
#' z <- data.frame(z_age, sex, diabetes)
#' colnames(z) <- c("age_50", "age_50_59", "age_70", "diabetes", "male")
#' data <- data.frame(time, death, z)
#' 
#' 
#' n <- nrow(data)
#' n_ext  <- floor(0.87 * n)
#' n_int  <- floor(0.03 * n)
#' n_test <- n - n_ext - n_int
#' idx     <- sample(seq_len(n))
#' idx_ext <- idx[       1:n_ext]
#' idx_int <- idx[(n_ext + 1):(n_ext + n_int)]
#' idx_test<- idx[(n_ext + n_int + 1):n]
#' 
#' external_data <- data[idx_ext, ]
#' internal_data <- data[idx_int, ]
#' test_data     <- data[idx_test, ]
#' 
#' library(survival)
#' ext_cox <- coxph(
#'   Surv(time, death) ~ age_50 + age_50_59 + age_70 + diabetes + male,
#'   data = external_data
#' )
#' beta_external <- coef(ext_cox)
#' 
#' result1 <- cv.coxkl(
#'   z        = internal_data[, c("age_50", "age_50_59", "age_70", "diabetes", "male")],
#'   delta    = internal_data$death,
#'   time     = internal_data$time,
#'   beta     = beta_external,
#'   eta_list = seq(0, 5, by = 1)
#' )
#' plot(result1$result)
#' }
"support"