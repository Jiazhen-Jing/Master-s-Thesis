library(semPlot)
library(lavaan)
library(blavaan)
library(sn)  
library(moments)  
library(simsem)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(semTools)
library(MASS)
library(gridExtra)
# Population ####
Lambda <- matrix(c(0.4, 0, 0, 
                   0.5, 0, 0,
                   0.6, 0, 0, 
                   0.7, 0, 0, 
                   0, 0.9, 0, 
                   0, 0.8, 0, 
                   0, 0.7, 0,  
                   0, 0.6, 0, 
                   0, 0, 0.5,  
                   0, 0, 0.7, 
                   0, 0, 0.6,  
                   0, 0, 0.7), 12, 3, byrow = T)


Phi <- matrix(c(1, 0.3, 0.3, 
                0.3, 1, 0.3, 
                0.3, 0.3, 1), 3, 3, byrow = T)

delta <- diag(12) - Lambda%*%Phi%*%t(Lambda)
diagT <- diag(delta)
Theta <- diag(diagT)
Sigma <- Theta + Lambda%*%Phi%*%t(Lambda)
truepara <- c(0.4, 0.5, 0.6, 0.7, 
              0.9, 0.8, 0.7, 0.6, 
              0.5, 0.7, 0.6, 0.7,
              0.3, 0.3, 0.3)

# Analysis model ####
analysisModel <- "
f1 =~ y1 + y2 + y3 + y4
f2 =~ y5 + y6 + y7 + y8
f3 =~ y9 + y10 + y11 + y12
"

# Functions ####
## Ordinal transfer ####
transform_to_ordered <- function(x, k, probs) {
  quantiles <- quantile(x, probs)
  ordered_data <- as.ordered(cut(x, breaks = c(-Inf, quantiles, Inf),
                                 labels = 1:k, include.lowest = TRUE))
  return(ordered_data)
}
## Convergence check ####
converge  <- function(Fit) {
  converged_results <- list()
  for (i in 1:length(Fit)) {
    con <- lavInspect(Fit[[i]], "converged")
    if (con == "TRUE") {
      converged_results[[i]] <- Fit[[i]]
    }else {
      converged_results[[i]] <- NULL
    }
  }
  Con <- converged_results[!sapply(converged_results, is.null)]
  return(Con)
}

## Admissible check ####
admissible_check <- function(Fit) {
  ad_results <- list()
  for (i in 1:length(Fit)) {
    result <- lavInspect(Fit[[i]], "est")
    if (sum(diag(result$theta) < 0) > 0|| sum(result$psi > 1) > 0 || sum(result$psi < -1) > 0) {
      ad_results[[i]] <- NULL  
    }else{
      ad_results[[i]] <- Fit[[i]]
    }
  }
  Ad <- ad_results[!sapply(ad_results, is.null)]
  return(Ad)
}

## Goodness of Fit values extract ####
gof_check <- function(Fit){
  gof <- matrix(0.0, nrow = length(Fit), ncol = 5)
  for (i in 1:length(Fit)) {
    measures <- fitmeasures(Fit[[i]])  
    gof[i, 1] <- measures["pvalue.scaled"]
    gof[i, 2] <- measures["rmsea"]
    gof[i, 3] <- measures["cfi.scaled"]
    gof[i, 4] <- measures["tli.scaled"]
    gof[i, 5] <- measures["srmr_bentler"]
    print(paste("Rep:", i))
  }
  colnames(gof) <- c("P_Value", "RMSEA", "CFI", "TLI", "SRMR")
  return(gof)
}

gof_check_bayes <- function(Fit){
  gof <- matrix(0.0, nrow = length(Fit), ncol = 4)
  for (i in 1:length(Fit)) {
    measures <- fitmeasures(Fit[[i]])  # Call fitmeasures once and store the result
    br <-  blavFitIndices(Fit[[i]])
    gof[i, 1] <- measures["ppp"]
    gof[i, 2] <- measures["bic"]
    gof[i, 3] <- measures["p_loo"]
    gof[i, 4] <- mean(br@indices$BRMSEA)
    print(paste("Rep:", i))
  }
  colnames(gof) <- c("PPP", "BIC", "P_LOO", "BRMSEA")
  return(gof)
}


## Extract estimates and se per Rep####
Estimate <- function(Fit){
  esti <- matrix(0.0, nrow = 15, ncol = 2)
  aa <- lavInspect(Fit, what = "std.all")
  esti[, 1] <- c(aa$lambda[1:4, 1], aa$lambda[5:8, 2], aa$lambda[9:12, 3], aa$psi[lower.tri(aa$psi, diag = FALSE)])
  bb <- lavInspect(Fit, what = "se")
  esti[, 2] <- c(bb$lambda[1:4, 1], bb$lambda[5:8, 2], bb$lambda[9:12, 3], bb$psi[lower.tri(aa$psi, diag = FALSE)])
  
  rownames(esti) <- c(paste0("lambda", 1:12), paste0("phi", 1:3))
  colnames(esti) <- c("std.est", "se")
  return(esti)
}

## ARB ARMSE calculation for parameter estimates####
AR.para <-function(fitt, true){
  rmse.lam <- numeric(length(fitt))
  rmse.phi <- numeric(length(fitt))
  rmse <- numeric(length(fitt))
  rb.lam <- numeric(length(fitt))
  rb.phi <- numeric(length(fitt))
  rb <- numeric(length(fitt))
  
  for (i in 1:length(fitt)) {
    para <- Estimate(fitt[[i]])
    estall <- para[, 1]
    all <- (estall - true)/true
    lam <- all[1:12]
    phi <- all[13:15]
    rb.lam[i] <- mean(lam)
    rb.phi[i] <- mean(phi)
    rb[i] <- mean(all)
    rmse.lam[i] <- sqrt(mean((lam^2)))
    rmse.phi[i] <- sqrt(mean((phi^2)))
    rmse[i] <- sqrt(mean((all^2)))
  }
  arb.lam <- mean(rb.lam)
  arb.phi <- mean(rb.phi)
  arb <- mean(rb)
  
  armse.lam <- mean(rmse.lam)
  armse.phi <- mean(rmse.phi)
  armse <- mean(rmse)
  
  return(list(ARB.lam = arb.lam,
              ARB.phi = arb.phi,
              ARB = arb,
              ARMSE.lam = armse.lam,
              ARMSE.phi = armse.phi,
              ARMSE = armse))
}

# Data generation ####
# dat 1 means sample size (1 = 100, 2=400, 3 = 1000)
# dat 2 means categories (1 = 3, 2 = 5, 3 = 7)
# dat 3 means sy or asy (1 = sy, 2 = asy)
# digits 4 means method (1=RML, 2=DWLS, 3=BNI, 4=BI)
set.seed(123)
## dat1 n = 100 ####
Rep <-500
dat1<- list()
for (i in 1:Rep) {
  simulated_data <-  mvrnorm(n = 100, mu = rep(0, 12), Sigma = Sigma)
  colnames(simulated_data) <- c("y1", "y2", "y3", "y4", 
                                "y5", "y6", "y7", "y8", 
                                "y9", "y10", "y11", "y12")
  dat1[[i]] <- data.frame(simulated_data)
}

### dat111 n =100, cate = 3, sy####
dat111 <- list()
for (i in 1:length(dat1)) {
  transformed_data <- apply(dat1[[i]], 2, transform_to_ordered, k = 3, probs = c(0.2, 0.8)) 
  dat111[[i]] <- transformed_data
  print(paste("Rep:", i))
}

### dat112 n =100, cate = 3, asy####
dat112 <- list()
for (i in 1:length(dat1)) {
  transformed_data <- apply(dat1[[i]], 2, transform_to_ordered, k = 3, probs = c(0.3, 0.9)) 
  dat112[[i]] <- transformed_data
}

save(dat1, dat111, dat112, file = "dat11.RData")

### dat121 n =100, cate = 5, sy####
dat121 <- list()
for (i in 1:length(dat1)) {
  transformed_data <- apply(dat1[[i]], 2, transform_to_ordered, k = 5, probs = c(0.1, 0.3, 0.7, 0.9)) 
  dat121[[i]] <- transformed_data
  print(paste("Rep:", i))
}
### dat122 n =100, cate = 5, asy####
dat122 <- list()
for (i in 1:length(dat1)) {
  transformed_data <- apply(dat1[[i]], 2, transform_to_ordered, k = 5, probs = c(0.25, 0.7, 0.85, 0.95)) 
  dat122[[i]] <- transformed_data
  print(paste("Rep:", i))
}
save(dat121, dat122, file = "dat12.RData")


### dat131 n =100, cate = 7, sy####
dat131 <- list()
for (i in 1:length(dat1)) {
  transformed_data <- apply(dat1[[i]], 2, transform_to_ordered, k = 7, probs = c(0.05, 0.17, 0.35, 0.65, 0.83, 0.95)) 
  dat131[[i]] <- transformed_data
  print(paste("Rep:", i))
}
### dat132 n =100, cate = 7, asy####
dat132 <- list()
for (i in 1:length(dat1)) {
  transformed_data <- apply(dat1[[i]], 2, transform_to_ordered, k = 7, probs = c(0.05, 0.2, 0.5, 0.77, 0.9, 0.97)) 
  dat132[[i]] <- transformed_data
  print(paste("Rep:", i))
}
save(dat131, dat132, file = "dat13.RData")


## dat2 n = 400 ####
Rep <-500
dat2<- list()
for (i in 1:Rep) {
  simulated_data <-  mvrnorm(n = 400, mu = rep(0, 12), Sigma = Sigma)
  colnames(simulated_data) <- c("y1", "y2", "y3", "y4", 
                                "y5", "y6", "y7", "y8", 
                                "y9", "y10", "y11", "y12")
  dat2[[i]] <- data.frame(simulated_data)
}



### dat211 n =400, cate = 3, sy####
dat211 <- list()
for (i in 1:length(dat2)) {
  transformed_data <- apply(dat2[[i]], 2, transform_to_ordered, k = 3, probs = c(0.2, 0.8)) 
  dat211[[i]] <- transformed_data
  print(paste("Rep:", i))
}
### dat212 n =400, cate = 3, asy####
dat212 <- list()
for (i in 1:length(dat2)) {
  transformed_data <- apply(dat2[[i]], 2, transform_to_ordered, k = 3, probs = c(0.3, 0.9)) 
  dat212[[i]] <- transformed_data
  print(paste("Rep:", i))
}
save(dat2, dat211, dat212, file = "dat21.RData")

### dat221 n =400, cate = 5, sy####
dat221 <- list()
for (i in 1:length(dat2)) {
  transformed_data <- apply(dat2[[i]], 2, transform_to_ordered, k = 5, probs = c(0.1, 0.3, 0.7, 0.9)) 
  dat221[[i]] <- transformed_data
}

### dat222 n =400, cate = 5, asy####
dat222 <- list()
for (i in 1:length(dat2)) {
  transformed_data <- apply(dat2[[i]], 2, transform_to_ordered, k = 5, probs = c(0.25, 0.7, 0.85, 0.95)) 
  dat222[[i]] <- transformed_data
}
save(dat2, dat221, dat222, file = "dat22.RData")

### dat231 n =400, cate = 7, sy####
dat231 <- list()
for (i in 1:length(dat2)) {
  transformed_data <- apply(dat2[[i]], 2, transform_to_ordered, k = 7, probs = c(0.05, 0.17, 0.35, 0.65, 0.83, 0.95)) 
  dat231[[i]] <- transformed_data
  print(paste("Rep:", i))
}
### dat232 n =400, cate = 7, asy####
dat232 <- list()
for (i in 1:length(dat2)) {
  transformed_data <- apply(dat2[[i]], 2, transform_to_ordered, k = 7, probs = c(0.05, 0.2, 0.5, 0.77, 0.9, 0.97)) 
  dat232[[i]] <- transformed_data
  print(paste("Rep:", i))
}
save(dat231, dat232, file = "dat23.RData")



## dat3 n = 1000 ####
Rep <-500
dat3<- list()
for (i in 1:Rep) {
  simulated_data <-  mvrnorm(n = 1000, mu = rep(0, 12), Sigma = Sigma)
  colnames(simulated_data) <- c("y1", "y2", "y3", "y4", 
                                "y5", "y6", "y7", "y8", 
                                "y9", "y10", "y11", "y12")
  dat3[[i]] <- data.frame(simulated_data)
}



### dat311 n =1000, cate = 3, sy####
dat311 <- list()
for (i in 1:length(dat3)) {
  transformed_data <- apply(dat3[[i]], 2, transform_to_ordered, k = 3, probs = c(0.2, 0.8)) 
  dat311[[i]] <- transformed_data
  print(paste("Rep:", i))
}
### dat312 n =1000, cate = 3, asy####
dat312 <- list()
for (i in 1:length(dat3)) {
  transformed_data <- apply(dat3[[i]], 2, transform_to_ordered, k = 3, probs = c(0.3, 0.9)) 
  dat312[[i]] <- transformed_data
  print(paste("Rep:", i))
}
save(dat3, dat311, dat312, file = "dat31.RData")

### dat321 n =1000, cate = 5, sy####
dat321 <- list()
for (i in 1:length(dat3)) {
  transformed_data <- apply(dat3[[i]], 2, transform_to_ordered, k = 5, probs = c(0.1, 0.3, 0.7, 0.9)) 
  dat321[[i]] <- transformed_data
}

### dat322 n =1000, cate = 5, asy####

dat322 <- list()
for (i in 1:length(dat3)) {
  transformed_data <- apply(dat3[[i]], 2, transform_to_ordered, k = 5, probs = c(0.25, 0.7, 0.85, 0.95)) 
  dat322[[i]] <- transformed_data
}
save(dat321, dat322, file = "dat32.RData")



### dat331 n =1000, cate = 7, sy####
dat331 <- list()
for (i in 1:length(dat3)) {
  transformed_data <- apply(dat3[[i]], 2, transform_to_ordered, k = 7, probs = c(0.05, 0.17, 0.35, 0.65, 0.83, 0.95)) 
  dat331[[i]] <- transformed_data
  print(paste("Rep:", i))
}
### dat332 n =1000, cate = 7, asy####
dat332 <- list()
for (i in 1:length(dat3)) {
  transformed_data <- apply(dat3[[i]], 2, transform_to_ordered, k = 7, probs = c(0.05, 0.2, 0.5, 0.77, 0.9, 0.97)) 
  dat332[[i]] <- transformed_data
  print(paste("Rep:", i))
}
save(dat331, dat332, file = "dat33.RData")




# dat111 dat112 ####
## Case 1 fit1111 ####
## n=100, cate=3, sy, RML
fit1111 <- list()
for (i in 1:length(dat111)) {
  dat <- dat111[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit1111[[i]] <- result
  print(paste("Rep:", i))
}

fit1111.con <- converge(fit1111)
(500 - length(fit1111.con))/500
fit1111.con.ad <- admissible_check(fit1111.con)
(length(fit1111.con) - length(fit1111.con.ad))/length(fit1111.con)
gof1111 <- gof_check(fit1111.con.ad)
sum(gof1111[,1]<0.05)/length(fit1111.con.ad)
AR.p1111 <- AR.para(fit1111.con.ad, truepara)
save(fit1111, fit1111.con.ad, gof1111, AR.p1111, file = "f1111.RData")
colMeans(gof1111)

## Case 2 fit1112 ####
## n=100, cate=3, sy, DWLS
fit1112 <- list()
for (i in 1:length(dat111)) {
  dat <- dat111[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit1112[[i]] <- result
  print(paste("Rep:", i))
}
fit1112.con <- converge(fit1112)
(500 - length(fit1112.con))/500
fit1112.con.ad <- admissible_check(fit1112.con)
(length(fit1112.con) - length(fit1112.con.ad))/length(fit1112.con)
gof1112 <- gof_check(fit1112.con.ad)
sum(gof1112[,1]<0.05)/length(fit1112.con.ad)
AR.p1112 <- AR.para(fit1112.con.ad, truepara)
save(fit1112, fit1112.con.ad, gof1112, AR.p1112, file = "f1112.RData")
colMeans(gof1112)

## Case 3 fit1113 ####
## n=100, cate=3, sy, Bayes NI
fit1113 <- list()
for (i in 1:length(dat111)) {
  dat <- dat111[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(2, 0.5)[prec]"),
                 n.chains = 2)
  fit1113[[i]] <- result
  print(paste("Rep:", i))
}
fit1113.con <- converge(fit1113)
(500 - length(fit1113.con))/500
fit1113.con.ad <- admissible_check(fit1113.con)
(length(fit1113.con) - length(fit1113.con.ad))/length(fit1113.con)
gof1113 <- gof_check_bayes(fit1113.con.ad)
sum(gof1113[,1]<0.05)/length(fit1113.con.ad)
AR.p1113 <- AR.para(fit1113.con.ad, truepara)
save(fit1113, fit1113.con.ad, gof1113, AR.p1113, file = "f1113.RData")

## Case 4 fit1114 ####
## n=100, cate=3, sy, Bayes I
fit1114 <- list()
for (i in 1:length(dat111)) {
  dat <- dat111[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit1114[[i]] <- result
  print(paste("Rep:", i))
}
fit1114.con <- converge(fit1114)
(500 - length(fit1114.con))/500
fit1114.con.ad <- admissible_check(fit1114.con)
(length(fit1114.con) - length(fit1114.con.ad))/length(fit1114.con)
gof1114 <- gof_check_bayes(fit1114.con.ad)
sum(gof1114[,1]<0.05)/length(fit1114.con.ad)
AR.p1114 <- AR.para(fit1114.con.ad, truepara)
save(fit1114, fit1114.con.ad, gof1114, AR.p1114, file = "f1114.RData")
colMeans(gof1114)


## Case 5 fit1121 ####
## n=100, cate=3, asy, RML
fit1121 <- list()
for (i in 1:length(dat112)) {
  dat <- dat112[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit1121[[i]] <- result
  print(paste("Rep:", i))
}

fit1121.con <- converge(fit1121)
(500 - length(fit1121.con))/500
fit1121.con.ad <- admissible_check(fit1121.con)
(length(fit1121.con) - length(fit1121.con.ad))/length(fit1121.con)
gof1121 <- gof_check(fit1121.con.ad)
sum(gof1121[,1]<0.05)/length(fit1121.con.ad)
AR.p1121 <- AR.para(fit1121.con.ad, truepara)
save(fit1121, fit1121.con.ad, gof1121, AR.p1121, file = "f1121.RData")
colMeans(gof1121)

## Case 6 fit1122 ####
## n=100, cate=3, asy, DWLS
fit1122 <- list()
for (i in 1:length(dat112)) {
  dat <- dat112[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit1122[[i]] <- result
  print(paste("Rep:", i))
}
fit1122.con <- converge(fit1122)
(500 - length(fit1122.con))/500
fit1122.con.ad <- admissible_check(fit1122.con)
(length(fit1122.con) - length(fit1122.con.ad))/length(fit1122.con)
gof1122 <- gof_check(fit1122.con.ad)
sum(gof1122[,1]<0.05)/length(fit1122.con.ad)
AR.p1122 <- AR.para(fit1122.con.ad, truepara)
save(fit1122, fit1122.con.ad, gof1122, AR.p1122, file = "f1122.RData")
colMeans(gof1122)

## Case 7 fit1123 ####
## n=100, cate=3, asy, Bayes NI
fit1123 <- list()
for (i in 1:length(dat112)) {
  dat <- dat112[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(1, 0.5)[prec]"),
                 n.chains = 2)
  fit1123[[i]] <- result
  print(paste("Rep:", i))
}
fit1123.con <- converge(fit1123)
(500 - length(fit1123.con))/500
fit1123.con.ad <- admissible_check(fit1123.con)
(length(fit1123.con) - length(fit1123.con.ad))/length(fit1123.con)
gof1123 <- gof_check_bayes(fit1123.con.ad)
sum(gof1123[,1]<0.05)/length(fit1123.con.ad)
AR.p1123 <- AR.para(fit1123.con.ad, truepara)
save(fit1123, fit1123.con.ad, gof1123, AR.p1123, file = "f1123.RData")





## Case 8 fit1124 ####
## n=100, cate=3, asy, Bayes I
fit1124 <- list()
for (i in 1:length(dat112)) {
  dat <- dat112[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit1124[[i]] <- result
  print(paste("Rep:", i))
}
fit1124.con <- converge(fit1124)
(500 - length(fit1124.con))/500
fit1124.con.ad <- admissible_check(fit1124.con)
(length(fit1124.con) - length(fit1124.con.ad))/length(fit1124.con)
gof1124 <- gof_check_bayes(fit1124.con.ad)
sum(gof1124[,1]<0.05)/length(fit1124.con.ad)
AR.p1124 <- AR.para(fit1124.con.ad, truepara)
save(fit1124, fit1124.con.ad, gof1124, AR.p1124, file = "f1124.RData")
colMeans(gof1124)



# dat121 dat122 ####
## Case 1 fit1211 ####
## n=100, cate=5, sy, RML
fit1211 <- list()
for (i in 1:length(dat121)) {
  dat <- dat121[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit1211[[i]] <- result
  print(paste("Rep:", i))
}

fit1211.con <- converge(fit1211)
(500-length(fit1211.con))/500

fit1211.con.ad <- admissible_check(fit1211.con)
(length(fit1211.con)- length(fit1211.con.ad))/length(fit1211.con)
gof1211 <- gof_check(fit1211.con.ad)
sum(gof1211[,1]<0.05)/length(fit1211.con.ad)
AR.p1211 <- AR.para(fit1211.con.ad, truepara)
save(fit1211, fit1211.con.ad, gof1211, AR.p1211, file = "f1211.RData")
colMeans(gof1211)

## Case 2 fit1212 ####
## n=100, cate=5, sy, DWLS
fit1212 <- list()
for (i in 1:length(dat121)) {
  dat <- dat121[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit1212[[i]] <- result
  print(paste("Rep:", i))
}

fit1212.con <- converge(fit1212)
(500 - length(fit1212.con))/500
fit1212.con.ad <- admissible_check(fit1212.con)
(length(fit1212.con) - length(fit1212.con.ad))/length(fit1212.con)
gof1212 <- gof_check(fit1212.con.ad)
sum(gof1212[,1]<0.05)/length(fit1212.con.ad)
AR.p1212 <- AR.para(fit1212.con.ad, truepara)
save(fit1212, fit1212.con.ad, gof1212, AR.p1212, file = "f1212.RData")
colMeans(gof1212)

## Case 3 fit1213 ####
## n=100, cate=5, sy, Bayes NI
fit1213 <- list()
for (i in 1:length(dat121)) {
  dat <- dat121[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(1, 0.5)[prec]"),
                 n.chains = 2)
  fit1213[[i]] <- result
  print(paste("Rep:", i))
}
fit1213.con <- converge(fit1213)
(500 - length(fit1213.con))/500
fit1213.con.ad <- admissible_check(fit1213.con)
(length(fit1213.con) - length(fit1213.con.ad))/length(fit1213.con)
gof1213 <- gof_check_bayes(fit1213.con.ad)
sum(gof1213[,1]<0.05)/length(fit1213.con.ad)
AR.p1213 <- AR.para(fit1213.con.ad, truepara)
save(fit1213, fit1213.con.ad, gof1213, AR.p1213, file = "f1213.RData")

## Case 4 fit1214 ####
## n=100, cate=5, sy, Bayes I
fit1214 <- list()
for (i in 1:length(dat121)) {
  dat <- dat121[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit1214[[i]] <- result
  print(paste("Rep:", i))
}
fit1214.con <- converge(fit1214)
(500 - length(fit1214.con))/500
fit1214.con.ad <- admissible_check(fit1214.con)
(length(fit1214.con) - length(fit1214.con.ad))/length(fit1214.con)
gof1214 <- gof_check_bayes(fit1214.con.ad)
sum(gof1214[,1]<0.05)/length(fit1214.con.ad)
AR.p1214 <- AR.para(fit1214.con.ad, truepara)
save(fit1214, fit1214.con.ad, gof1214, AR.p1214, file = "f1214.RData")
colMeans(gof1214)



## Case 5 fit1221 ####
## n=100, cate=5, asy, RML
fit1221 <- list()
for (i in 1:length(dat122)) {
  dat <- dat122[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit1221[[i]] <- result
  print(paste("Rep:", i))
}

fit1221.con <- converge(fit1221)
(500-length(fit1221.con))/500

fit1221.con.ad <- admissible_check(fit1221.con)
(length(fit1221.con)- length(fit1221.con.ad))/length(fit1221.con)
gof1221 <- gof_check(fit1221.con.ad)
sum(gof1221[,1]<0.05)/length(fit1221.con.ad)
AR.p1221 <- AR.para(fit1221.con.ad, truepara)
save(fit1221, fit1221.con.ad, gof1221, AR.p1221, file = "f1221.RData")
colMeans(gof1221)

## Case 6 fit1222 ####
## n=100, cate=5, asy, DWLS
fit1222 <- list()
for (i in 1:length(dat122)) {
  dat <- dat122[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit1222[[i]] <- result
  print(paste("Rep:", i))
}

fit1222.con <- converge(fit1222)
(500 - length(fit1222.con))/500
fit1222.con.ad <- admissible_check(fit1222.con)
(length(fit1222.con) - length(fit1222.con.ad))/length(fit1222.con)
gof1222 <- gof_check(fit1222.con.ad)
sum(gof1222[,1]<0.05)/length(fit1222.con.ad)
AR.p1222 <- AR.para(fit1222.con.ad, truepara)
save(fit1222, fit1222.con.ad, gof1222, AR.p1222, file = "f1222.RData")
colMeans(gof1222)


## Case 7 fit1223 ####
## n=100, cate=5, asy, Bayes NI
fit1223 <- list()
for (i in 1:length(dat122)) {
  dat <- dat122[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(1, 0.5)[prec]"),
                 n.chains = 2)
  fit1223[[i]] <- result
  print(paste("Rep:", i))
}
fit1223.con <- converge(fit1223)
(500 - length(fit1223.con))/500
fit1223.con.ad <- admissible_check(fit1223.con)
(length(fit1223.con) - length(fit1223.con.ad))/length(fit1223.con)
gof1223 <- gof_check_bayes(fit1223.con.ad)
sum(gof1223[,1]<0.05)/length(fit1223.con.ad)
AR.p1223 <- AR.para(fit1223.con.ad, truepara)
save(fit1223, fit1223.con.ad, gof1223, AR.p1223, file = "f1223.RData")


## Case 8 fit1224 ####
## n=100, cate=5, asy, Bayes I
fit1224 <- list()
for (i in 1:length(dat122)) {
  dat <- dat122[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit1224[[i]] <- result
  print(paste("Rep:", i))
}
fit1224.con <- converge(fit1224)
(500 - length(fit1224.con))/500
fit1224.con.ad <- admissible_check(fit1224.con)
(length(fit1224.con) - length(fit1224.con.ad))/length(fit1224.con)
gof1224 <- gof_check_bayes(fit1224.con.ad)
sum(gof1224[,1]<0.05)/length(fit1224.con.ad)
AR.p1224 <- AR.para(fit1224.con.ad, truepara)
save(fit1224, fit1224.con.ad, gof1224, AR.p1224, file = "f1224.RData")
colMeans(gof1224)



# dat131 dat132 ####
## Case 1 fit1311 ####
## n=100, cate=7, sy, RML
fit1311 <- list()
for (i in 1:length(dat131)) {
  dat <- dat131[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit1311[[i]] <- result
  print(paste("Rep:", i))
}

fit1311.con <- converge(fit1311)
(500-length(fit1311.con))/500

fit1311.con.ad <- admissible_check(fit1311.con)
(length(fit1311.con)- length(fit1311.con.ad))/length(fit1311.con)
gof1311 <- gof_check(fit1311.con.ad)
sum(gof1311[,1]<0.05)/length(fit1311.con.ad)
AR.p1311 <- AR.para(fit1311.con.ad, truepara)
save(fit1311, fit1311.con.ad, gof1311, AR.p1311, file = "f1311.RData")
colMeans(gof1311)

## Case 2 fit1312 ####
## n=100, cate=7, sy, DWLS
fit1312 <- list()
for (i in 1:length(dat131)) {
  dat <- dat131[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit1312[[i]] <- result
  print(paste("Rep:", i))
}

fit1312.con <- converge(fit1312)
(500 - length(fit1312.con))/500
fit1312.con.ad <- admissible_check(fit1312.con)
(length(fit1312.con) - length(fit1312.con.ad))/length(fit1312.con)
gof1312 <- gof_check(fit1312.con.ad)
sum(gof1312[,1]<0.05)/length(fit1312.con.ad)
AR.p1312 <- AR.para(fit1312.con.ad, truepara)
save(fit1312, fit1312.con.ad, gof1312, AR.p1312, file = "f1312.RData")
colMeans(gof1312)

## Case 3 fit1313 ####
## n=100, cate=7, sy, Bayes NI
fit1313 <- list()
for (i in 1:length(dat131)) {
  dat <- dat131[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(1, 0.5)[prec]"),
                 n.chains = 2)
  fit1313[[i]] <- result
  print(paste("Rep:", i))
}
fit1313.con <- converge(fit1313)
(500 - length(fit1313.con))/500
fit1313.con.ad <- admissible_check(fit1313.con)
(length(fit1313.con) - length(fit1313.con.ad))/length(fit1313.con)
gof1313 <- gof_check_bayes(fit1313.con.ad)
sum(gof1313[,1]<0.05)/length(fit1313.con.ad)
AR.p1313 <- AR.para(fit1313.con.ad, truepara)
save(fit1313, fit1313.con.ad, gof1313, AR.p1313, file = "f1313.RData")

## Case 4 fit1314 ####
## n=100, cate=7, sy, Bayes I
fit1314 <- list()
for (i in 1:length(dat131)) {
  dat <- dat131[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit1314[[i]] <- result
  print(paste("Rep:", i))
}
fit1314.con <- converge(fit1314)
(500 - length(fit1314.con))/500
fit1314.con.ad <- admissible_check(fit1314.con)
(length(fit1314.con) - length(fit1314.con.ad))/length(fit1314.con)
gof1314 <- gof_check_bayes(fit1314.con.ad)
sum(gof1314[,1]<0.05)/length(fit1314.con.ad)
AR.p1314 <- AR.para(fit1314.con.ad, truepara)
save(fit1314, fit1314.con.ad, gof1314, AR.p1314, file = "f1314.RData")
colMeans(gof1314)



## Case 5 fit1321 ####
## n=100, cate=7, asy, RML
fit1321 <- list()
for (i in 1:length(dat132)) {
  dat <- dat132[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit1321[[i]] <- result
  print(paste("Rep:", i))
}

fit1321.con <- converge(fit1321)
(500-length(fit1321.con))/500

fit1321.con.ad <- admissible_check(fit1321.con)
(length(fit1321.con)- length(fit1321.con.ad))/length(fit1321.con)
gof1321 <- gof_check(fit1321.con.ad)
sum(gof1321[,1]<0.05)/length(fit1321.con.ad)
AR.p1321 <- AR.para(fit1321.con.ad, truepara)
save(fit1321, fit1321.con.ad, gof1321, AR.p1321, file = "f1321.RData")
colMeans(gof1321)



## Case 6 fit1322 ####
## n=100, cate=7, asy, DWLS
fit1322 <- list()
for (i in 1:length(dat132)) {
  dat <- dat132[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit1322[[i]] <- result
  print(paste("Rep:", i))
}

fit1322.con <- converge(fit1322)
(500 - length(fit1322.con))/500
fit1322.con.ad <- admissible_check(fit1322.con)
(length(fit1322.con) - length(fit1322.con.ad))/length(fit1322.con)
gof1322 <- gof_check(fit1322.con.ad)
sum(gof1322[,1]<0.05)/length(fit1322.con.ad)
AR.p1322 <- AR.para(fit1322.con.ad, truepara)
save(fit1322, fit1322.con.ad, gof1322, AR.p1322, file = "f1322.RData")
colMeans(gof1322)


## Case 7 fit1323 ####
## n=100, cate=7, asy, Bayes NI
fit1323 <- list()
for (i in 1:length(dat132)) {
  dat <- dat132[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit1323[[i]] <- result
  print(paste("Rep:", i))
}
fit1323.con <- converge(fit1323)
(500 - length(fit1323.con))/500
fit1323.con.ad <- admissible_check(fit1323.con)
(length(fit1323.con) - length(fit1323.con.ad))/length(fit1323.con)
gof1323 <- gof_check_bayes(fit1323.con.ad)
sum(gof1323[,1]<0.05)/length(fit1323.con.ad)
AR.p1323 <- AR.para(fit1323.con.ad, truepara)
save(fit1323, fit1323.con.ad, gof1323, AR.p1323, file = "f1323.RData")


## Case 8 fit1324 ####
## n=100, cate=7, asy, Bayes I
fit1324 <- list()
for (i in 1:length(dat132)) {
  dat <- dat132[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit1324[[i]] <- result
  print(paste("Rep:", i))
}
fit1324.con <- converge(fit1324)
(500 - length(fit1324.con))/500
fit1324.con.ad <- admissible_check(fit1324.con)
(length(fit1324.con) - length(fit1324.con.ad))/length(fit1324.con)
gof1324 <- gof_check_bayes(fit1324.con.ad)
sum(gof1324[,1]<0.05)/length(fit1324.con.ad)
AR.p1324 <- AR.para(fit1324.con.ad, truepara)
save(fit1324, fit1324.con.ad, gof1324, AR.p1324, file = "f1324.RData")
colMeans(gof1324)


# dat211 dat212 ####
## Case 1 fit2111 ####
## n=400, cate=3, sy, RML
fit2111 <- list()
for (i in 1:length(dat211)) {
  dat <- dat211[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit2111[[i]] <- result
  print(paste("Rep:", i))
}

fit2111.con <- converge(fit2111)
(500-length(fit2111.con))/500

fit2111.con.ad <- admissible_check(fit2111.con)
(length(fit2111.con)- length(fit2111.con.ad))/length(fit2111.con)
gof2111 <- gof_check(fit2111.con.ad)
sum(gof2111[,1]<0.05)/length(fit2111.con.ad)
AR.p2111 <- AR.para(fit2111.con.ad, truepara)
save(fit2111, fit2111.con.ad, gof2111, AR.p2111, file = "f2111.RData")
colMeans(gof2111)

## Case 2 fit2112 ####
## n=400, cate=3, sy, DWLS
fit2112 <- list()
for (i in 1:length(dat211)) {
  dat <- dat211[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit2112[[i]] <- result
  print(paste("Rep:", i))
}

fit2112.con <- converge(fit2112)
(500 - length(fit2112.con))/500
fit2112.con.ad <- admissible_check(fit2112.con)
(length(fit2112.con) - length(fit2112.con.ad))/length(fit2112.con)
gof2112 <- gof_check(fit2112.con.ad)
sum(gof2112[,1]<0.05)/length(fit2112.con.ad)
AR.p2112 <- AR.para(fit2112.con.ad, truepara)
save(fit2112, fit2112.con.ad, gof2112, AR.p2112, file = "f2112.RData")
colMeans(gof2112)

## Case 3 fit2113 ####
## n=400, cate=3, sy, Bayes NI
fit2113 <- list()
for (i in 1:length(dat211)) {
  dat <- dat211[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit2113[[i]] <- result
  print(paste("Rep:", i))
}
fit2113.con <- converge(fit2113)
(500 - length(fit2113.con))/500
fit2113.con.ad <- admissible_check(fit2113.con)
(length(fit2113.con) - length(fit2113.con.ad))/length(fit2113.con)
gof2113 <- gof_check_bayes(fit2113.con.ad)
sum(gof2113[,1]<0.05)/length(fit2113.con.ad)
AR.p2113 <- AR.para(fit2113.con.ad, truepara)
save(fit2113, fit2113.con.ad, gof2113, AR.p2113, file = "f2113.RData")

## Case 4 fit2114 ####
## n=400, cate=3, sy, Bayes I
fit2114 <- list()
for (i in 1:length(dat211)) {
  dat <- dat211[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit2114[[i]] <- result
  print(paste("Rep:", i))
}
fit2114.con <- converge(fit2114)
(500 - length(fit2114.con))/500
fit2114.con.ad <- admissible_check(fit2114.con)
(length(fit2114.con) - length(fit2114.con.ad))/length(fit2114.con)
gof2114 <- gof_check_bayes(fit2114.con.ad)
sum(gof2114[,1]<0.05)/length(fit2114.con.ad)
AR.p2114 <- AR.para(fit2114.con.ad, truepara)
save(fit2114, fit2114.con.ad, gof2114, AR.p2114, file = "f2114.RData")
colMeans(gof2114)


## Case 5 fit2121 ####
## n=400, cate=3, asy, RML
fit2121 <- list()
for (i in 1:length(dat212)) {
  dat <- dat212[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit2121[[i]] <- result
  print(paste("Rep:", i))
}

fit2121.con <- converge(fit2121)
(500-length(fit2121.con))/500
fit2121.con.ad <- admissible_check(fit2121.con)
(length(fit2121.con)- length(fit2121.con.ad))/length(fit2121.con)
gof2121 <- gof_check(fit2121.con.ad)
sum(gof2121[,1]<0.05)/length(fit2121.con.ad)
AR.p2121 <- AR.para(fit2121.con.ad, truepara)
save(fit2121, fit2121.con.ad, gof2121, AR.p2121, file = "f2121.RData")
colMeans(gof2121)

## Case 6 fit2122 ####
## n=400, cate=3, asy, DWLS
fit2122 <- list()
for (i in 1:length(dat212)) {
  dat <- dat212[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit2122[[i]] <- result
  print(paste("Rep:", i))
}

fit2122.con <- converge(fit2122)
(500 - length(fit2122.con))/500
fit2122.con.ad <- admissible_check(fit2122.con)
(length(fit2122.con) - length(fit2122.con.ad))/length(fit2122.con)
gof2122 <- gof_check(fit2122.con.ad)
sum(gof2122[,1]<0.05)/length(fit2122.con.ad)
AR.p2122 <- AR.para(fit2122.con.ad, truepara)
save(fit2122, fit2122.con.ad, gof2122, AR.p2122, file = "f2122.RData")
colMeans(gof2122)


## Case 7 fit2123 ####
## n=400, cate=3, asy, Bayes NI
fit2123 <- list()
for (i in 1:length(dat212)) {
  dat <- dat212[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit2123[[i]] <- result
  print(paste("Rep:", i))
}
fit2123.con <- converge(fit2123)
(500 - length(fit2123.con))/500
fit2123.con.ad <- admissible_check(fit2123.con)
(length(fit2123.con) - length(fit2123.con.ad))/length(fit2123.con)
gof2123 <- gof_check_bayes(fit2123.con.ad)
sum(gof2123[,1]<0.05)/length(fit2123.con.ad)
AR.p2123 <- AR.para(fit2123.con.ad, truepara)
save(fit2123, fit2123.con.ad, gof2123, AR.p2123, file = "f2123.RData")

## Case 8 fit2124 ####
## n=400, cate=3, asy, Bayes I
fit2124 <- list()
for (i in 1:length(dat212)) {
  dat <- dat212[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit2124[[i]] <- result
  print(paste("Rep:", i))
}
fit2124.con <- converge(fit2124)
(500 - length(fit2124.con))/500
fit2124.con.ad <- admissible_check(fit2124.con)
(length(fit2124.con) - length(fit2124.con.ad))/length(fit2124.con)
gof2124 <- gof_check_bayes(fit2124.con.ad)
sum(gof2124[,1]<0.05)/length(fit2124.con.ad)
AR.p2124 <- AR.para(fit2124.con.ad, truepara)
save(fit2124, fit2124.con.ad, gof2124, AR.p2124, file = "f2124.RData")
colMeans(gof2124)




# dat221 dat222 ####
## Case 1 fit2211 ####
## n=400, cate=5, sy, RML
fit2211 <- list()
for (i in 1:length(dat221)) {
  dat <- dat221[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit2211[[i]] <- result
  print(paste("Rep:", i))
}

fit2211.con <- converge(fit2211)
(500-length(fit2211.con))/500

fit2211.con.ad <- admissible_check(fit2211.con)
(length(fit2211.con)- length(fit2211.con.ad))/length(fit2211.con)
gof2211 <- gof_check(fit2211.con.ad)
sum(gof2211[,1]<0.05)/length(fit2211.con.ad)
AR.p2211 <- AR.para(fit2211.con.ad, truepara)
save(fit2211, fit2211.con.ad, gof2211, AR.p2211, file = "f2211.RData")
colMeans(gof2211)


## Case 2 fit2212 ####
## n=400, cate=5, sy, DWLS
fit2212 <- list()
for (i in 1:length(dat221)) {
  dat <- dat221[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit2212[[i]] <- result
  print(paste("Rep:", i))
}

fit2212.con <- converge(fit2212)
(500 - length(fit2212.con))/500
fit2212.con.ad <- admissible_check(fit2212.con)
(length(fit2212.con) - length(fit2212.con.ad))/length(fit2212.con)
gof2212 <- gof_check(fit2212.con.ad)
sum(gof2212[,1]<0.05)/length(fit2212.con.ad)
AR.p2212 <- AR.para(fit2212.con.ad, truepara)
save(fit2212, fit2212.con.ad, gof2212, AR.p2212, file = "f2212.RData")
colMeans(gof2212)


## Case 3 fit2213 ####
## n=400, cate=5, sy, Bayes NI
fit2213 <- list()
for (i in 1:length(dat221)) {
  dat <- dat221[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit2213[[i]] <- result
  print(paste("Rep:", i))
}
fit2213.con <- converge(fit2213)
(500 - length(fit2213.con))/500
fit2213.con.ad <- admissible_check(fit2213.con)
(length(fit2213.con) - length(fit2213.con.ad))/length(fit2213.con)
gof2213 <- gof_check_bayes(fit2213.con.ad)
sum(gof2213[,1]<0.05)/length(fit2213.con.ad)
AR.p2213 <- AR.para(fit2213.con.ad, truepara)
save(fit2213, fit2213.con.ad, AR.p2213, gof2213, file = "f2213.RData")

## Case 4 fit2214 ####
## n=400, cate=5, sy, Bayes I
fit2214 <- list()
for (i in 1:length(dat221)) {
  dat <- dat221[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit2214[[i]] <- result
  print(paste("Rep:", i))
}
fit2214.con <- converge(fit2214)
(500 - length(fit2214.con))/500
fit2214.con.ad <- admissible_check(fit2214.con)
(length(fit2214.con) - length(fit2214.con.ad))/length(fit2214.con)
gof2214 <- gof_check_bayes(fit2214.con.ad)
sum(gof2214[,1]<0.05)/length(fit2214.con.ad)
AR.p2214 <- AR.para(fit2214.con.ad, truepara)
save(fit2214, fit2214.con.ad, gof2214, AR.p2214, file = "f2214.RData")
colMeans(gof2214)


## Case 5 fit2221 ####
## n=400, cate=5, asy, RML
fit2221 <- list()
for (i in 1:length(dat222)) {
  dat <- dat222[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit2221[[i]] <- result
  print(paste("Rep:", i))
}

fit2221.con <- converge(fit2221)
(500-length(fit2221.con))/500

fit2221.con.ad <- admissible_check(fit2221.con)
(length(fit2221.con)- length(fit2221.con.ad))/length(fit2221.con)
gof2221 <- gof_check(fit2221.con.ad)
sum(gof2221[,1]<0.05)/length(fit2221.con.ad)
AR.p2221 <- AR.para(fit2221.con.ad, truepara)
save(fit2221, fit2221.con.ad, gof2221, AR.p2221, file = "f2221.RData")
colMeans(gof2221)


## Case 6 fit2222 ####
## n=400, cate=5, asy, DWLS
fit2222 <- list()
for (i in 1:length(dat222)) {
  dat <- dat222[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit2222[[i]] <- result
  print(paste("Rep:", i))
}

fit2222.con <- converge(fit2222)
(500 - length(fit2222.con))/500
fit2222.con.ad <- admissible_check(fit2222.con)
(length(fit2222.con) - length(fit2222.con.ad))/length(fit2222.con)
gof2222 <- gof_check(fit2222.con.ad)
sum(gof2222[,1]<0.05)/length(fit2222.con.ad)
AR.p2222 <- AR.para(fit2222.con.ad, truepara)
save(fit2222, fit2222.con.ad, gof2222, AR.p2222, file = "f2222.RData")
colMeans(gof2222)


## Case 7 fit2223 ####
## n=400, cate=5, asy, Bayes NI
fit2223 <- list()
for (i in 1:length(dat222)) {
  dat <- dat222[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit2223[[i]] <- result
  print(paste("Rep:", i))
}
fit2223.con <- converge(fit2223)
(500 - length(fit2223.con))/500
fit2223.con.ad <- admissible_check(fit2223.con)
(length(fit2223.con) - length(fit2223.con.ad))/length(fit2223.con)
gof2223 <- gof_check_bayes(fit2223.con.ad)
sum(gof2223[,1]<0.05)/length(fit2223.con.ad)
AR.p2223 <- AR.para(fit2223.con.ad, truepara)
save(fit2223, fit2223.con.ad, gof2223, AR.p2223, file = "f2223.RData")

## Case 8 fit2224 ####
## n=400, cate=5, asy, Bayes I
fit2224 <- list()
for (i in 1:length(dat222)) {
  dat <- dat222[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit2224[[i]] <- result
  print(paste("Rep:", i))
}
fit2224.con <- converge(fit2224)
(500 - length(fit2224.con))/500
fit2224.con.ad <- admissible_check(fit2224.con)
(length(fit2224.con) - length(fit2224.con.ad))/length(fit2224.con)
gof2224 <- gof_check_bayes(fit2224.con.ad)
sum(gof2224[,1]<0.05)/length(fit2224.con.ad)
AR.p2224 <- AR.para(fit2224.con.ad, truepara)
save(fit2224, fit2224.con.ad, gof2224, AR.p2224, file = "f2224.RData")
colMeans(gof2224)



# dat231 dat232 ####
## Case 1 fit2311 ####
## n=400, cate=7, sy, RML
fit2311 <- list()
for (i in 1:length(dat231)) {
  dat <- dat231[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit2311[[i]] <- result
  print(paste("Rep:", i))
}

fit2311.con <- converge(fit2311)
(500-length(fit2311.con))/500

fit2311.con.ad <- admissible_check(fit2311.con)
(length(fit2311.con)- length(fit2311.con.ad))/length(fit2311.con)
gof2311 <- gof_check(fit2311.con.ad)
sum(gof2311[,1]<0.05)/length(fit2311.con.ad)
AR.p2311 <- AR.para(fit2311.con.ad, truepara)
save(fit2311, fit2311.con.ad, gof2311, AR.p2311, file = "f2311.RData")
colMeans(gof2311)

## Case 2 fit2312 ####
## n=400, cate=7, sy, DWLS
fit2312 <- list()
for (i in 1:length(dat231)) {
  dat <- dat231[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit2312[[i]] <- result
  print(paste("Rep:", i))
}

fit2312.con <- converge(fit2312)
(500 - length(fit2312.con))/500
fit2312.con.ad <- admissible_check(fit2312.con)
(length(fit2312.con) - length(fit2312.con.ad))/length(fit2312.con)
gof2312 <- gof_check(fit2312.con.ad)
sum(gof2312[,1]<0.05)/length(fit2312.con.ad)
AR.p2312 <- AR.para(fit2312.con.ad, truepara)
save(fit2312, fit2312.con.ad, gof2312, AR.p2312, file = "f2312.RData")
colMeans(gof2312)

## Case 3 fit2313 ####
## n=400, cate=7, sy, Bayes NI
fit2313 <- list()
for (i in 1:length(dat231)) {
  dat <- dat231[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit2313[[i]] <- result
  print(paste("Rep:", i))
}
fit2313.con <- converge(fit2313)
(500 - length(fit2313.con))/500
fit2313.con.ad <- admissible_check(fit2313.con)
(length(fit2313.con) - length(fit2313.con.ad))/length(fit2313.con)
gof2313 <- gof_check_bayes(fit2313.con.ad)
sum(gof2313[,1]<0.05)/length(fit2313.con.ad)
AR.p2313 <- AR.para(fit2313.con.ad, truepara)
save(fit2313, fit2313.con.ad, gof2313, AR.p2313, file = "f2313.RData")

## Case 4 fit2314 ####
## n=400, cate=7, sy, Bayes I
fit2314 <- list()
for (i in 1:length(dat231)) {
  dat <- dat231[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit2314[[i]] <- result
  print(paste("Rep:", i))
}
fit2314.con <- converge(fit2314)
(500 - length(fit2314.con))/500
fit2314.con.ad <- admissible_check(fit2314.con)
(length(fit2314.con) - length(fit2314.con.ad))/length(fit2314.con)
gof2314 <- gof_check_bayes(fit2314.con.ad)
sum(gof2314[,1]<0.05)/length(fit2314.con.ad)
AR.p2314 <- AR.para(fit2314.con.ad, truepara)
save(fit2314, fit2314.con.ad, gof2314, AR.p2314, file = "f2314.RData")
colMeans(gof2314)



## Case 5 fit2321 ####
## n=400, cate=7, asy, RML
fit2321 <- list()
for (i in 1:length(dat232)) {
  dat <- dat232[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit2321[[i]] <- result
  print(paste("Rep:", i))
}

fit2321.con <- converge(fit2321)
(500-length(fit2321.con))/500
fit2321.con.ad <- admissible_check(fit2321.con)
(length(fit2321.con)- length(fit2321.con.ad))/length(fit2321.con)
gof2321 <- gof_check(fit2321.con.ad)
sum(gof2321[,1]<0.05)/length(fit2321.con.ad)
AR.p2321 <- AR.para(fit2321.con.ad, truepara)
save(fit2321, fit2321.con.ad, gof2321, AR.p2321, file = "f2321.RData")
colMeans(gof2321)


## Case 6 fit2322 ####
## n=400, cate=7, asy, DWLS
fit2322 <- list()
for (i in 1:length(dat232)) {
  dat <- dat232[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit2322[[i]] <- result
  print(paste("Rep:", i))
}

fit2322.con <- converge(fit2322)
(500 - length(fit2322.con))/500
fit2322.con.ad <- admissible_check(fit2322.con)
(length(fit2322.con) - length(fit2322.con.ad))/length(fit2322.con)
gof2322 <- gof_check(fit2322.con.ad)
sum(gof2322[,1]<0.05)/length(fit2322.con.ad)
AR.p2322 <- AR.para(fit2322.con.ad, truepara)
save(fit2322, fit2322.con.ad, gof2322, AR.p2322, file = "f2322.RData")
colMeans(gof2322)


## Case 7 fit2323 ####
## n=400, cate=7, asy, Bayes NI
fit2323 <- list()
for (i in 1:length(dat232)) {
  dat <- dat232[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit2323[[i]] <- result
  print(paste("Rep:", i))
}
fit2323.con <- converge(fit2323)
(500 - length(fit2323.con))/500
fit2323.con.ad <- admissible_check(fit2323.con)
(length(fit2323.con) - length(fit2323.con.ad))/length(fit2323.con)
gof2323 <- gof_check_bayes(fit2323.con.ad)
sum(gof2323[,1]<0.05)/length(fit2323.con.ad)
AR.p2323 <- AR.para(fit2323.con.ad, truepara)
save(fit2323, fit2323.con.ad, gof2323, AR.p2323, file = "f2323.RData")

## Case 8 fit2324 ####
## n=400, cate=7, asy, Bayes I
fit2324 <- list()
for (i in 1:length(dat232)) {
  dat <- dat232[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit2324[[i]] <- result
  print(paste("Rep:", i))
}
fit2324.con <- converge(fit2324)
(500 - length(fit2324.con))/500
fit2324.con.ad <- admissible_check(fit2324.con)
(length(fit2324.con) - length(fit2324.con.ad))/length(fit2324.con)
gof2324 <- gof_check_bayes(fit2324.con.ad)
sum(gof2324[,1]<0.05)/length(fit2324.con.ad)
AR.p2324 <- AR.para(fit2324.con.ad, truepara)
save(fit2324, fit2324.con.ad, gof2324, AR.p2324, file = "f2324.RData")
colMeans(gof2324)


# dat311 dat312 ####
## Case 1 fit3111 ####
## n=1000, cate=3, sy, RML
fit3111 <- list()
for (i in 1:length(dat311)) {
  dat <- dat311[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit3111[[i]] <- result
  print(paste("Rep:", i))
}

fit3111.con <- converge(fit3111)
(500-length(fit3111.con))/500
fit3111.con.ad <- admissible_check(fit3111.con)
(length(fit3111.con)- length(fit3111.con.ad))/length(fit3111.con)
gof3111 <- gof_check(fit3111.con.ad)
sum(gof3111[,1]<0.05)/length(fit3111.con.ad)
AR.p3111 <- AR.para(fit3111.con.ad, truepara)
save(fit3111, fit3111.con.ad, gof3111, AR.p3111, file = "f3111.RData")
colMeans(gof3111)



## Case 2 fit3112 ####
## n=1000, cate=3, sy, DWLS
fit3112 <- list()
for (i in 1:length(dat311)) {
  dat <- dat311[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit3112[[i]] <- result
  print(paste("Rep:", i))
}

fit3112.con <- converge(fit3112)
(500 - length(fit3112.con))/500
fit3112.con.ad <- admissible_check(fit3112.con)
(length(fit3112.con) - length(fit3112.con.ad))/length(fit3112.con)
gof3112 <- gof_check(fit3112.con.ad)
sum(gof3112[,1]<0.05)/length(fit3112.con.ad)
AR.p3112 <- AR.para(fit3112.con.ad, truepara)
save(fit3112, fit3112.con.ad, gof3112, AR.p3112, file = "f3112.RData")
colMeans(gof3112)



## Case 3 fit3113 ####
## n=1000, cate=3, sy, Bayes NI
fit3113 <- list()
for (i in 1:length(dat311)) {
  dat <- dat311[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit3113[[i]] <- result
  print(paste("Rep:", i))
}
fit3113.con <- converge(fit3113)
(500 - length(fit3113.con))/500
fit3113.con.ad <- admissible_check(fit3113.con)
(length(fit3113.con) - length(fit3113.con.ad))/length(fit3113.con)
gof3113 <- gof_check_bayes(fit3113.con.ad)
sum(gof3113[,1]<0.05)/length(fit3113.con.ad)
AR.p3113 <- AR.para(fit3113.con.ad, truepara)
save(fit3113, fit3113.con.ad, AR.p3113, gof3113, file = "f3113.RData")

## Case 4 fit3114 ####
## n=1000, cate=3, sy, Bayes I
fit3114 <- list()
for (i in 1:length(dat311)) {
  dat <- dat311[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit3114[[i]] <- result
  print(paste("Rep:", i))
}
fit3114.con <- converge(fit3114)
(500 - length(fit3114.con))/500
fit3114.con.ad <- admissible_check(fit3114.con)
(length(fit3114.con) - length(fit3114.con.ad))/length(fit3114.con)
gof3114 <- gof_check_bayes(fit3114.con.ad)
sum(gof3114[,1]<0.05)/length(fit3114.con.ad)
AR.p3114 <- AR.para(fit3114.con.ad, truepara)
save(fit3114, fit3114.con.ad, gof3114, AR.p3114, file = "f3114.RData")
colMeans(gof3114)


## Case 5 fit3121 ####
## n=1000, cate=3, asy, RML
fit3121 <- list()
for (i in 1:length(dat312)) {
  dat <- dat312[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit3121[[i]] <- result
  print(paste("Rep:", i))
}

fit3121.con <- converge(fit3121)
(500-length(fit3121.con))/500

fit3121.con.ad <- admissible_check(fit3121.con)
(length(fit3121.con)- length(fit3121.con.ad))/length(fit3121.con)
gof3121 <- gof_check(fit3121.con.ad)
sum(gof3121[,1]<0.05)/length(fit3121.con.ad)
AR.p3121 <- AR.para(fit3121.con.ad, truepara)
save(fit3121, fit3121.con.ad, gof3121, AR.p3121, file = "f3121.RData")
colMeans(gof3121)



## Case 6 fit3122 ####
## n=1000, cate=3, asy, DWLS
fit3122 <- list()
for (i in 1:length(dat312)) {
  dat <- dat312[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit3122[[i]] <- result
  print(paste("Rep:", i))
}

fit3122.con <- converge(fit3122)
(500 - length(fit3122.con))/500
fit3122.con.ad <- admissible_check(fit3122.con)
(length(fit3122.con) - length(fit3122.con.ad))/length(fit3122.con)
gof3122 <- gof_check(fit3122.con.ad)
sum(gof3122[,1]<0.05)/length(fit3122.con.ad)
AR.p3122 <- AR.para(fit3122.con.ad, truepara)
save(fit3122, fit3122.con.ad, gof3122, AR.p3122, file = "f3122.RData")
colMeans(gof3122)


## Case 7 fit3123 ####
## n=1000, cate=3, asy, Bayes NI
fit3123 <- list()
for (i in 1:length(dat312)) {
  dat <- dat312[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit3123[[i]] <- result
  print(paste("Rep:", i))
}
fit3123.con <- converge(fit3123)
(500 - length(fit3123.con))/500
fit3123.con.ad <- admissible_check(fit3123.con)
(length(fit3123.con) - length(fit3123.con.ad))/length(fit3123.con)
gof3123 <- gof_check_bayes(fit3123.con.ad)
sum(gof3123[,1]<0.05)/length(fit3123.con.ad)
AR.p3123 <- AR.para(fit3123.con.ad, truepara)
save(fit3123, fit3123.con.ad, AR.p3123, gof3123, file = "f3123.RData")

## Case 8 fit3124 ####
## n=1000, cate=3, asy, Bayes I
fit3124 <- list()
for (i in 1:length(dat312)) {
  dat <- dat312[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit3124[[i]] <- result
  print(paste("Rep:", i))
}
fit3124.con <- converge(fit3124)
(500 - length(fit3124.con))/500
fit3124.con.ad <- admissible_check(fit3124.con)
(length(fit3124.con) - length(fit3124.con.ad))/length(fit3124.con)
gof3124 <- gof_check_bayes(fit3124.con.ad)
sum(gof3124[,1]<0.05)/length(fit3124.con.ad)
AR.p3124 <- AR.para(fit3124.con.ad, truepara)
save(fit3124, fit3124.con.ad, gof3124, AR.p3124, file = "f3124.RData")
colMeans(gof3124)


# dat321 dat322 ####
## Case 1 fit3211 ####
## n=1000, cate=5, sy, RML
fit3211 <- list()
for (i in 1:length(dat321)) {
  dat <- dat321[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit3211[[i]] <- result
  print(paste("Rep:", i))
}

fit3211.con <- converge(fit3211)
(500-length(fit3211.con))/500
fit3211.con.ad <- admissible_check(fit3211.con)
(length(fit3211.con)- length(fit3211.con.ad))/length(fit3211.con)
gof3211 <- gof_check(fit3211.con.ad)
sum(gof3211[,1]<0.05)/length(fit3211.con.ad)
AR.p3211 <- AR.para(fit3211.con.ad, truepara)
save(fit3211, fit3211.con.ad, gof3211, AR.p3211, file = "f3211.RData")
colMeans(gof3211)

## Case 2 fit3212 ####
## n=1000, cate=5, sy, DWLS
fit3212 <- list()
for (i in 1:length(dat321)) {
  dat <- dat321[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit3212[[i]] <- result
  print(paste("Rep:", i))
}

fit3212.con <- converge(fit3212)
(500 - length(fit3212.con))/500
fit3212.con.ad <- admissible_check(fit3212.con)
(length(fit3212.con) - length(fit3212.con.ad))/length(fit3212.con)
gof3212 <- gof_check(fit3212.con.ad)
sum(gof3212[,1]<0.05)/length(fit3212.con.ad)
AR.p3212 <- AR.para(fit3212.con.ad, truepara)
save(fit3212, fit3212.con.ad, gof3212, AR.p3212, file = "f3212.RData")
colMeans(gof3212)


## Case 3 fit3213 ####
## n=1000, cate=5, sy, Bayes NI
fit3213 <- list()
for (i in 1:length(dat321)) {
  dat <- dat321[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit3213[[i]] <- result
  print(paste("Rep:", i))
}
fit3213.con <- converge(fit3213)
(500 - length(fit3213.con))/500
fit3213.con.ad <- admissible_check(fit3213.con)
(length(fit3213.con) - length(fit3213.con.ad))/length(fit3213.con)
gof3213 <- gof_check_bayes(fit3213.con.ad)
sum(gof3213[,1]<0.05)/length(fit3213.con.ad)
AR.p3213 <- AR.para(fit3213.con.ad, truepara)
save(fit3213, fit3213.con.ad, AR.p3213, gof3213, file = "f3213.RData")

## Case 4 fit3214 ####
## n=1000, cate=5, sy, Bayes I
fit3214 <- list()
for (i in 1:length(dat321)) {
  dat <- dat321[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit3214[[i]] <- result
  print(paste("Rep:", i))
}
fit3214.con <- converge(fit3214)
(500 - length(fit3214.con))/500
fit3214.con.ad <- admissible_check(fit3214.con)
(length(fit3214.con) - length(fit3214.con.ad))/length(fit3214.con)
gof3214 <- gof_check_bayes(fit3214.con.ad)
sum(gof3214[,1]<0.05)/length(fit3214.con.ad)
AR.p3214 <- AR.para(fit3214.con.ad, truepara)
save(fit3214, fit3214.con.ad, gof3214, AR.p3214, file = "f3214.RData")
colMeans(gof3214)

## Case 5 fit3221 ####
## n=1000, cate=5, asy, RML
fit3221 <- list()
for (i in 1:length(dat322)) {
  dat <- dat322[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit3221[[i]] <- result
  print(paste("Rep:", i))
}

fit3221.con <- converge(fit3221)
(500-length(fit3221.con))/500
fit3221.con.ad <- admissible_check(fit3221.con)
(length(fit3221.con)- length(fit3221.con.ad))/length(fit3221.con)
gof3221 <- gof_check(fit3221.con.ad)
sum(gof3221[,1]<0.05)/length(fit3221.con.ad)
AR.p3221 <- AR.para(fit3221.con.ad, truepara)
save(fit3221, fit3221.con.ad, gof3221, AR.p3221, file = "f3221.RData")
colMeans(gof3221)


## Case 6 fit3222 ####
## n=1000, cate=5, asy, DWLS
fit3222 <- list()
for (i in 1:length(dat322)) {
  dat <- dat322[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit3222[[i]] <- result
  print(paste("Rep:", i))
}

fit3222.con <- converge(fit3222)
(500 - length(fit3222.con))/500
fit3222.con.ad <- admissible_check(fit3222.con)
(length(fit3222.con) - length(fit3222.con.ad))/length(fit3222.con)
gof3222 <- gof_check(fit3222.con.ad)
sum(gof3222[,1]<0.05)/length(fit3222.con.ad)
AR.p3222 <- AR.para(fit3222.con.ad, truepara)
save(fit3222, fit3222.con.ad, gof3222, AR.p3222, file = "f3222.RData")
colMeans(gof3222)

## Case 7 fit3223 ####
## n=1000, cate=5, asy, Bayes NI
fit3223 <- list()
for (i in 1:length(dat322)) {
  dat <- dat322[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit3223[[i]] <- result
  print(paste("Rep:", i))
}
fit3223.con <- converge(fit3223)
(500 - length(fit3223.con))/500
fit3223.con.ad <- admissible_check(fit3223.con)
(length(fit3223.con) - length(fit3223.con.ad))/length(fit3223.con)
gof3223 <- gof_check_bayes(fit3223.con.ad)
sum(gof3223[,1]<0.05)/length(fit3223.con.ad)
AR.p3223 <- AR.para(fit3223.con.ad, truepara)
save(fit3223, fit3223.con.ad, AR.p3223, gof3223, file = "f3223.RData")

## Case 8 fit3224 ####
## n=1000, cate=5, asy, Bayes I
fit3224 <- list()
for (i in 1:length(dat322)) {
  dat <- dat322[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit3224[[i]] <- result
  print(paste("Rep:", i))
}
fit3224.con <- converge(fit3224)
(500 - length(fit3224.con))/500
fit3224.con.ad <- admissible_check(fit3224.con)
(length(fit3224.con) - length(fit3224.con.ad))/length(fit3224.con)
gof3224 <- gof_check_bayes(fit3224.con.ad)
sum(gof3224[,1]<0.05)/length(fit3224.con.ad)
AR.p3224 <- AR.para(fit3224.con.ad, truepara)
save(fit3224, fit3224.con.ad, gof3224, AR.p3224, file = "f3224.RData")
colMeans(gof3224)


# dat331 dat332 ####
## Case 1 fit3311 ####
## n=1000, cate=7, sy, RML
fit3311 <- list()
for (i in 1:length(dat331)) {
  dat <- dat331[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit3311[[i]] <- result
  print(paste("Rep:", i))
}

fit3311.con <- converge(fit3311)
(500-length(fit3311.con))/500

fit3311.con.ad <- admissible_check(fit3311.con)
(length(fit3311.con)- length(fit3311.con.ad))/length(fit3311.con)
gof3311 <- gof_check(fit3311.con.ad)
sum(gof3311[,1]<0.05)/length(fit3311.con.ad)
AR.p3311 <- AR.para(fit3311.con.ad, truepara)
save(fit3311, fit3311.con.ad, gof3311, AR.p3311, file = "f3311.RData")
colMeans(gof3311)

## Case 2 fit3312 ####
## n=1000, cate=7, sy, DWLS
fit3312 <- list()
for (i in 1:length(dat331)) {
  dat <- dat331[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit3312[[i]] <- result
  print(paste("Rep:", i))
}

fit3312.con <- converge(fit3312)
(500 - length(fit3312.con))/500
fit3312.con.ad <- admissible_check(fit3312.con)
(length(fit3312.con) - length(fit3312.con.ad))/length(fit3312.con)
gof3312 <- gof_check(fit3312.con.ad)
sum(gof3312[,1]<0.05)/length(fit3312.con.ad)
AR.p3312 <- AR.para(fit3312.con.ad, truepara)
save(fit3312, fit3312.con.ad, gof3312, AR.p3312, file = "f3312.RData")
colMeans(gof3312)


## Case 3 fit3313 ####
## n=1000, cate=7, sy, Bayes NI
fit3313 <- list()
for (i in 1:length(dat331)) {
  dat <- dat331[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit3313[[i]] <- result
  print(paste("Rep:", i))
}
fit3313.con <- converge(fit3313)
(500 - length(fit3313.con))/500
fit3313.con.ad <- admissible_check(fit3313.con)
(length(fit3313.con) - length(fit3313.con.ad))/length(fit3313.con)
gof3313 <- gof_check_bayes(fit3313.con.ad)
sum(gof3313[,1]<0.05)/length(fit3313.con.ad)
AR.p3313 <- AR.para(fit3313.con.ad, truepara)
save(fit3313, fit3313.con.ad, gof3313, AR.p3313, file = "f3313.RData")

## Case 4 fit3314 ####
## n=1000, cate=7, sy, Bayes I
fit3314 <- list()
for (i in 1:length(dat331)) {
  dat <- dat331[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit3314[[i]] <- result
  print(paste("Rep:", i))
}
fit3314.con <- converge(fit3314)
(500 - length(fit3314.con))/500
fit3314.con.ad <- admissible_check(fit3314.con)
(length(fit3314.con) - length(fit3314.con.ad))/length(fit3314.con)
gof3314 <- gof_check_bayes(fit3314.con.ad)
sum(gof3314[,1]<0.05)/length(fit3314.con.ad)
AR.p3314 <- AR.para(fit3314.con.ad, truepara)
save(fit3314, fit3314.con.ad, gof3314, AR.p3314, file = "f3314.RData")
colMeans(gof3314)

## Case 5 fit3321 ####
## n=1000, cate=7, asy, RML
fit3321 <- list()
for (i in 1:length(dat332)) {
  dat <- dat332[[i]]
  result <- cfa(
    analysisModel,
    data = dat,
    mimic = 'Mplus',
    estimator = 'MLMV',
    se = "robust",
    std.lv = TRUE
  )
  fit3321[[i]] <- result
  print(paste("Rep:", i))
}

fit3321.con <- converge(fit3321)
(500-length(fit3321.con))/500
fit3321.con.ad <- admissible_check(fit3321.con)
(length(fit3321.con)- length(fit3321.con.ad))/length(fit3321.con)
gof3321 <- gof_check(fit3321.con.ad)
sum(gof3321[,1]<0.05)/length(fit3321.con.ad)
AR.p3321 <- AR.para(fit3321.con.ad, truepara)
save(fit3321, fit3321.con.ad, gof3321, AR.p3321, file = "f3321.RData")
colMeans(gof3321)



## Case 6 fit3322 ####
## n=1000, cate=7, asy, DWLS
fit3322 <- list()
for (i in 1:length(dat332)) {
  dat <- dat332[[i]]
  result <- cfa(
    analysisModel,
    data  = dat,
    mimic = 'Mplus',
    test = "satorra.bentler",
    se = "robust",
    std.lv = TRUE,
    ordered = TRUE 
  )
  fit3322[[i]] <- result
  print(paste("Rep:", i))
}

fit3322.con <- converge(fit3322)
(500 - length(fit3322.con))/500
fit3322.con.ad <- admissible_check(fit3322.con)
(length(fit3322.con) - length(fit3322.con.ad))/length(fit3322.con)
gof3322 <- gof_check(fit3322.con.ad)
sum(gof3322[,1]<0.05)/length(fit3322.con.ad)
AR.p3322 <- AR.para(fit3322.con.ad, truepara)
save(fit3322, fit3322.con.ad, gof3322, AR.p3322, file = "f3322.RData")
colMeans(gof3322)



## Case 7 fit3323 ####
## n=1000, cate=7, asy, Bayes NI
fit3323 <- list()
for (i in 1:length(dat332)) {
  dat <- dat332[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(theta = "gamma(3, 0.1)[sd]",
                              ibpsi = "wishart(5, iden)"),
                 n.chains = 2)
  fit3323[[i]] <- result
  print(paste("Rep:", i))
}
fit3323.con <- converge(fit3323)
(500 - length(fit3323.con))/500
fit3323.con.ad <- admissible_check(fit3323.con)
(length(fit3323.con) - length(fit3323.con.ad))/length(fit3323.con)
gof3323 <- gof_check_bayes(fit3323.con.ad)
sum(gof3323[,1]<0.05)/length(fit3323.con.ad)
AR.p3323 <- AR.para(fit3323.con.ad, truepara)
save(fit3323, fit3323.con.ad, AR.p3323, gof3323, file = "f3323.RData")

## Case 8 fit3324 ####
## n=1000, cate=7, asy, Bayes I
fit3324 <- list()
for (i in 1:length(dat332)) {
  dat <- dat332[[i]]
  result <- bcfa(analysisModel,
                 data  = dat,
                 mcmcfile = TRUE,
                 std.lv = TRUE,
                 cp = "srs",
                 dp = dpriors(lambda = "normal(0.6, 0.2)", # 0.65+/-1.96*sqrt(0.02)
                              theta = "gamma(1, 2)[prec]",
                              rho = "beta(3, 8)"),
                 n.chains = 2)
  fit3324[[i]] <- result
  print(paste("Rep:", i))
}
fit3324.con <- converge(fit3324)
(500 - length(fit3324.con))/500
fit3324.con.ad <- admissible_check(fit3324.con)
(length(fit3324.con) - length(fit3324.con.ad))/length(fit3324.con)
gof3324 <- gof_check_bayes(fit3324.con.ad)
sum(gof3324[,1]<0.05)/length(fit3324.con.ad)
AR.p3324 <- AR.para(fit3324.con.ad, truepara)
save(fit3324, fit3324.con.ad, gof3324, AR.p3324, file = "f3324.RData")
colMeans(gof3324)


