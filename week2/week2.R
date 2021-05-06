############# Week 2 ##############
## Programmed by: Mats Ahrenshop ##
####### Date: 04 May 2020 #########

# RANDOMIZATION INFERENCE #

library(tidyverse)
library(data.table)
library(gtools)
library(lubridate)
library(AER)
library(xtable)
library(pBrackets)
library(Hmisc)
library(car)
library(ri)
library(psych)
library(ggpubr)

set.seed(103648)

#---- PRELIMINARIES ----

## function to remove NAs

replace_NA <- function(data, varnames, NAvals) {
  numvars <- length(varnames)
  numNA <- length(NAvals)
  data <- as.data.frame(data)
  for (i in 1:numvars){
    varname <- varnames[i]
    for (j in 1:numNA) {
      data[,varname] <- ifelse(data[,varname]==NAvals[j], NA, data[,varname])
    }
  }
  return(data)
}

st_out <- function(x, data = dat, treat = 'treat_all') {
  (x - mean(x[data[,treat]==0],na.rm=T))/sd(x[data[,treat]==0],na.rm=T)
}

robust <- function(model) {
  X <- model.matrix(model)
  n <- dim(X)[1]
  k <- dim(X)[2]
  u <- matrix(resid(model))
  meat1 <- t(X) %*% diag(diag(crossprod(t(u)))) %*% X
  dfc <- n/(n-k) 
  se <- sqrt(dfc*diag(solve(crossprod(X)) %*% meat1 %*% solve(crossprod(X))))
  return(se)
}

#---- RI WITH TOY DATA ----

## let's start with a quick example with simulated data

Y <- sample(x = 1:10, size = 18, replace = TRUE)
Z <- sample(x = 0:1, size = 18, replace = TRUE)
block <- c(rep(1, 4), rep(2, 6), rep(3, 8))

perms <- genperms(Z, blockvar = block) # enumerate all possible permutations given assignment vector (permutation matrix)
probs <- genprobexact(Z, blockvar = block) # compute probability of treatment per block
ate <- estate(Y, Z, prob = probs) # estimate observed ATE

Ys <- genouts(Y, Z, ate = 0) # generate full schedule of potential outcomes under exact H0
distout <- gendist(Ys, perms, prob = probs) # generate randomization distribution
p <- mean(abs(distout) >= abs(ate))
dispdist(distout, ate)

## other simulation alternative:

N <- 50
m <- 25

d <- ifelse(1:N %in% sample(1:N, m), 1, 0)
Y0 <- runif(N, 0, 1)
Y1 <- Y0 + rnorm(N, 2, 2)
Y <- Y1*d + Y0*(1-d)

cbind(Y0, Y1, d, Y)	# look at your data

## Conduct analysis of actual experiment
## Estimate the ATE

# nonparametric
mean(Y[d == 1]) - mean(Y[d == 0])

# or fitting data to ols
lm(Y ~ d)

Z <- d
probs <- genprobexact(Z)
ate <- estate(Y, Z, prob = probs)
  
perms <- genperms(Z, maxiter = 10000)

# Create potential outcomes UNDER THE SHARP NULL OF NO EFFECT FOR ANY UNIT
Ys <- genouts(Y, Z, ate = 0)
  
# Generate the sampling distribution based on schedule of potential outcome
# implied by the sharp null hypothesis
distout <- gendist(Ys, perms, prob = probs)
  
sum(distout >= ate)                 # one-tailed comparison used to calculate p-value (greater than)
sum(abs(distout) >= abs(ate))       # two-tailed comparison used to calculate p-value
  
dispdist(distout, ate)               # display p-values, 95% confidence interval, standard error under the null, and graph the sampling distribution under the null

# estimation of confidence intervals assuming ATE=estimated ATE

Ys <- genouts(Y, Z, ate = ate)            # create potential outcomes UNDER THE ASSUMPTION THAT ATE=ESTIMATED ATE
distout <- gendist(Ys, perms, prob = probs)  # generate the sampling distribution  based on the schedule of potential outcomes implied by the null hypothesis
dispdist(distout, ate)               # display p-values, 95% confidence interval, standard error under the null, and graph the sampling distribution under the null



#---- YOUNG 2019 ----

## experiment data

dat <- read.csv("young19.csv")

## create subset of actual wristband days

datw <- subset(dat, dat$wristband_real == 1)

## models include weights as 1/p and heteroskedastic robust SEs (treatment groups are of unequal size)

## ATEs and robust SEs

mod1g <- lm(prob_act_st ~ treat_assign, data = dat[dat$treat_assign != "TP", ], 
            weights = dat$TG_inv[dat$treat_assign != "TP"])
mod1g$se <- robust(mod1g)

mod1p <- lm(prob_act_st ~ treat_assign, data = dat[dat$treat_assign != "TG", ], 
            weights = dat$TP_inv[dat$treat_assign != "TG"])
mod1p$se <- robust(mod1p)

mod2g <- lm(wristband ~ treat_assign, data = datw[datw$treat_assign != "TP", ], 
            weights = datw$TG_inv[datw$treat_assign != "TP"])
mod2g$se <- robust(mod2g)

mod2p <- lm(wristband ~ treat_assign, data = datw[datw$treat_assign != "TG", ], 
            weights = datw$TP_inv[datw$treat_assign != "TG"])
mod2p$se <- robust(mod2p)


## ATE estimates differ slightly from the ones reported in the paper; in paper estimated as part of RI routine
## but SEs are exactly reproduced

## now let's quickly revisit results in Young 2019 in light of randomization inference

## hypotheical -- general fear

t <- subset(dat, complete.cases(dat[, 'prob_act_st']) & dat$treat_assign %in% c('C', 'TG'))
mod1g$n <- dim(t)[1] # extract nrow
Z <- t$treat_all
Y <- t$prob_act_st
block <- t$block
probs <- genprobexact(Z = Z, blockvar = block) # generate probabilities of treatment by block
mod1g$ate <- estate(Y = Y, Z = Z, prob = probs) # calculate ate
perms <- genperms(Z, maxiter = 100000) # enumerate all possible ways of random assignment (permutations)
Ys <- genouts(Y = Y, Z = Z, ate = 0) # generate schedule of potential outcomes under exact H0 (ate = 0)
distout <- gendist(Ys, perms) # generate distribution of ATE's (100,000)
mod1g$p <- mean(abs(distout) >= abs(mod1g$ate)) # p-value 2-sided
dispdist(distout, mod1g$ate)

## behavioural -- general fear

t <- subset(datw, complete.cases(datw[, 'wristband']) & datw$treat_assign %in% c('C', 'TG'))
mod2g$n <- dim(t)[1]
Z <- t$treat_all
Y <- t$wristband
block <- t$block
probs <- genprobexact(Z, blockvar = block)
mod2g$ate <- estate(Y, Z, prob = probs)
perms <- genperms(Z, maxiter = 100000)
Ys <- genouts(Y, Z, ate = 0)
distout <- gendist(Ys, perms)
mod2g$p <- mean(abs(distout) >= abs(mod2g$ate))


#---- COOPERMAN 2017 ----

rm(list = ls(all = TRUE))

library(plm) # Created using version 1.6-5
library(dplyr) # Created using version 0.5.0
library(tidyr) # Created using version 0.6.1
library(lme4) # Created using version 1.1-12
library(lfe) # Created using version 2.5-1998
library(reshape2) # Created using version 1.4.2
library(Matrix) # Created using version 1.2-7.1
library(Formula) # Created using version 1.2-1


# The RI procedure generates 1000 permutations of the treatment assignment vector (ie. 1000 potential weather assignments ) for each specification and then calculates the ATE.

# Load Data --------
load("cooperman17.Rdata") 

# Subset data to drop 7 counties in Gomez et al. (2007) dataset that do not have out-of-sample rainfall data
data2 <- dataset[!is.na(dataset$Rainfall12), ]
data2 <- data2 %>% arrange(FIPS.County, Year)

nsims = 10 # Use nsims = 10 to run RI efficiently on 10 permutations of treatment assignment.
# nsims = 1000 # Paper uses n of 1000. To replicate full results, uncomment beginning of line. 

# Generate permutations of the treatment assignment vector.
# They assume all counties are independent.

assign_independently_inch <- function(df){
  df_sim <-
    df %>%
    group_by(FIPS.County) %>%
    mutate(Z_sim = sample(Rainfall12, n(), replace = T)) %>%
    arrange(FIPS.County)
  return(df_sim$Z_sim)
}

# 2. Output needed for Table 3 and Figure 4 -----
##################################################################

# Conduct RI procedure for Rainfall (inches)
# COUNTY ---------
# 1. Create potential weather assignments using out-of-sample rainfall data, assume counties independent
data2 <- data2 %>% arrange(FIPS.County, Year)
set.seed(1234)

# zindep is a matrix in which the columns correspond to different permutations of the treatment assignment vector generated by the function assign_independently_inch()
zindep <- matrix(NA, nrow = nrow(data2), ncol = nsims)
colnames(zindep) <- paste0("V", 1:nsims)
for (i in 1:nsims) {
  zindep[,i] <- assign_independently_inch(data2)
  print(i)
}

data.county <- cbind(data2, zindep)

# 2. Calculate estimates for ATE under different permutations of the treatment assignment vector
#load("Data_county.Rdata")
zsims.cols <- colnames(data.county)[grep("V", colnames(data.county), fixed = TRUE)]

# This loop estimates the random effects model from Gomez et al. (2007) utilizing each permutation of the treatment assignment vector. 
# It returns the estimated ATE under each permutation.
ate.sims.indep <- rep(NA, nsims)
for (i in 1:nsims) {
  model <- lmer(Turnout ~ get(zsims.cols[i]) + Snow  + ZPcntHSGrad + AdjIncome + PcntBlack + FarmsPerCap + Closing + Motor + Property + Literacy + PollTax  + GubElection + SenElection + Turnout.Lag + (1|FIPS.County) + as.factor(Year), data=data.county)
  ate.sims.indep[i] <- fixef(model)["get(zsims.cols[i])"]
  print(i)
}

p <- mean(abs(ate.sims.indep) >= abs(-1.052)) # result Table 3 row 1



