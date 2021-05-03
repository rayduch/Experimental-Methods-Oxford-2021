############# Week 1 ##############
## Programmed by: Mats Ahrenshop ##
###### Date: 29 April 2020 ########

# Verification of Young (2019) #

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
library(tidyverse)
library(ggpubr)

set.seed(103648)

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


## experiment data

dat <- read.csv("young19.csv")

## create subset of actual wristband days

datw <- subset(dat, dat$wristband_real == 1)



#----- Difference-in-means estimator -----

dim_hg <- mean(dat$prob_act_st[dat$treat_assign == "TG"],
               na.rm = T) - mean(dat$prob_act_st[dat$treat_assign == "C"],
                                 na.rm = T)

dim_hp <- mean(dat$prob_act_st[dat$treat_assign == "TP"],
               na.rm = T) - mean(dat$prob_act_st[dat$treat_assign == "C"],
                                 na.rm = T)

dim_bg <- mean(dat$wristband[dat$treat_assign == "TG"],
               na.rm = T) - mean(dat$wristband[dat$treat_assign == "C"],
                                 na.rm = T)

dim_bp <- mean(dat$wristband[dat$treat_assign == "TP"],
               na.rm = T) - mean(dat$wristband[dat$treat_assign == "C"],
                                 na.rm = T)

options(scipen = 999, digits = 6)

#---- OLS ----
m_hg <- lm(prob_act_st ~ treat_assign, data = dat[dat$treat_assign != "TP", ])
summary(m_hg) # exact same

m_hp <- lm(prob_act_st ~ treat_assign, data = dat[dat$treat_assign != "TG", ])
summary(m_hp) # exact same

m_bg <- lm(wristband ~ treat_assign, data = dat[dat$treat_assign != "TP", ])
summary(m_bg) # exact same

m_bp <- lm(wristband ~ treat_assign, data = dat[dat$treat_assign != "TG", ])
summary(m_bp) # exact same


#---- Reproduction of Table 3 ----

## which includes weights as 1/p and heteroskedastic robust SEs (treatment groups are of unequal size)

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
## but SEs are exactly reproduced; we'll do that next week

## let's collect all estimates from the different estimators and compare

ates <- tibble(
  outcome = rep(c("hypothetical", "behavioral"), each = 8),
  treatment = rep(rep(c("general", "political"), each = 4), times = 2),
  approach = rep(c("diff-in-means", "ols", "ols_weight", "table3"), times = 4),
  ate = c(round(dim_hg, 3), round(m_hg$coefficients[2], 3), round(mod1g$coefficients[2], 3), -0.0545,
          round(dim_hp, 3), round(m_hp$coefficients[2], 3), round(mod1p$coefficients[2], 3), -0.773,
          round(dim_bg, 3), round(m_bg$coefficients[2], 3), round(mod2g$coefficients[2], 3), -0.104,
          round(dim_bp, 3), round(m_bp$coefficients[2], 3), round(mod2p$coefficients[2], 3), -0.189),
  se = c("NA", round(coef(summary(m_hg))[2, 2], 3), round(mod1g$se[2], 3), 0.077,
         "NA", round(coef(summary(m_hp))[2, 2], 3), round(mod1p$se[2], 3), 0.080,
         "NA", round(coef(summary(m_bg))[2, 2], 3), round(mod2g$se[2], 3), 0.050,
         "NA", round(coef(summary(m_bp))[2, 2], 3), round(mod2p$se[2], 3), 0.053)
)


#---- Visualization of diff in means ----

# jitter plot with mean and CI's
p1 <- ggplot(data = dat,
       mapping = aes(x = treat_assign,
                     y = prob_act_st)) +
  geom_jitter(shape = 1, alpha = .45, height = 0, width = .25) +
  geom_point(stat = "summary", fun.y = "mean") +
  stat_summary(fun.data = mean_cl_normal, fun.args = list(conf.int = 0.95),
               geom = "errorbar", width = .1, size = .8) +
  xlab("") +
  theme_minimal()
p1

# CDFs with group means
con <- mean(dat$prob_act_st[dat$treat_assign == "C"],
            na.rm = T)
tg <- mean(dat$prob_act_st[dat$treat_assign == "TG"],
           na.rm = T)
tp <- mean(dat$prob_act_st[dat$treat_assign == "TP"],
           na.rm = T)

p2 <- ggplot(data = dat) +
  stat_density(mapping = aes(x = prob_act_st, group = treat_assign, color = treat_assign,
                             linetype = treat_assign), geom = "line", position = "identity",
               size = 1) +
  geom_vline(xintercept = con, color = "red", linetype = 3, size = .8) +
  geom_vline(xintercept = tg, color = "green", linetype = 3, size = .8) +
  geom_vline(xintercept = tp, color = "blue", linetype = 3, size = .8) +
  theme_minimal()
p2

# boxplots per group
p3 <- ggplot(data = dat) +
  geom_boxplot(mapping = aes(y = prob_act_st, x = treat_assign)) +
  xlab("") +
  ylab("") +
  theme_minimal()
p3

# arrange with panels
ggarrange(p2,
          ggarrange(p1, p3, ncol = 2),
          nrow = 2)


## FIGURE 2: Proportion of respondents very likely or sure to dissent during an election period

acts <- c('prob_shirt_elec_bin', 'prob_joke_elec_bin', 'prob_meeting_elec_bin',
          'prob_vh_elec_bin', 'prob_pungwe_elec_bin', 'prob_testify_elec_bin')

p <- melt(dat[,c('X', 'surveyor_id', 'treat_assign', acts)],
          id=c('X', 'surveyor_id', 'treat_assign'))

p$variable <- gsub("prob_", "", p$variable)
p$variable <- gsub("_elec_bin", "", p$variable)
p$variable <- Recode(p$variable, "'joke'='Joke about President';
                     'meeting'='Attend Opposition Meeting';'pungwe'='Refuse ZANU Meeting';
                     'shirt'='Wear Opposition Shirt';'testify'='Testify in Trial';
                     'vh'='Reveal Opposition to Leader'",
                     levels = c('Wear Opposition Shirt', 'Joke about President', 
                                'Attend Opposition Meeting', 'Refuse ZANU Meeting',
                                'Reveal Opposition to Leader', 'Testify in Trial'))

p <- subset(p, is.na(p$variable) == F)

ggplot(p, aes(treat_assign, value)) +
  geom_bar(stat = "summary", fun.y = "mean", fill = rep(c('grey', 'lightblue', 'darkred'), 6)) +
  stat_summary(fun.data = mean_cl_normal, fun.args = list(conf.int = 0.95), geom = "errorbar", 
               position = position_dodge(width = 0.90), width = 0.2) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion") + xlab("Treatment") +
  facet_wrap(~variable, ncol = 3)


