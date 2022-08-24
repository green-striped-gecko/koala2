
library(tidyverse)
library(jagsUI)

unlogit <- function(x) exp(x) / (1 + exp(x))

# Import data
dat <- read.csv("loctrees_sub_trees_08072021.csv")
glimpse(dat)

# using SITE_VISIT_ID
sub.dat <- dat %>%
  dplyr::rename(site = SITE_VISIT_ID,
                kfp = KFP,
                species = SPECIES,
                tree_num = TREE_NUM) %>%
  filter(!is.na(period) & !is.na(kfp)) %>%
  select(site, year, subarea, tree_num, kfp, species, period) %>%
  mutate(plot = paste(site, year))

glimpse(sub.dat)

# number of occupied sites per period
no <- sub.dat %>%
  group_by(plot, period) %>%
  summarise(kfp = sum(kfp)) %>%
  mutate(kfp = ifelse(kfp > 0, 1, 0)) %>%
  group_by(period) %>%
  summarise(occ = sum(kfp),
            n = n()) %>%
  mutate(prop.occ = occ/n)
no

# number of trees with kfps detected by plot
kfp.dat <- sub.dat %>%
  mutate(kfp.pres = ifelse(kfp > 0, 1, 0)) %>%
  group_by(plot, subarea, period) %>%
  summarise(n.search = n(),
            n.det = sum(kfp.pres))

glimpse(kfp.dat)

# number of trees searched per plot
table(kfp.dat$n.search)

# barplot of number of detections at each plot where kfps were detected
barplot(table(factor(kfp.dat$n.det, levels = 1:15)))

# data for JAGS model
y <- kfp.dat$n.det
period <- kfp.dat$period
subarea <- kfp.dat$subarea
N <- length(y)

########################################################################################
# Royle-Nichols model with negative binomial for period only

mod <- "model {
  for(i in 1:N) {
    y[i] ~ dbin(p[i], 30)
    p[i] <- 1 - (1 - r) ^ K[i]
    K[i] ~ dnegbin(prob_nb[i], theta)
    prob_nb[i] <- theta / (theta + mu[i])
    log(mu[i]) <- b[subarea[i], period[i]]
  }

  for(i in 1:5) {
    for(j in 1:3) {
      b[i, j] ~ dnorm(0, 0.1)
    }
  }
  logit(r) <- lr
  lr ~ dnorm(0, 0.1)
  theta ~ dunif(0, 100)

}"  

write(mod, "model2.txt")

# Sample model in JAGS with three chains
mod.jags <- jags(model = "model2.txt",
                 data = list(y = y, N = N, period = period, subarea = subarea),
                 inits = function() list(K = y), 
                 param = c("b", "lr", "theta", "K"),
                 n.chains = 3,
                 n.iter =20000,
                 n.burnin = 5000,
                 parallel = T)  

out <- mod.jags$summary
out[1:20, c(1, 2, 3, 7, 8, 9)]


##################################################################################  
# by period
# proportion occupied (if K = 0 then there were no koalas)
# extract all the K values from the MCMC runs
tout <- mod.jags$sims.list$K
dim(tout)

# convert to presence/absence of koala
tout <- ifelse(tout > 0, 1, 0)

# number of plots per period 
nsites <- table(period)
nsites

occ <- matrix(nrow = nrow(tout), ncol = 3)
prop.occ <- matrix(nrow = nrow(tout), ncol = 3)

# calculate number of occupied plots for each period  
for(i in 1:nrow(tout)) {
  occ[i, ] <- tapply(tout[i, ], period, sum)
  prop.occ[i, ] <- occ[i, ] / nsites
}

# convert to a useable form
po <- data.frame(prop.occ = apply(prop.occ, 2, mean),
                 lcl = apply(prop.occ, 2, function(x) quantile(x, 0.025)),
                 ucl = apply(prop.occ, 2, function(x) quantile(x, 0.975)))

po$period <- 1:3

# plot with raw occupancy in red  
ggplot(po, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.4) +
  theme_classic()  

##################################################################################  
# by subarea and period

# number of plots per period per subarea  
nsites <- table(subarea, period)
nsites

occ <- array(dim = c(5, 3, nrow(tout)))
prop.occ <- array(dim = c(5, 3, nrow(tout)))

# calculate number of occupied plots for each subarea and period  
for(i in 1:nrow(tout)) {
  occ[, , i] <- tapply(tout[i, ], list(subarea, period), sum)
  prop.occ[, , i] <- occ[, , i] / nsites
}

# convert to a useable form
po <- data.frame(occ = apply(prop.occ, c(1, 2), mean),
                 lcl = apply(prop.occ, c(1, 2), function(x) quantile(x, 0.025)),
                 ucl = apply(prop.occ, c(1, 2), function(x) quantile(x, 0.975)))

glimpse(po)
po


prop.occ <- c(po$occ.1, po$occ.2, po$occ.3)
lcl <- c(po$lcl.1, po$lcl.2, po$lcl.3)
ucl <- c(po$ucl.1, po$ucl.2, po$ucl.3)

eo <- data.frame(prop.occ = prop.occ,
                 lcl = lcl, 
                 ucl = ucl,
                 subarea = rep(1:5, 3),
                 period = rep(1:3, each = 5))
eo


# raw occupancy  
ro <- sub.dat %>%
  group_by(plot, subarea, period) %>%
  summarise(kfp = sum(kfp)) %>%
  mutate(kfp = ifelse(kfp > 0, 1, 0)) %>%
  group_by(subarea, period) %>%
  summarise(occ = sum(kfp),
            n = n()) %>%
  mutate(prop.occ = occ/n)
ro


ggplot(eo, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  geom_point(data = ro, col = "red", size = 3) +
  facet_wrap(~subarea) +
  theme_bw()

