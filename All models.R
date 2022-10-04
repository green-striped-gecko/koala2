
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
  mutate(plot = paste(site, year)) %>%
  arrange(plot, tree_num)

glimpse(sub.dat)
head(sub.dat)

table(sub.dat$kfp)

# number of plots and trees searched per plot
table(table(sub.dat$plot))

#######################################################################################
# first model
# using all data

# number of trees with kfps detected by plot
kfp.dat <- sub.dat %>%
  mutate(kfp.pres = ifelse(kfp > 0, 1, 0)) %>%
  group_by(plot, subarea, period) %>%
  summarise(n.search = n(),
            n.det = sum(kfp.pres),
            n.tree = n())

glimpse(kfp.dat)

# number of trees searched per plot
table(kfp.dat$n.search)

# barplot of number of detections at each plot where kfps were detected
barplot(table(factor(kfp.dat$n.det, levels = 1:15)))

# data for JAGS model
y <- kfp.dat$n.det
period <- kfp.dat$period
n.tree <- kfp.dat$n.tree
N <- length(y)
N

# no of detections
table(ifelse(y > 0, 1, 0))

#-----------------------------------------------------------------------------------------
# Royle-Nichols model with negative binomial 

mod1 <- "model {
  for(i in 1:N) {
    y[i] ~ dbin(p[i], n.tree[i])
    p[i] <- 1 - (1 - r) ^ K[i]
    K[i] ~ dnegbin(prob_nb[i], theta)
    prob_nb[i] <- theta / (theta + mu[i])
    log(mu[i]) <- b[period[i]]
  }

  for(i in 1:3) {
    b[i] ~ dnorm(0, 0.1)
  }
  logit(r) <- lr
  lr ~ dnorm(0, 0.1)
  theta ~ dunif(0, 100)

}"  

write(mod1, "model1.txt")

# Sample model in JAGS with three chains
mod1.jags <- jags(model = "model1.txt",
                 data = list(y = y, N = N, period = period, n.tree = n.tree),
                 inits = function() list(K = y), 
                 param = c("b", "lr", "theta", "K"),
                 n.chains = 3,
                 n.iter =10000,
                 n.burnin = 2000,
                 parallel = T)  

out1 <- mod1.jags$summary
out1[1:20, c(1, 2, 3, 7, 8, 9)]

#-----------------------------------------------------------------------------------------
# proportion occupied (if K = 0 then there were no koalas)
# extract all the K values from the MCMC runs
tout <- mod1.jags$sims.list$K
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
po1 <- data.frame(prop.occ = apply(prop.occ, 2, mean),
                 lcl = apply(prop.occ, 2, function(x) quantile(x, 0.025)),
                 ucl = apply(prop.occ, 2, function(x) quantile(x, 0.975)))

po1$period <- 1:3

# observed number of occupied sites per period
no <- sub.dat %>%
  group_by(plot, period) %>%
  summarise(kfp = sum(kfp)) %>%
  mutate(kfp = ifelse(kfp > 0, 1, 0)) %>%
  group_by(period) %>%
  summarise(occ = sum(kfp),
            n = n()) %>%
  mutate(prop.occ = occ/n)
no

# plot with raw occupancy in red  
ggplot(po1, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.5) +
  theme_classic()  

#################################################################################
#################################################################################
# fit the same model as above but using only preferred trees

a <- table(sub.dat$species, sub.dat$kfp)
prop <- round(a[, 2] / rowSums(a), 4)
a <- cbind(a, prop)
a <- a[order(-a[, 2]), ]
a

# top 5 trees
spp <- rownames(a)[1:5]
spp

# number of kfp found under preferred and less-preferred trees
kfp.dat <- sub.dat %>%
  mutate(pref = ifelse(species %in% spp, 1, 0)) %>%
  group_by(plot, subarea, period) %>%
  summarise(n.pref = sum(pref),
            n.det = sum(kfp),
            n.det.pref = sum(kfp[pref == 1])) %>%
  mutate(n.det.nopref = n.det - n.det.pref)

glimpse(kfp.dat)

# frequency distribution of preferred trees per plot
barplot(table(factor(kfp.dat$n.pref, levels = 0:30)))

# detections under preferred and less-preferred trees
table(kfp.dat$n.det.pref, kfp.dat$n.det.nopref)

# probability of finding kfp as a function of number of preferred trees
pd <- kfp.dat %>%
  group_by(n.pref) %>%
  summarise(prob.det = sum(ifelse(n.det.pref > 0, 1, 0)/n()))

ggplot(pd, aes(y = prob.det, x = n.pref)) +
  geom_point(size = 3) +
  theme_bw()

# barplot of number of detections at each plot where kfps were detected
barplot(table(factor(kfp.dat$n.det.pref, levels = 1:15)))

#-----------------------------------------------------------------------------------------
# data for JAGS model
y <- kfp.dat$n.det.pref
period <- kfp.dat$period
N <- length(y)
N

# number of preferred trees searched per plot
n.tree <- kfp.dat$n.pref

# Sample model in JAGS with three chains
mod2.jags <- jags(model = "model1.txt",
                  data = list(y = y, N = N, period = period, n.tree = n.tree),
                  inits = function() list(K = y), 
                  param = c("b", "lr", "theta", "K"),
                  n.chains = 3,
                  n.iter =10000,
                  n.burnin = 2000,
                  parallel = T)  

out2 <- mod2.jags$summary
out2[1:20, c(1, 2, 3, 7, 8, 9)]

################################################################################
# proportion occupied (if K = 0 then there were no koalas)
# extract all the K values from the MCMC runs
tout <- mod2.jags$sims.list$K
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
po2 <- data.frame(prop.occ = apply(prop.occ, 2, mean),
                  lcl = apply(prop.occ, 2, function(x) quantile(x, 0.025)),
                  ucl = apply(prop.occ, 2, function(x) quantile(x, 0.975)))

po2$period <- 1:3

# plot with raw occupancy in red  
ggplot(po1, aes(y = prop.occ, x = period)) +
  geom_point(data = po1, aes(x = period), size = 6) +
  geom_line(data = po1, aes(x = period)) +
  geom_errorbar(data = po1, aes(x = period, ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  
  geom_point(data = po2, aes(x = period + 0.1), size = 6, colour = "blue") +
  geom_line(data = po2, aes(x = period + 0.1), colour = "blue") +
  geom_errorbar(data = po2, aes(x = period + 0.1, ymin = lcl, ymax = ucl), width = 0, lwd = 1.1, colour = "blue") +
  
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.5) +
  theme_classic()  

#################################################################################
#################################################################################
# Models with stopping at first detection
# first using all data

# list of plots
plot_list <- sub.dat %>%
  select(site, year, subarea, period, plot) %>%
  unique()

# identify the first kfp occurrence tree and number of trees kfps were found under
frst <- sub.dat %>%
  group_by(plot) %>%
  mutate(ind = kfp == 1 & !duplicated(kfp == 1)) %>%
  summarise(first_tree = tree_num[ind],
            n_tree = sum(kfp))

# number of trees searched
nt <- sub.dat %>%
  group_by(plot) %>%
  summarise(n_search = n())

# join to plot list
ft <- full_join(plot_list, frst) %>%
  full_join(nt)
glimpse(ft)
head(ft)

# number of plots with observed detection
table(is.na(ft$first_tree))

# number of trees versus first tree
mft <- ft %>%
  group_by(first_tree) %>%
  summarise(mnt = mean(n_tree),
            n = n()) %>%
  drop_na()

mft

par(mfrow = c(2, 2))
plot(ft$n_tree ~ ft$first_tree, pch = 19, col = "grey")
  points(mft$mnt ~ mft$first_tree, pch = 19, cex = 1.2)
  
hist(ft$first_tree, breaks = 0:30, xlab = "first_tree")


#-----------------------------------------------------------------------------------------
# Data for JAGS model
Tmax <- ft$n_search
ttd <- ft$first_tree

# Non-detections considered censored
isCensored <- ifelse(is.na(ttd) == T, 1, 0)

# Occupany only known if detection occured
z <- ifelse(is.na(ttd) == F, 1, NA)

# N
N <- length(ttd)

# TTD of non-detections initialised greater than search maximum
ttdst <- Tmax + 1
ttdst[isCensored == 0] <- NA

period <- ft$period

# Weibull occupancy-detection model with censoring ------------------------

mod2 <- "model {
  for(i in 1:N) {
    # Occupancy as Bernoulli variable
    z[i] ~ dbern(psi[period[i]])

    # Detection times follow Weibull distribution
    ttd[i] ~ dweib(shape, lambda)
 
    # Model for censoring 
    isCensored[i] ~ dbern(theta[i])
    theta[i] <- z[i] * step(ttd[i] - Tmax[i]) + (1 - z[i])

  }

  # Uninformative priors
  for(i in 1:3) {
    psi[i] ~ dunif(0, 1)
  }
  lambda ~ dgamma(0.01, 0.01)
  shape ~ dgamma(0.01, 0.01)

  # Estimated detection probability within search period (if 30 trees searched)
  D.hat <- 1 - exp(-lambda * 30 ^ shape)

  # Estimated number of occupancies
  N_occ <- sum(z[])

}"

write(mod2, "model2.txt")

# Sample model in JAGS with three chains
mod3.jags <- jags(model = "model2.txt",
                 data = list(ttd = ttd, 
                             isCensored = isCensored, 
                             Tmax = Tmax, 
                             N = N, 
                             z = z,
                             period = period),
                 inits = function() list(ttd = ttdst),
                 param = c("lambda", 
                           "shape", 
                           "psi", 
                           "N_occ", 
                           "D.hat", 
                           "z"),
                 n.chains = 3,
                 n.iter = 20000,
                 n.burnin = 5000,
                 parallel = T)  

out3 <- mod3.jags$summary
out3[1:20, c(1, 2, 3, 7, 8, 9)]

#-----------------------------------------------------------------------------------------
# proportion occupied 
tout <- mod3.jags$sims.list$z
dim(tout)

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
po3 <- data.frame(prop.occ = apply(prop.occ, 2, mean),
                 lcl = apply(prop.occ, 2, function(x) quantile(x, 0.025)),
                 ucl = apply(prop.occ, 2, function(x) quantile(x, 0.975)))

po3$period <- 1:3

# plot with raw occupancy in red  
ggplot(po1, aes(y = prop.occ, x = period)) +
  geom_point(data = po1, aes(x = period), size = 6) +
  geom_line(data = po1, aes(x = period)) +
  geom_errorbar(data = po1, aes(x = period, ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  
  geom_point(data = po3, aes(x = period + 0.1), size = 6, colour = "blue") +
  geom_line(data = po3, aes(x = period + 0.1), colour = "blue") +
  geom_errorbar(data = po3, aes(x = period + 0.1, ymin = lcl, ymax = ucl), width = 0, lwd = 1.1, colour = "blue") +
  
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.5) +
  theme_classic()  

#################################################################################
#################################################################################
# Models with stopping at first detection
# only preferred trees

# list of plots
plot_list <- sub.dat %>%
  filter(species )
  select(site, year, subarea, period, plot) %>%
  unique()

# identify the first kfp occurrence tree and number of trees kfps were found under
frst <- sub.dat %>%
  group_by(plot) %>%
  mutate(ind = kfp == 1 & !duplicated(kfp == 1)) %>%
  summarise(first_tree = tree_num[ind],
            n_tree = sum(kfp))

# number of trees searched
nt <- sub.dat %>%
  group_by(plot) %>%
  summarise(n_search = n())

# join to plot list
ft <- full_join(plot_list, frst) %>%
  full_join(nt)
glimpse(ft)
head(ft)

# number of plots with observed detection
table(is.na(ft$first_tree))

# number of trees versus first tree
mft <- ft %>%
  group_by(first_tree) %>%
  summarise(mnt = mean(n_tree),
            n = n()) %>%
  drop_na()

mft

par(mfrow = c(2, 2))
plot(ft$n_tree ~ ft$first_tree, pch = 19, col = "grey")
points(mft$mnt ~ mft$first_tree, pch = 19, cex = 1.2)

hist(ft$first_tree, breaks = 0:30, xlab = "first_tree")


#-----------------------------------------------------------------------------------------
# Data for JAGS model
Tmax <- ft$n_search
ttd <- ft$first_tree

# Non-detections considered censored
isCensored <- ifelse(is.na(ttd) == T, 1, 0)

# Occupany only known if detection occured
z <- ifelse(is.na(ttd) == F, 1, NA)

# N
N <- length(ttd)

# TTD of non-detections initialised greater than search maximum
ttdst <- Tmax + 1
ttdst[isCensored == 0] <- NA

period <- ft$period

