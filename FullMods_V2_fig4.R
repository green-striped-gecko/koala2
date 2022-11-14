
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
n.tree <- kfp.dat$n.search
N <- length(y)
N

########################################################################################
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
                  n.iter = 10000,
                  n.burnin = 2000,
                  parallel = T)  

out1 <- mod1.jags$summary
out1[1:20, c(1, 2, 3, 7, 8, 9)]

################################################################################
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

# plot with raw occupancy in red  
ggplot(po1, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.4) +
  theme_classic()  

##############################################################################################
##############################################################################################
# same model but searching only the top 5 trees

a <- table(sub.dat$species, sub.dat$kfp)
prop <- round(a[, 2] / rowSums(a), 4)
a <- cbind(a, prop)
a <- a[order(-a[, 2]), ]
a

# top 5 trees
spp <- rownames(a)[1:5]
spp

# observed detections

# number of trees searched and number with kfps detected by plot
kfp.dat <- sub.dat %>%
  mutate(pref = ifelse(species %in% spp, 1, 0),
         det.pref = ifelse(species %in% spp &  kfp == 1, 1, 0)) %>%
  group_by(plot, period) %>%
  summarise(n.search = sum(pref),
            n.det = sum(det.pref),
            total.det = sum(kfp))

glimpse(kfp.dat)

# what was missed?
table(kfp.dat$n.det, kfp.dat$total.det)

# barplot of number of detections at each plot where kfps were detected
barplot(table(factor(kfp.dat$n.det, levels = 1:15)))

# number of detections as a function of number of trees searched
ggplot(kfp.dat, aes(y = n.det, x = n.search)) +
  geom_count() +
  theme_bw()

# distribution of number of trees searched
barplot(table(factor(kfp.dat$n.search, levels = 0:30)))

# proportion of plots detected as a function of number searched
pd <- kfp.dat %>%
  mutate(det = ifelse(n.det > 0, 1, 0)) %>%
  group_by(n.search) %>%
  summarise(prop.det = sum(det) / n(),
            n = n())

ggplot(pd, aes(y = prop.det, x = n.search)) +
  geom_point(aes(size = n)) +
  theme_bw()


# data for JAGS model
y <- kfp.dat$n.det
period <- kfp.dat$period
n.tree <- kfp.dat$n.search
N <- length(y)
N

# Sample model in JAGS with three chains
mod2.jags <- jags(model = "model1.txt",
                  data = list(y = y, N = N, period = period, n.tree = n.tree),
                  inits = function() list(K = y), 
                  param = c("b", "lr", "theta", "K"),
                  n.chains = 3,
                  n.iter = 10000,
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
ggplot(po2, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.4) +
  theme_classic()  

# compare models
ggplot(po1, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  
  geom_point(data = po2, aes(x = period + 0.1), size = 6, colour = "blue") +
  geom_line(data = po2, aes(x = period + 0.1), colour = "blue") +
  geom_errorbar(data = po2, aes(ymin = lcl, ymax = ucl, x = period + 0.1), width = 0, lwd = 1.1, colour = "blue") +
  
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.4) +
  theme_classic()  

################################################################################
#################################################################################
# Model with stopping all data

# list of plots
plot_list <- sub.dat %>%
  select(site, year, subarea, period, plot) %>%
  unique()

# number of trees searched
nt <- sub.dat %>%
  group_by(plot) %>%
  summarise(n.search = n())

# identify the first kfp occurrence tree and number of trees kfps were found under
frst <- sub.dat %>%
  group_by(plot) %>%
  mutate(ind = kfp == 1 & !duplicated(kfp == 1)) %>%
  summarise(first_tree = tree_num[ind],
            n_tree = sum(kfp))

# join to plot list
ft <- full_join(plot_list, frst) %>%
  full_join(nt)

glimpse(ft)

# Data for JAGS model
Tmax <- ft$n.search
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
                            "z"),
                  n.chains = 3,
                  n.iter = 20000,
                  n.burnin = 5000,
                  parallel = T)  

out3 <- mod3.jags$summary
out3[1:20, c(1, 2, 3, 7, 8, 9)]

#########################################################################################
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
ggplot(po3, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.5) +
  theme_classic()  

##############################################################################################
# compare occupancy probabilities

ggplot(po1, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  
  geom_point(data = po2, aes(x = period + 0.1), size = 6, colour = "blue") +
  geom_line(data = po2, aes(x = period + 0.1), colour = "blue") +
  geom_errorbar(data = po2, aes(ymin = lcl, ymax = ucl, x = period + 0.1), width = 0, lwd = 1.1, colour = "blue") +
  
  geom_point(data = po3, aes(x = period + 0.2), size = 6, colour = "red") +
  geom_line(data = po3, aes(x = period + 0.2), colour = "red") +
  geom_errorbar(data = po3, aes(ymin = lcl, ymax = ucl, x = period + 0.2), width = 0, lwd = 1.1, colour = "red") +
  
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.5) +
  theme_classic()  

################################################################################
#################################################################################
# Model with stopping preferred trees only

# list of plots
plot_list <- sub.dat %>%
  select(site, year, subarea, period, plot) %>%
  unique()

# number of preferred trees searched
nt <- sub.dat %>%
  mutate(pref = ifelse(species %in% spp, 1, 0)) %>%
  group_by(plot) %>%
  summarise(n.search = sum(pref))

# identify the first kfp occurrence tree and number of trees kfps were found under
frst <- sub.dat %>%
  filter(species %in% spp) %>%
  group_by(plot) %>%
  mutate(tree_num = row_number(),
         ind = kfp == 1 & !duplicated(kfp == 1)) %>%
  summarise(first_tree = tree_num[ind],
            n_tree = sum(kfp))

# join to plot list
ft <- full_join(plot_list, frst) %>%
  full_join(nt)

glimpse(ft)

# Data for JAGS model
Tmax <- ft$n.search
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

# Sample model in JAGS with three chains
mod4.jags <- jags(model = "model2.txt",
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
                            "z"),
                  n.chains = 3,
                  n.iter = 20000,
                  n.burnin = 5000,
                  parallel = T)  

out4 <- mod4.jags$summary
out4[1:20, c(1, 2, 3, 7, 8, 9)]

#########################################################################################
# proportion occupied 
tout <- mod4.jags$sims.list$z
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
po4 <- data.frame(prop.occ = apply(prop.occ, 2, mean),
                  lcl = apply(prop.occ, 2, function(x) quantile(x, 0.025)),
                  ucl = apply(prop.occ, 2, function(x) quantile(x, 0.975)))

po4$period <- 1:3


# plot with raw occupancy in red  
ggplot(po4, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.6) +
  theme_classic()  

##############################################################################################
# compare occupancy probabilities

ggplot(po1, aes(y = prop.occ, x = period)) +
  geom_point(size = 6) +
  geom_line() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, lwd = 1.1) +
  
  geom_point(data = po2, aes(x = period + 0.1), size = 6, colour = "blue") +
  geom_line(data = po2, aes(x = period + 0.1), colour = "blue") +
  geom_errorbar(data = po2, aes(ymin = lcl, ymax = ucl, x = period + 0.1), width = 0, lwd = 1.1, colour = "blue") +
  
  geom_point(data = po3, aes(x = period + 0.2), size = 6, colour = "red") +
  geom_line(data = po3, aes(x = period + 0.2), colour = "red") +
  geom_errorbar(data = po3, aes(ymin = lcl, ymax = ucl, x = period + 0.2), width = 0, lwd = 1.1, colour = "red") +
  
  geom_point(data = po4, aes(x = period + 0.3), size = 6, colour = "orange") +
  geom_line(data = po4, aes(x = period + 0.3), colour = "orange") +
  geom_errorbar(data = po4, aes(ymin = lcl, ymax = ucl, x = period + 0.3), width = 0, lwd = 1.1, colour = "orange") +
  
  geom_point(data = no, col = "red", size = 3) +
  ylim(0, 0.6) +
  theme_classic()  


#Code to generate Figure 4, Report Koala2
#From Ewen Lawler 

####### All four models occupancy 

# add method column to each df
po1$method <- "RGb-SAT"
po2$method <- "Top5-SAT"
po3$method <- "Rapid-SAT"
po4$method <- "Top5 & Rapid-SAT"

# create df for merge from raw data
no$method <- "Raw data"
no_formerge <- data.frame(prop.occ = no$prop.occ, lcl = 0, ucl = 0,
                          period = no$period, method = no$method)

# merge to one df (ordered to match with other plot)
all_mods <- rbind(po1, po3, po2, po4, no_formerge)
all_mods
# factor order
all_mods$method <- factor(all_mods$method,levels = c("RGb-SAT", 
                                                     "Rapid-SAT",
                                                     "Top5-SAT",
                                                     "Top5 & Rapid-SAT", 
                                                     "Raw data"))
# add a periods object for x axis labels
periods <- c("2009-2011", "2012-2015", "2016-2019")

occu.plot <- ggplot(data = all_mods, aes(y = prop.occ, 
                                         x = period, 
                                         colour = method,
                                         shape = method))+
  geom_point(position = position_dodge(width = 0.4), size = 5)+
  geom_line(position = position_dodge(width = 0.4), size = 1)+
  geom_errorbar(aes(ymin = lcl, ymax =ucl, width = 0), size = 1.5,
                position = position_dodge(width = 0.4))+
  labs(x = "Period", y = "Proportion occupied")+
  theme_classic()+
  scale_x_continuous(breaks=seq(1,3,1),labels=periods)+
  scale_colour_manual(name = "Survey method",
                      labels = c("RGb-SAT",
                                 "Rapid-SAT",
                                 "Top5-SAT",
                                 "Top5 & Rapid-SAT", 
                                 "Raw data"),
                      values = c("black", "red", "blue", "purple", "grey54")) +   
  scale_shape_manual(name = "Survey method",
                     labels = c("RGb-SAT",
                                "Rapid-SAT",
                                "Top5-SAT",
                                "Top5 & Rapid-SAT", 
                                "Raw data"),
                     values = c(16,16,16,16,18))+
  theme(legend.position = "bottom")
print(occu.plot)

###plot 5


