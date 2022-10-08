
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
# top 5 trees
a <- table(sub.dat$species, sub.dat$kfp)
prop <- round(a[, 2] / rowSums(a), 4)
a <- cbind(a, prop)
a <- a[order(-a[, 2]), ]
a
spp <- rownames(a)[1:5]
spp

# number of trees with kfps detected by plot
kfp.dat <- sub.dat %>%
  mutate(pref = ifelse(species %in% spp, 1, 0)) %>%
  group_by(plot, subarea, period) %>%
  summarise(n.search = n(),
            n.det = sum(kfp),
            n.pref = sum(pref), 
            n.pref.det = sum(kfp[pref == 1]))

glimpse(kfp.dat)

# number of kfp detected under preferred trees versus number detected under all trees
plot(kfp.dat$n.pref.det ~ kfp.dat$n.det)
abline(0, 1)

# total number of trees searched per plot
table(kfp.dat$n.search)

# number of preferred trees searched per plot
table(kfp.dat$n.pref)

# proportion of trees searched if only preferred trees are searcherd
sum(kfp.dat$n.pref) / sum(kfp.dat$n.search)

# barplot of number of detections at each plot 
barplot(table(factor(kfp.dat$n.det, levels = 0:15)))

# data for JAGS model
y <- kfp.dat$n.det
period <- kfp.dat$period
n.tree <- kfp.dat$n.search
n.pref <- kfp.dat$n.pref
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
                  param = c("b", "r", "theta", "K"),
                  n.chains = 3,
                  n.iter = 10000,
                  n.burnin = 2000,
                  parallel = T)  

out1 <- mod1.jags$summary
out1[1:20, c(1, 2, 3, 7, 8, 9)]

#----------------------------------------------------------------------------------------------------
# get parameters from model

# detection probability for a tree
r <- out1[4, 1]
# parameters of the K negative binomial distribution
size <- out1[5, 1]
mu <- exp(mean(out1[1:3, 1]))

# a plot is unoccupied when K = 0
# for a negative binomial distribution, the probability of a zero value, and hence
# the probability of non-occupancy is:
# (size / (size + mu))^size
# occupancy probability is one minus this:
# occ <- 1 - (size / (size + mu))^size

# Hence, for specified occupancy probability and size: mu = (size / (1-occ)^(1/size)) - size

#--------------------------------------------------------------------------------------
# function to generate observed data for a plot with n.tree searched given
# specified occupancy probability, r and size

gen.plot <- function(n.tree, occ, r, size) {
  # mu value given specified occupancy
  mu <- (size / (1-occ)^(1/size)) - size
  
  # generate K value for the plot
  K <- rnbinom(n = 1, size = size, mu = mu)
  
  # probability of detection
  p <- 1 - (1-r)^K
  
  # observed outcome for the n.tree searched
  y <- rbinom(n = n.tree, size = 1, prob = p)
  return(y)
}

#--------------------------------------------------------------------------------------

# generate data for N plots
N <- 500

# occ = 0.3
gen.out <- replicate(N, gen.plot(n.tree = 30, occ = 0.3, r = r, size = size))
dim(gen.out)

# number of kfps per plot
kfp <- colSums(gen.out)
table(kfp)

# occ = 0.1
gen.out <- replicate(N, gen.plot(n.tree = 30, occ = 0.1, r = r, size = size))
kfp <- colSums(gen.out)
table(kfp)
