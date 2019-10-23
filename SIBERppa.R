rm(list=ls())


library(SIBER)
library(ggplot2)
library(magrittr) # to enable piping
library(dplyr)
library(rjags)
library(nicheROVER)
library(mvtnorm)


setwd("C:/Users/Uli/Desktop/Análisis dieta/02. Análisis isótopos/00. SIBER/SIBER PPA")

data <- read.csv ("CNData_PPA_SIBER_clusters.csv", sep = ";")
head(data)

# create the siber object
siber.example<- createSiberObject(data)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 3, lwd = 3)
group.ellipses.args  <- list(n = 100, p.interval = NULL , lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'),
                
)


# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siber.example) 
print(community.ML)

group.ML <- groupMetricsML(siber.example)
group.ML

write.csv(group.ML)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains
# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3


ellipses.posterior <- siberMVN(siber.example, parms, priors)

SEA.B <- siberEllipses(ellipses.posterior)
SEA.B

siberDensityPlot(SEA.B, xticklabels = c("North", "South"),
                 ylims = c(0,20),
                 xlab = c(""),
                 ylab = expression("Standard Ellipse Area  " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
)
data


# Calculate some credible intervals 
cr.p <- c(0.95) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

SEA.B.credibles
# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)



#posterior comparisons

Pg1.lt.g2 <- sum( SEA.B[,2] > SEA.B[,3] ) / nrow(SEA.B)
print(Pg1.lt.g2)

#OVERLAP

# In this example, I will calculate the overlap between ellipses for groups 2
# and 3 in community 1 (i.e. the green and yellow open circles of data).

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# So in this example: community 1, group 2
ellipse1 <- "2014.2014" 

# Ellipse two is similarly defined: community 1, group3
ellipse2 <- "2016.2016"

ellipse3 <- "2017.2017"

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse2, ellipse3, siber.example, 
                                   p.interval = 0.95, n = 100)

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])

prop.95.over

# and the corresponding Bayesian estimates for the overlap between the 
# 95% ellipses is given by:



bayes95.overlap <- bayesianOverlap(ellipse2, ellipse3, ellipses.posterior,
                                   draws = 100, p.interval = 0.95, n = 100)


bayes.prop.95.over <- (bayes95.overlap[,3] / (bayes95.overlap[,2] + 
                                                bayes95.overlap[,1] -
                                                bayes95.overlap[,3])
)

bayes.prop.95.over

hist(bayes.prop.95.over, 10)

median(bayes.prop.95.over)
quantile(bayes.prop.95.over, probs = c(0.025, 0.5, 0.975))
