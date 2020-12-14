# Climate indices - this gives MEI impact on overlap

library(knitr)
library(readr)
library(tidyverse)
library(R2jags)
library(abind)
library(boot)
library(MCMCvis)
library(reshape2)
library(truncnorm)

library(lattice)
library(knitr)
library(kableExtra)

mei.data <- read.csv("mei_overlap.csv", header=F, stringsAsFactors = F)
colnames(mei.data)<- c("year","island","mei","overlap")
mei.data$island[mei.data$island=="midway"] <- 1
mei.data$island[mei.data$island=="tern"] <- 2
mei.data <- mei.data[c(2,1,3,4,5,7,6,8,9,10),]
data_by_group <- mei.data %>% group_split(island)
year=c(1,2,3,4,5,1,2,3,4,5)

jags.directory="/usr/local/bin"

JAGS<-"climateANOVA.jags"

cat("
model {
  # likelihood
    for (i in 1:10){
  		y[i] ~ dbeta(a[i], b[i])
  		a[i] <- mu[i]*phi
  		b[i] <- (1-mu[i])*phi
  		logit(mu[i]) <- alpha + beta_island[island[i]]*mei[i]+beta_year[year[i]]*mei[i]
    }
    for (j in 1:2){
  	beta_island[j] ~ dnorm(slope_hyper_island,sd_hyper_island)
  	}
  	for (j in 1:10){
  	beta_year[j] ~ dnorm(slope_hyper_year,sd_hyper_year)
  	}
	
	  # contrasts
    	for (j in 1:5){
    	  delta_year[j] <- beta_year[j]-beta_year[1]
    	}
    	for (j in 6:10){
    	  delta_year[j] <- beta_year[j]-beta_year[6]
    	}
    	for (j in 1:2){
    		delta_ternminusmid[j] <- beta_island[j]-beta_island[1]
      }

	# priors
	alpha ~ dnorm(0,0.001)
	slope_hyper_island ~ dnorm(0,0.25)
	slope_hyper_year ~ dnorm(0,0.25)
	sd_hyper_island ~ dnorm(0,0.25)T(0,)
	sd_hyper_year ~ dnorm(0,0.25)T(0,)
	phi ~ dgamma(0.1,0.1) # variance term (proportional reciprocal of the variance, called 'phi' here)
}",fill = TRUE, file=JAGS)

Dat <- list(
  y=mei.data$overlap,
  mei=mei.data$mei,
  island = mei.data$island,
  year=year
)


InitStage <- list(list(slope_hyper_year=-2,slope_hyper_island=-2,phi=0.0001, sd_hyper_island= 0.1,sd_hyper_year= 0.1), 
                  list(slope_hyper_year=0,slope_hyper_island=0,phi=1, sd_hyper_island= 1,sd_hyper_year= 1), 
                  list(slope_hyper_year=2,slope_hyper_island=2,phi=10, sd_hyper_island= 5,sd_hyper_year= 5))

# ParsStage <- c("mu","sigma", "beta_island", "slope_hyper", "alpha", "delta_ternminusmid")
ParsStage <- c("delta_year","beta_island","beta_year", "delta_ternminusmid","mu")
ni <- 120000  # number of draws from the posterior
nt <- 1   # thinning rate
nb <- 20000  # number to discard for burn-in
nc <- 3  # number of chains

results = jags(inits=InitStage,
               n.chains=nc,
               model.file=JAGS,
               working.directory=getwd(),
               data=Dat,
               parameters.to.save=ParsStage,
               n.thin=nt,
               n.iter=ni,
               n.burnin=nb,
               DIC=T)

results
pretty<- MCMCsummary(results, params="mu")
MCMCtrace(results, pdf=T)


kable(pretty,
      digits=3, booktabs=TRUE,
      caption='Summary of Posteriors for model parameters') %>%
  kable_styling(latex_options = "hold_position") 

MCMCplot(results, params="means")
MCMCtrace(results, pdf=T)
MCMCtrace(results, pdf=F, type="density", params = "delta_ternminusmid[2]")

ternminusmid_slopepparam<- results$BUGSoutput$sims.matrix[,11]
save(ternminusmid_slopepparam, file="ternminusmid_slopeparam")
hist(ternminusmid_slopepparam, breaks=200)