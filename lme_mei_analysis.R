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
		y[i] ~ dnorm(mu[i], tau)T(0,)
		mu[i] <- alpha + alpha_year[year[i]]+alpha_island[island[i]]+beta_island[island[i]]*mei[i]+beta_year[year[i]]*mei[i]
  }
  for (j in 1:2){
		delta_ternminusmid[j] <- beta_island[j]-beta_island[1]
  }
	for (j in 1:5){
	  delta_year[j] <- beta_year[j]-beta_year[1]
	}
	for (j in 6:10){
	  delta_year[j] <- beta_year[j]-beta_year[6]
	}

	# priors
	for (j in 1:2){
	alpha_island[j] ~ dnorm(0,0.001)
	beta_island[j] ~ dnorm(slope_hyper_island,0.001)
	}
	for (j in 1:10){
	alpha_year[j] ~ dnorm(0,0.001)
	beta_year[j] ~ dnorm(slope_hyper_year,0.001)
	}
	alpha ~ dnorm(0,0.001)
	slope_hyper_island ~ dnorm(0,0.001)
	slope_hyper_year ~ dnorm(0,0.001)
	sigma ~ dunif(0,100)
	tau <- pow(sigma, -2)
}",fill = TRUE, file=JAGS)

Dat <- list(
  y=mei.data$overlap,
  mei=mei.data$mei,
  island = mei.data$island,
  year=year
)


InitStage <- list(list(slope_hyper_year=-1,sigma=1), list(slope_hyper_year=0,sigma=10), list(slope_hyper_year=1,sigma=20))

# ParsStage <- c("mu","sigma", "beta_island", "slope_hyper", "alpha_island", "delta_ternminusmid", "delta_midway_year")
ParsStage <- c("delta_year","beta_island","beta_year", "mu")
ni <- 11000  # number of draws from the posterior
nt <- 1    # thinning rate
nb <- 1000  # number to discard for burn-in
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
pretty<- MCMCsummary(results)
MCMCtrace(results, pdf=T)


kable(pretty,
      digits=3, booktabs=TRUE,
      caption='Summary of Posteriors for model parameters') %>%
  kable_styling(latex_options = "hold_position") %>%
  footnote(general='Estimates based on final 5,000 samples.')

MCMCplot(results, params="means")
MCMCtrace(results, pdf=T)
MCMCtrace(results, pdf=F, type="density", params = "delta_ternminusmid[2]")

ternminusmid_slopepparam<- results$BUGSoutput$sims.matrix[,11]
save(ternminusmid_slopepparam, file="ternminusmid_slopeparam")
hist(ternminusmid_slopepparam, breaks=200)


# maybe helpful: 

# for footnote
# kable(pretty,
#       digits=3, booktabs=TRUE,
#       caption='Summary of Posteriors for model parameters') %>%
#   kable_styling(latex_options = "hold_position") %>%
#   footnote(general='Estimates based on final 5,000 samples.')
# files <- paste0("simulated_path",1:50,".txt") 
# for( i in 1:length(files))  {
#   data <- all.paths[[i]]
#   filename <-  sub(".txt",".jpg",files[i])
#   jpeg(file=filename)
#   plot(data$x, data$y, type="l", asp =1)
#   points(data$x[1],data$y[1], col="green")
#   points(data$x[165],data$y[165], col="red")                         
#   dev.off()
# }
