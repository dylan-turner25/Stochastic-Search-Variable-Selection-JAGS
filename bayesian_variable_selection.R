rm(list = ls()) # clear console
options(scipen = 999) # forces R to avoid exponential notation
system.info <- Sys.info()
setwd("/home/dylan/Dropbox/Research/Sustainability Special Issue")


# load libraries 
library(runjags) # need to install JAGS first, which is independent of R
library(tidyverse)

# load cleaned data
data <- read_csv("./Data/cleaned_data_09_20_20.csv") %>% drop_na()
data$risk <- exp(data$ln_risk)
data <- data[,-which(colnames(data) == "ln_risk")]

regressors <- data[,-c(1,2)] # use all regressors
regressors <- data[,c("prem_fc","hhinc","probcat31","expdamage","risk","sfha",
                      "flood_exp","newtocoast","stormrisk","retreat","armor","affect_homeloss",
                      "eloan","eia","hhsize","age","higheredu")] # use some regressors

# put data into jags format
Y <- data$floodins
n <- length(Y)  
X <- as.matrix(regressors)
b <- ncol(X) # number of parameters under consideration

#put data into list format that JAGS likes
data.list <- list(Y = Y, n = n, X = X, b = b)
data.runjags <- dump.format(data.list)

#define the model string 

model <- "model{
    for(i in 1:n){
  Y[i] ~ dbern(p[i])
  probit(p[i]) <- mu[i]
  mu[i] <- cons + inprod(X[i,],delta)

}

for (j in 1:b) {
ind[j]~dbern(pind)
deltaT[j]~dnorm(0,taub)
delta[j]<-ind[j]*deltaT[j]
}
cons~dnorm(0,0.0001)
tau~dgamma(1,1)
taub~dgamma(1,1)
pind~dbeta(2,8)

}"

monitor <- c("cons","delta","ind","tau","taub","pind")
inits <- list()

# set initial values for each chain
chains <- 3
for (m in 1:chains) {
  assign(paste("inits", m , sep = ""), list(tau = 1,taub = 1,cons = .6 ,deltaT=rep(0,b),ind=rep(0,b),.RNG.name = "base::Wichmann-Hill", .RNG.seed = sample(1:1000,1))) 
  inits[[m]] <- eval(parse(text = paste("inits",m, sep = "")))      
}


# run the model and assign it to a numbered object "jags.output"
jags.output <- run.jags(model = model, method = "parallel",monitor = monitor , 
                        thin = 10, inits = inits, data = data.runjags , adapt = 10000, 
                        n.chains = chains, burnin = 10000, sample = 1000, keep.jags.files = F)

jags.df <- data.frame(summary(jags.output))
jags.df$var_names <- c("cons",colnames(regressors),colnames(regressors),"tau","taub","pind" )
jags.df <- jags.df[,c(ncol(jags.df),1:ncol(jags.df)-1)]

while (1 == 1 ) {
  # save image in case of a crash
  save.image("./Bayesian Variable Selection/Some Regressors/SSVS_Results.Rdata")
  
  # save intermediate results to a csv file
  write.csv(jags.df,paste("./Bayesian Variable Selection/Some Regressors/SSVS_Results_",jags.output$sample,".csv",sep=""), row.names = T)
  
  # extend the chain
  jags.output <- extend.jags(jags.output, sample = 10000, adapt = 10000, thin = 10)
  
  # put output in data.frame
  jags.df <- summary(jags.output)
  jags.df$var_names <- c("cons",colnames(regressors),colnames(regressors),"tau","taub","pind" )
  jags.df <- jags.df[,c(ncol(jags.df),1:ncol(jags.df)-1)]
}



