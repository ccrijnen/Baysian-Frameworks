##############
#### Data ####
##############

N = 29

Y = c(1.12, 1.12, 0.99, 1.03, 0.92, 0.90, 0.81, 0.83, 0.65, 0.67, 0.60, 
      0.59, 0.51, 0.44, 0.43, 0.43, 0.33, 0.30, 0.25, 0.24, 0.13, -0.01, 
      -0.13,  -0.14, -0.30, -0.33, -0.46,  -0.43, -0.65)

x = c(-1.39, -1.39, -1.08, -1.08, -0.94, -0.80, -0.63, -0.63, -0.25, 
      -0.25, -0.12, -0.12, 0.01, 0.11, 0.11, 0.11,  0.25, 0.25, 0.34, 
      0.34, 0.44, 0.59, 0.70, 0.70, 0.85, 0.85,  0.99, 0.99, 1.19)

inits1 = list(chain1 = list(alpha = 0.2, beta = c(-0.45, -1.0), tau = 5, k = 16), 
              chain2 = list(alpha = 0.6, beta = c(-0.45, -1.0), tau = 5, k = 8))

inits2 = list(chain1 = list(alpha = 0.47, beta = c(-0.45, -1.0), tau = 5, x.change = 0.5), 
              chain2 = list(alpha = 0.47, beta = c(-0.45, -1.0), tau = 5, x.change = 0.5))




##########################################
#### Visualize Data and BUGS solution ####
##########################################

alpha.mean = 0.538
beta1.mean = -0.417
beta2.mean = -1.015
sigma.mean = 0.02214
x.change.mean = 0.02557

line1 = alpha.mean + beta1.mean * (x - x.change.mean)
line2 = alpha.mean + beta2.mean * (x - x.change.mean)

plot(x, Y, main="Stagnant Data", xlab="x", ylab="Y")
lines(x, line1, col="blue", lwd=2, lty=2)
lines(x, line2, col="red", lwd=2, lty=2)
points(c(x.change.mean), (alpha.mean), pch=23, cex=1.5, col="black", bg="purple")




plot.results = function(chain1, chain2, monitors, main) {
   par(mfrow = c(length(monitors), 1), mar = c(4.1, 4.1, 0.5, 0.5), 
       mgp = c(2.3, 1, 0), oma = c(0, 0, 3, 0))
   
   for(i in 1:length(monitors)) {
      plot(chain1[[monitors[i]]], type = "l", col = "red", xlab = "iteration", 
           ylab = monitors[i], ylim = range(c(chain1[[monitors[i]]], chain2[[monitors[i]]])))
      lines(chain2[[monitors[i]]], col = "blue")
   }
   
   title(main = main, outer = T, line = 1)
   
   par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0),
       oma = c(0, 0, 0, 0))
}




##################
#### OpenBUGS ####
##################

library(R2OpenBUGS)
data = list("N", "Y", "x")


#### Model 1

cat("model {
   	for( i in 1 : N ) {
   		Y[i] ~ dnorm(mu[i],tau)
   		mu[i] <- alpha + beta[J[i]] * (x[i] - x[k])
   		J[i] <- 1 + step(i - k - 0.5)
   		punif[i] <- 1/N
   	}
   	tau ~ dgamma(0.001,0.001)
   	alpha ~ dnorm(0.0,1.0E-6)
   	for( j in 1 : 2 ) {
   		beta[j] ~ dnorm(0.0,1.0E-6)
   	}
   	k ~ dcat(punif[])
   	sigma <- 1 / sqrt(tau)
   }", file = "stagnant_model1.txt")

par1 = c("alpha", "k")

openbugs.results1 = bugs(data = data,
                         inits = inits1,
                         parameters.to.save = par1,
                         model = "stagnant_model1.txt",
                         n.chains = 2,
                         n.iter = 10000,
                         n.thin = 1,
                         n.burnin = 0,
                         debug = FALSE,
                         DIC = TRUE)

sims = openbugs.results1$sims.array

plot.results(list("alpha" = sims[,1,"alpha"],
                  "k" = sims[,1,"k"]),
             list("alpha" = sims[,2,"alpha"],
                  "k" = sims[,2,"k"]),
             monitors = c("alpha", "k"),
             main = "OpenBUGS Model 1 simulations")

# Print summary
print(openbugs.results1)


#### Model 2

cat("model {
   	for(i in 1 : N) {
   		Y[i] ~ dnorm(mu[i], tau)
   		mu[i] <- alpha + beta[J[i]] * (x[i] - x.change)		
   		J[i] <- 1 + step(x[i] - x.change)
   	}
   	tau ~ dgamma(0.001, 0.001)
   	alpha ~ dnorm(0.0,1.0E-6)
   	for(j in 1 : 2) {
   		beta[j] ~ dnorm(0.0,1.0E-6)
   	}
   	sigma <- 1 / sqrt(tau)
   	x.change ~ dunif(-1.3,1.1)
   }", file = "stagnant_model2.txt")

par2 = c("alpha", "beta", "sigma", "x.change")

openbugs.results2 = bugs(data = data,
                         inits = inits2,
                         parameters.to.save = par2,
                         model = "stagnant_model2.txt",
                         n.chains = 2,
                         n.iter = 10000,
                         n.thin = 1,
                         n.burnin = 0,
                         debug = FALSE,
                         DIC = TRUE)

sims = openbugs.results2$sims.array

plot.results(list("alpha" = sims[,1,"alpha"],
                  "beta[1]" = sims[,1,"beta[1]"],
                  "beta[2]" = sims[,1,"beta[2]"],
                  "sigma" = sims[,1,"sigma"],
                  "x.change" = sims[,1,"x.change"]),
             list("alpha" = sims[,2,"alpha"],
                  "beta[1]" = sims[,2,"beta[1]"],
                  "beta[2]" = sims[,2,"beta[2]"],
                  "sigma" = sims[,2,"sigma"],
                  "x.change" = sims[,2,"x.change"]),
             monitors = c("alpha", "beta[1]", "beta[2]", "sigma", "x.change"),
             main = "OpenBUGS Model 2 simulations")

index = sample(1:10000, 500)
plot(sims[index, 1, "alpha"], sims[index, 1, "x.change"], xlab = "alpha", ylab = "x.change", 
     main = "alpha vs x.change")

# Print summary
print(openbugs.results2)




##############
#### JAGS ####
##############

library(R2jags)
data = list("N", "Y", "x")


#### Model 1

par1 = c("alpha", "k")

jags.result1 = jags(data = data,
                    inits = inits1,
                    parameters.to.save = par1,
                    model = "stagnant_model1.txt", 
                    n.chains = 2,
                    n.iter = 10000,  
                    n.thin = 1, 
                    n.burnin = 0,
                    DIC = TRUE)

sims = jags.result1$BUGSoutput$sims.array

plot.results(list("alpha" = sims[,1,"alpha"],
                  "k" = sims[,1,"k"]),
             list("alpha" = sims[,2,"alpha"],
                  "k" = sims[,2,"k"]),
             monitors = c("alpha", "k"),
             main = "JAGS Model 1 simulations")

# Print summary
print(jags.result1)


#### Model 2

par2 = c("alpha", "beta", "sigma", "x.change")

jags.result2 = jags(data = data,
                    inits = inits2,
                    parameters.to.save = par2,
                    model = "stagnant_model2.txt",
                    n.chains = 2,
                    n.iter = 10000,
                    n.thin = 1,
                    n.burnin = 0,
                    DIC = TRUE)

sims = jags.result2$BUGSoutput$sims.array

plot.results(list("alpha" = sims[,1,"alpha"],
                  "beta[1]" = sims[,1,"beta[1]"],
                  "beta[2]" = sims[,1,"beta[2]"],
                  "sigma" = sims[,1,"sigma"],
                  "x.change" = sims[,1,"x.change"]),
             list("alpha" = sims[,2,"alpha"],
                  "beta[1]" = sims[,2,"beta[1]"],
                  "beta[2]" = sims[,2,"beta[2]"],
                  "sigma" = sims[,2,"sigma"],
                  "x.change" = sims[,2,"x.change"]),
             monitors = c("alpha", "beta[1]", "beta[2]", "sigma", "x.change"),
             main = "JAGS Model 2 simulations")

# Print summary
print(jags.result2)




################
#### Nimble ####
################

library(nimble, warn.conflicts = FALSE)

dataNimble = list(x = x, Y = Y)
constants = list(N = N)

#### Model 1

monitors.model1 = list("alpha", "k")

# model configuration
code.model1 = nimbleCode({
   for(i in 1 : N) {
      Y[i] ~ dnorm(mu[i],tau)
      mu[i] <- alpha + beta[J[i]] * (x[i] - x[k])
      J[i] <- 1 + step(i - k - 0.5)
      punif[i] <- 1/N
   }
   tau ~ dgamma(0.001,0.001)
   alpha ~ dnorm(0.0,1.0E-6)
   for( j in 1 : 2 ) {
      beta[j] ~ dnorm(0.0,1.0E-6)
   }
   k ~ dcat(punif[1:N])
   sigma <- 1 / sqrt(tau)
})

nimble.model1 = nimbleModel(code.model1,
                            constants = constants,
                            data = dataNimble,
                            inits = inits1$chain1)

# MCMC configuration and building
mcmc.conf.model1 = configureMCMC(nimble.model1, monitors = monitors.model1)
mcmc.model1 = buildMCMC(mcmc.conf.model1, monitors = monitors.model1)

# compile to C++, and run
compiled.model1 = compileNimble(nimble.model1)
compiled.mcmc.model1 = compileNimble(mcmc.model1, project = nimble.model1)

# run multiple MCMC chains
nimble.results1 = runMCMC(compiled.mcmc.model1, 
                          nburnin = 0,
                          inits = inits1,
                          niter = 10000,
                          nchains = 2,
                          summary = TRUE)

sims = nimble.results1$samples

plot.results(chain1 = list("alpha" = sims$chain1[,"alpha"],
                           "k" = sims$chain1[,"k"]),
             chain2 = list("alpha" = sims$chain2[,"alpha"],
                           "k" = sims$chain2[,"k"]),
             monitors = c("alpha", "k"),
             main = "Nimble Model 1 simulations")

nimble.results1$summary


#### Model 2

monitors.model2 = list("alpha", "beta", "sigma", "x.change")

# model configuration
code.model2 = nimbleCode({
   for(i in 1 : N) {
      Y[i] ~ dnorm(mu[i],  tau)
      mu[i] <- alpha + beta[J[i]] * (x[i] - x.change)		
      J[i] <- 1 + step(x[i] - x.change)
   }
   tau ~ dgamma(0.001, 0.001)
   alpha ~ dnorm(0.0,1.0E-6)
   
   for(j in 1:2){
      beta[j] ~ dnorm(0.0,1.0E-6)
   }
   sigma <-  1 / sqrt(tau)
   x.change ~ dunif(-1.3,1.1)
})

nimble.model2 = nimbleModel(code.model2,
                     constants = constants,
                     data = dataNimble,
                     inits = inits2$chain1)

# MCMC configuration and building
mcmc.conf.model2 = configureMCMC(nimble.model2, monitors = monitors.model2)
mcmc.model2 = buildMCMC(mcmc.conf.model2, monitors = monitors.model2)

# compile to C++, and run
compiled.model2 = compileNimble(nimble.model2)
compiled.mcmc.model2 = compileNimble(mcmc.model2, project = nimble.model2)

# run multiple MCMC chains
nimble.results2 = runMCMC(compiled.mcmc.model2, 
                          nburnin = 0,
                          inits = inits2,
                          niter = 10000,
                          nchains = 2,
                          summary = TRUE)

sims = nimble.results2$samples

plot.results(chain1 = list("alpha" = sims$chain1[,"alpha"],
                           "beta[1]" = sims$chain1[,"beta[1]"],
                           "beta[2]" = sims$chain1[,"beta[2]"],
                           "sigma" = sims$chain1[,"sigma"],
                           "x.change" = sims$chain1[,"x.change"]),
             chain2 = list("alpha" = sims$chain2[,"alpha"],
                           "beta[1]" = sims$chain2[,"beta[1]"],
                           "beta[2]" = sims$chain2[,"beta[2]"],
                           "sigma" = sims$chain2[,"sigma"],
                           "x.change" = sims$chain2[,"x.change"]),
             monitors = c("alpha", "beta[1]", "beta[2]", "sigma", "x.change"),
             main = "Nimble Model 2 simulations")

nimble.results2$summary




##############
#### Stan ####
##############

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


#### Model 1

stan.results1 = stan('stagnant.stan', 
                     data = c("N","x","Y"),
                     chains = 4, 
                     iter = 2000)

cat(get_stancode(stan.results1))

sims = as.array(stan.results1)

# plot.results(list("alpha" = sims[,1,"alpha"],
#                   "beta[1]" = sims[,1,"beta[1]"],
#                   "beta[2]" = sims[,1,"beta[2]"],
#                   "sigma" = sims[,1,"sigma"]),
#              list("alpha" = sims[,2,"alpha"],
#                   "beta[1]" = sims[,2,"beta[1]"],
#                   "beta[2]" = sims[,2,"beta[2]"],
#                   "sigma" = sims[,2,"sigma"]),
#              monitors = c("alpha", "beta[1]", "beta[2]", "sigma"),
#              main = "Stan Model 1 results")

print(stan.results1)

traceplot(stan.results1, pars = c("alpha", "beta", "sigma"))


#### Model 2

stan.results2 = stan('stagnant2.stan', 
                     data = c("N","x","Y"), 
                     chains = 2, 
                     iter = 2000)

cat(get_stancode(stan.results2))

sims = as.array(stan.results2)

plot.results(list("alpha" = sims[,1,"alpha"],
                  "beta[1]" = sims[,1,"beta[1]"],
                  "beta[2]" = sims[,1,"beta[2]"],
                  "sigma" = sims[,1,"sigma"],
                  "x.change" = sims[,1,"x_change"]),
             list("alpha" = sims[,2,"alpha"],
                  "beta[1]" = sims[,2,"beta[1]"],
                  "beta[2]" = sims[,2,"beta[2]"],
                  "sigma" = sims[,2,"sigma"],
                  "x.change" = sims[,2,"x_change"]),
             monitors = c("alpha", "beta[1]", "beta[2]", "sigma", "x.change"),
             main = "Stan Model 2 simulations")

# Print summary
print(stan.results2)

# traceplot(stan.results2, pars = c("alpha", "beta", "sigma", "x_change"))
