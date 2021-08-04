library(nimble); library(spdep); library(raster)

modelCode <- nimbleCode({
    for(i in 1:n){
        Beta[i] <- exp(beta0[class[i]] + beta1[class[i]] * Hum.log[i]^2)
        P[i] <- 1 - exp(-Beta[i] * inprod(indices[], S[i,]))
        log(H.bit[i]) <- log(P[i] * Hum.tot[i]) + rho[i]
        P.r[i] <- r[class[i]]/(r[class[i]] + H.bit[i])
        Hum.bit[i] ~ dnegbin(P.r[i], r[class[i]])
    }

    #Priors
    rho[1:n] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n], tau, zero_mean = 0)
    sigma ~ dunif(0, 10)   # prior for variance components
    tau <- 1 / sigma^2
    for(j in 1:n.classes){
        beta0[j] ~ dnorm(E.b0[j], sd = sd.b0[j])
        beta1[j] ~ dnorm(E.b1[j], sd = sd.b1[j])
        r[j] ~ T(dnorm(E.r[j], sd = sd.r[j],), min = 0, max = 50)
    }
    for(k in 1:n.species){
        indices[k] ~ T(dnorm(E.ind[k], sd = sd.ind[k]), min = 0,  max =10)
    }
})

#Reading and formatting data

incid.data <- read.csv("Data/Incid-data.csv")[,-1]

indices <- read.csv("Data/Agressivenes-indices.csv")

snake.pars <- read.csv("Data/Parameters-questions.csv")

rel.abund <- with(snake.pars, Density_5k)

S <- as.matrix(incid.data[, 3:9])
max.s <- apply(S, 2, max)
S <- sweep(S, 2, max.s, "/")
S <- sweep(S, 2, indices$Agressiveness/10 * rel.abund, "*")

Hum.log <- log(incid.data$hum.pop)
Hum.bit <- sapply(incid.data$bites * incid.data$hum.pop, function(x){qpois(0.5, x)})

# Spatial neighbourhood

W.nb <- dnearneigh(as.matrix(incid.data[, c("x", "y")]), d1 = 6000, d2 = 10000, longlat = F)

## Determine neighborhood/adjacency information needed for neighborhood-based CAR model
nbInfo <- nb2WB(W.nb)

data <- list(Hum.bit = Hum.bit)

E.b0 <- c(-12.9, -12.84, -13.18, -11.91, -11.27); sd.b0 <- c(0.19, 0.2, 0.21, 0.26, 0.61)
E.b1 <- c(-0.004, -0.005, 0.0, -0.017, -0.028); sd.b1 <- c(0.001, 0.001, 0.002, 0.002, 0.002)
E.r <- c(23.70, 16.82, 9.32, 7.18, 12.10); sd.r <- c(2.95, 1.13, 1.16, 0.69, 12.10)

E.ind <- c(3.48, 5.3, 1.17, 2.99, 0.036, 8.48, 0.17)
sd.ind <- c(0.56, 2.85, 1.17, 1.07, 0.02, 1.29, 0.16)

rho.r <- raster("Data/Rho-mean.asc")
rho.df <- data.frame(rasterToPoints(rho.r))
E.rho <- round(rho.df$Rho.mean, 2)

constants <- list(class = incid.data$land.cover,
                  n.classes = 5,
                  S = S,
                  n.species = 7,
                  Hum.log = Hum.log,
                  Hum.tot = incid.data$hum.pop,
                  n = nrow(incid.data),
                  L = length(nbInfo$adj),
                  adj = nbInfo$adj,
                  weights = nbInfo$weights,
                  num = nbInfo$num,
                  E.b0 = E.b0, E.b1 = E.b1,
                  sd.b0 = sd.b0, sd.b1 = sd.b1,
                  E.r = E.r, sd.r = sd.r,
                  E.ind = E.ind, sd.ind = sd.ind)

inits <- list(beta0 = E.b0, beta1 = E.b1,
              r = E.r, indices = E.ind,
              rho = E.rho)

model <- nimbleModel(modelCode, constants = constants, data = data, inits = inits)

cmodel <- compileNimble(model)

conf.model <- configureMCMC(model, monitors = c("beta0", "beta1", "r", "indices", "tau", "rho"))

MCMC <- buildMCMC(conf.model)

cMCMC <- compileNimble(MCMC, project = cmodel)

set.seed(4687)
library(doParallel)
registerDoParallel(cores = 3)
samples <- foreach(i = 1:3) %dopar% {
    runMCMC(cMCMC, niter = 1000000, nburnin = 100000, thin = 900, samplesAsCodaMCMC = T)
}
samples <- coda::as.mcmc.list(samples)

saveRDS(samples, "Model-results/Snakebites-mass-action-NIMBLE-model.rds")
