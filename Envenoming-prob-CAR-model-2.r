library(nimble); library(raster); library(spdep)

modelCode <- nimbleCode({
      #Likelihood
      for(i in 1:n){
            logit(P[i]) <- inter[class[i]] + inprod(s.effects[], S[i,])
            log(env[i]) <- log(P[i] * Hum.bit[i]) + rho[i]
            P.r[i] <- r/(r + env[i])
            Hum.env[i] ~ dnegbin(P.r[i], r)
      }

      #Priors
      rho[1:n] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n], tau, zero_mean = 0)
      sigma ~ dunif(0, 10)   # prior for variance components
      tau <- 1 / sigma^2
      for(j in 1:n.species){
            s.effects[j] ~ dnorm(E.s[j], sd = sd.s[j])
      }

      for(k in 1:n.classes){inter[k] ~ dnorm(E.int[k], sd = sd.int[k])}

      r ~ dunif(0, 50)
})

incid.data <- read.csv("Data/Incid-data.csv")[, -1]

indices <- read.csv("Data/Agressivenes-indices.csv")

snake.pars <- read.csv("Data/Parameters-questions.csv")

rel.abund <- with(snake.pars, Density_5k)


model.bites <- raster("Model-results/Spatial/Snakebites/Incidence.asc")

Hum.log <- log(incid.data$hum.pop)
bit.inc <- extract(model.bites, incid.data[, c("x", "y")])
env.inc <- incid.data$envenomings
env.inc <- ifelse(env.inc > bit.inc, bit.inc, env.inc)

Hum.env <- ceiling(env.inc * incid.data$hum.pop)
Hum.bit <- ceiling(bit.inc * incid.data$hum.pop)

S <- as.matrix(incid.data[, 3:9])
max.s <- apply(S, 2, max)
S <- sweep(S, 2, max.s, "/")
S <- sweep(S, 2, rel.abund, "*")
S.tot <- rowSums(S)
S.prop <- sweep(S, 1, S.tot, "/")

# Spatial neighbourhood

W.nb <- dnearneigh(as.matrix(incid.data[, c("x", "y")]), d1 = 6000, d2 = 10000, longlat = F)

## Determine neighborhood/adjacency information needed for neighborhood-based CAR model
nbInfo <- nb2WB(W.nb)

data <- list(Hum.env = Hum.env)

E.s <- c(280.14, -8.77, -97.34, 5.49, -24.55, -290.34, -2.88)
sd.s <- c(13.99, 32.12, 28.4, 31.91, 23.74, 29.84, 31.8)
E.int <- c(0.255, -0.39, 0.056, -0.72, 0.59)
sd.int <- c(0.07, 0.05, 0.07, 0.04, 0.17)

constants <- list(
      S = S.prop,
      n.species = 7,
      Hum.bit = Hum.bit,
      n = length(Hum.env),
      n.classes = 5,
      L = length(nbInfo$adj),
      adj = nbInfo$adj,
      weights = nbInfo$weights,
      num = nbInfo$num,
      class = incid.data$land.cover,
      E.s = E.s, sd.s = sd.s,
      E.int = E.int, sd.int = sd.int
)

inits <- list(s.effects = E.s, inter = E.int, rho = rnorm(nrow(incid.data)))

model <- nimbleModel(modelCode, constants = constants, data = data,
                      inits = inits)

cmodel <- compileNimble(model)

conf.model <- configureMCMC(model, monitors = c("s.effects", "inter", "r", "rho"))

MCMC <- buildMCMC(conf.model)

cMCMC <- compileNimble(MCMC, project = cmodel)

set.seed(4687)

library(doParallel)
registerDoParallel(cores = 3)

samples <- foreach(i = 1:3) %dopar% {
  runMCMC(cMCMC, niter = 250000, nburnin = 25000, thin = 225,
    samplesAsCodaMCMC = T)
}

samples <- coda::as.mcmc.list(samples)

saveRDS(samples, "Model-results/Envenoming-NIMBLE-model.rds")
