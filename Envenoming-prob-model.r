library(R2jags); library(raster); library(parallel)

incid.data <- read.csv("Incid-data.csv")

indices <- read.csv("Agressivenes-indices.csv")

snake.pars <- read.csv("Parameters-questions.csv")

rel.abund <- with(snake.pars, Density_5k)

bites.model <- readRDS("Snakebites-mass-action-JAGS-model.rds")

bites.hum.pars <- bites.model$BUGSoutput$summary

beta0 <- bites.hum.pars[paste0("beta0[", 1:5, "]"), "50%"]
beta1 <- bites.hum.pars[paste0("beta1[", 1:5, "]"), "50%"]

model.bites <- raster("Incidence-model-SI.asc")

Hum.log <- log(incid.data$hum.pop)
bit.inc <- extract(model.bites, incid.data[, c("x", "y")])
env.inc <- incid.data$envenomings
env.inc <- ifelse(env.inc > bit.inc, bit.inc, env.inc)

ones <- which((env.inc/bit.inc) == 1)

Hum.env <- ceiling(env.inc[-ones] * incid.data$hum.pop[-ones])
Hum.bit <- ceiling(bit.inc[-ones] * incid.data$hum.pop[-ones])

Beta.i <- exp(beta0[incid.data$land.cover] + beta1[incid.data$land.cover] * Hum.log^2)[-ones]

S <- as.matrix(incid.data[, 3:9])
max.s <- apply(S, 2, max)
S <- sweep(S, 2, max.s, "/")
S <- sweep(S, 2, indices$Severity/10 * rel.abund, "*")
S <- S[-ones,]
S <- sweep(S, 1, Beta.i, "*")

##Try binomial regression

modelString <- "
model{
      #Likelihood
      for(i in 1:n){
            logit(P[i]) <- int[class[i]] + inprod(s.effects[], S[i,])
            env[i] <- P[i] * Hum.bit[i]
            P.r[i] <- r/(r + env[i])
            Hum.env[i] ~ dnegbin(P.r[i], r)
      }

      #Priors
      for(j in 1:n.species){
            s.effects[j] ~ dnorm(0, 1.0E-3)
      }

      for(k in 1:n.classes){int[k] ~ dnorm(0, 1.0E-3)}

      r ~ dunif(0, 50)
}
"

writeLines(modelString, "Envenoming-model.jags")



jags.data <- list(
      S = S,
      n.species = 7,
      Hum.bit = Hum.bit,
      Hum.env = Hum.env,
      n = length(Hum.env),
      n.classes = 5,
      class = incid.data$land.cover[-ones]
)

params <- c("s.effects", "r", "int")

model <- jags.parallel(data = jags.data,
                       parameters.to.save = params,
                       model.file = "Envenoming-model.jags",
                       n.chains = 2,
                       n.burnin = 10000,
                       n.iter = 100000,
                       n.thin = 90,
                       n.cluster = 2)

saveRDS(model, "Envenoming-JAGS-model.rds")
