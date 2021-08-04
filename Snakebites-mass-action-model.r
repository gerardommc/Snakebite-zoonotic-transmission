library(R2jags)

incid.data <- read.csv("Incid-data.csv")[,-1]
incid.data$land.cover <- incid.data$land.cover

indices <- read.csv("Agressivenes-indices.csv")

snake.pars <- read.csv("Parameters-questions.csv")

rel.abund <- with(snake.pars, Density_5k)

modelString <- "model{
      #Likelihood
      for(i in 1:n){
            Beta[i] <- exp(beta0[class[i]] + beta1[class[i]] * Hum.log[i]^2)
            S.sum[i] <- inprod(indices[], S[i,])
            P[i] <- 1 - exp(-Beta[i] * S.sum[i])
            H.bit[i] <- P[i] * Hum.tot[i]
            P.r[i] <- r[class[i]]/(r[class[i]] + H.bit[i])
            Hum.bit[i] ~ dnegbin(P.r[i], r[class[i]])
      }
      #Priors
      for(j in 1:n.classes){
            beta0[j] ~ dnorm(0, 1.0E-3)
            beta1[j] ~ dnorm(0, 1.0E-3)
            r[j] ~ dunif(0, 100)
      }
      for(k in 1:n.species){
            indices[k] ~ dunif(0,10)
      }
}"

writeLines(modelString, "Snakebites-simple-mass-action.jags")

S <- as.matrix(incid.data[, 3:9])
max.s <- apply(S, 2, max)
S <- sweep(S, 2, max.s, "/")
S <- sweep(S, 2, indices$Agressiveness/10 * rel.abund, "*")

Hum.log <- log(incid.data$hum.pop)
Hum.bit <- sapply(incid.data$bites * incid.data$hum.pop, function(x){qpois(0.5, x)})

jags.data <- list(
      class = incid.data$land.cover,
      n.classes = 5,
      S = S,
      n.species = 7,
      Hum.log = Hum.log,
      Hum.tot = incid.data$hum.pop,
      Hum.bit = Hum.bit,
      n = nrow(incid.data)
)

params <- c("beta0", "beta1",
            "qh", "r", "indices")

model <- jags.parallel(data = jags.data,
                       parameters.to.save = params,
                       model.file = "Snakebites-simple-mass-action.jags",
                       n.chains = 2,
                       n.burnin = 100000,
                       n.iter = 1000000,
                       n.thin = 900,
                       n.cluster = 2)

saveRDS(model, "Snakebites-mass-action-JAGS-model.rds")
