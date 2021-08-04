x <- readRDS("Model-results/Envenoming-NIMBLE-model.rds")

library(coda); library(raster); library(foreach)
incid.data <- read.csv("Data/Incid-data.csv")

rho <- x[[1]][, paste0("rho[", 1:3057, "]")]

rho.r <- rasterFromXYZ(data.frame(incid.data[, c("x", "y")], colMeans(rho)))
plot(exp(rho.r))

x1 <- lapply(x, function(x1){x1[, c(paste0("inter[", 1:5, "]"), 
                                    paste0("s.effects[", 1:7, "]"), "r")]})

traceplot(as.mcmc.list(x1))
gelman.diag(as.mcmc.list(x1))

sum.inc.pars <- foreach(i = 1:13, .combine = rbind) %do% {
    med <- median(c(x1[[1]][, i], x1[[2]][, i], x1[[3]][, i]))
    std.dev <- sd(c(x1[[1]][, i], x1[[2]][, i], x1[[3]][, i]))
    cr.int <- HPDinterval(as.mcmc(x1[[1]][, i], x1[[2]][, i], x1[[3]][, i]), 0.95)
    return(c(med, std.dev, cr.int))
}

sum.inc.pars <- data.frame(sum.inc.pars)

names(sum.inc.pars) <- c("Median", "Std.Dev", "Cr.025", "Cr.975")
rownames(sum.inc.pars) <- colnames(x1[[1]])

write.csv(sum.inc.pars, "Model-results/Envenoming-estimates.csv")

x.rho <- lapply(x, function(x1){x1[, paste0("rho[", 1:3057, "]")]})

diag.rho <- geweke.diag(x.rho)

diag.rho.r <- rasterFromXYZ(data.frame(incid.data[, c("x", "y")], diag.rho$z))
plot(diag.rho.r)
##
incid.data <- read.csv("Data/Incid-data.csv")[, -1]

indices <- read.csv("Data/Agressivenes-indices.csv")

snake.pars <- read.csv("Data/Parameters-questions.csv")

rel.abund <- with(snake.pars, Density_5k)

bites.model <- readRDS("Model-results/Snakebites-mass-action-NIMBLE-model.rds")

bites.hum.pars <- bites.model[[1]]

beta0 <- apply(bites.hum.pars[, paste0("beta0[", 1:5, "]")], 2, function(x){
    coda::HPDinterval(coda::as.mcmc(x), 0.5)
})[1,]
beta1 <- apply(bites.hum.pars[, paste0("beta1[", 1:5, "]")], 2, function(x){
    coda::HPDinterval(coda::as.mcmc(x), 0.5)
})[1, ]

model.bites <- raster("Model-results/Spatial/Snakebites/Incidence.asc")

Hum.log <- log(incid.data$hum.pop)
bit.inc <- extract(model.bites, incid.data[, c("x", "y")])
env.inc <- incid.data$envenomings
env.inc <- ifelse(env.inc > bit.inc, bit.inc, env.inc)

Hum.env <- ceiling(env.inc * incid.data$hum.pop)
Hum.bit <- ceiling(bit.inc * incid.data$hum.pop)

Beta.i <- exp(beta0[incid.data$land.cover] + beta1[incid.data$land.cover] * Hum.log^2)

S <- as.matrix(incid.data[, 3:9])
max.s <- apply(S, 2, max)
S <- sweep(S, 2, max.s, "/")
S <- sweep(S, 2, indices$Severity/10 * rel.abund, "*")
S <- sweep(S, 1, Beta.i, "*")

library(data.table)

x.full <- rbind(x[[1]], x[[2]], x[[3]])

x.int <- x.full[, paste0("inter[", 1:5, "]")]
x.s <- x.full[, paste0("s.effects[", 1:7, "]")]
x.rho <- x.full[, paste0("rho[", 1:3057, "]")]

s.preds <- apply(x.s, 1, function(x){
    rowSums(sweep(S, MARGIN = 2, x, FUN = "*"))
})

s.effs <- s.preds + t(x.int[, incid.data$land.cover])

inc.preds <- exp(s.effs)/(1+exp(s.effs)) * exp(t(x.rho))

inc.sum <- data.frame(t(apply(inc.preds, 1, function(x1){HPDinterval(as.mcmc(x1), 0.5)})))
names(inc.sum) <- c("incid", "incid.1")

prob.env <- rasterFromXYZ(data.frame(incid.data[, c("x", "y")], inc.sum$incid))
env.inc.mod <- prob.env*model.bites

env.inc <- rasterFromXYZ(incid.data[, c("x", "y", "envenomings")])

plot(stack(env.inc, env.inc.mod))
pairs(stack(env.inc, env.inc.mod))

hum.pop.r <- rasterFromXYZ(incid.data[, c("x", "y", "hum.pop")])
no.envs <- hum.pop.r * env.inc.mod
plot(no.envs)

rho.med <- rasterFromXYZ(data.frame(incid.data[c("x", "y")], apply(x.rho, 2, median)))
plot(rho.med)

writeRaster(env.inc.mod, "Model-results/Spatial/Envenomings/Incidence", "ascii", overwrite = T)
writeRaster(no.envs, "Model-results/Spatial/Envenomings/No-envenomings", "ascii", overwrite = T)
writeRaster(rho.med, "Model-results/Spatial/Envenomings/Rho", "ascii", overwrite = T)
