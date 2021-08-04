x <- readRDS("Snakebites-mass-action-NIMBLE-model.rds")

library(coda); library(raster); library(foreach)

incid.data <- read.csv("Incid-data.csv")

rho <- x[[1]][,23:3079]

rho.r <- rasterFromXYZ(data.frame(incid.data[, c("x", "y")], colMeans(rho)))
plot(rho.r)

x1 <- lapply(x, function(x1){x1[, 1:22]})

traceplot(as.mcmc.list(x1))
gelman.diag(as.mcmc.list(x1))

sum.inc.pars <- foreach(i = 1:22, .combine = rbind) %do% {
    med <- median(c(x1[[1]][, i], x1[[2]][, i], x1[[3]][, i]))
    std.dev <- sd(c(x1[[1]][, i], x1[[2]][, i], x1[[3]][, i]))
    cr.int <- HPDinterval(as.mcmc(x1[[1]][, i], x1[[2]][, i], x1[[3]][, i]), 0.95)
    return(c(med, std.dev, cr.int))
}

sum.inc.pars <- data.frame(sum.inc.pars)

names(sum.inc.pars) <- c("Median", "Std.Dev", "Cr.025", "Cr.975")
rownames(sum.inc.pars) <- colnames(x1[[1]])

write.csv(sum.inc.pars, "Mass-action-estimates.csv")

x.rho <- lapply(x, function(x1){x1[, 23:3079]})

diag.rho <- geweke.diag(x.rho)

diag.rho.r <- rasterFromXYZ(data.frame(incid.data[, c("x", "y")], diag.rho$z))
plot(diag.rho.r)
##
incid.data <- read.csv("Incid-data.csv")[,-1]

indices <- read.csv("Agressivenes-indices.csv")

snake.pars <- read.csv("Parameters-questions.csv")

rel.abund <- with(snake.pars, Density_5k)

S <- as.matrix(incid.data[, 3:9])
max.s <- apply(S, 2, max)
S <- sweep(S, 2, max.s, "/")
S <- sweep(S, 2, indices$Agressiveness/10 * rel.abund, "*")

Hum.log <- log(incid.data$hum.pop)
Hum.bit <- sapply(incid.data$bites * incid.data$hum.pop, function(x){qpois(0.5, x)})

library(data.table)

x.full <- x[[1]]

x.beta0 <- x.full[, paste0("beta0[", 1:5, "]")]
x.beta1 <- x.full[, paste0("beta1[", 1:5, "]")]
x.indices <- x.full[, paste0("indices[", 1:7, "]")]
x.rho <- x.full[, paste0("rho[", 1:3057, "]")]

s.preds <- apply(x.indices, 1, function(x){
    rowSums(sweep(S, MARGIN = 2, x, FUN = "*"))
})

library(foreach)
hu.ef <- foreach(i = 1:1000, .combine = cbind) %do%{
        exp(x.beta0[i, incid.data$land.cover] + x.beta1[i, incid.data$land.cover] * Hum.log^2)
}
    
inc.preds <- (1 - exp(-s.preds * hu.ef)) * exp(t(x.rho))

inc.sum <- data.frame(t(apply(inc.preds, 1, function(x){HPDinterval(as.mcmc(x), prob = 0.5)})))
names(inc.sum) <- c("id", "incid")

inc.r <- rasterFromXYZ(data.frame(incid.data[, c("x", "y")], inc.sum$incid))

plot(inc.r)

writeRaster(inc.r, "Incidence-NIMBLE-model-SI.asc", overwrite = T)