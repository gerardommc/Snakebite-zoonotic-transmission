library(raster)

incid.data <- readRDS("Data/Incidence models data-Apr2020.rds")

envs.r <- raster("Model-results/Spatial/Envenomings/No-envenomings.asc")
env.inc.r <- raster("Model-results/Spatial/Envenomings/Incidence.asc")

envs <- data.frame(rasterToPoints(envs.r))
env.inc <- data.frame(rasterToPoints(env.inc.r))

data.gg <- data.frame(land.cover = as.factor(incid.data[, "land.cover"]),
                      envenomings = envs$No.envenomings,
                      incidence = env.inc$Incidence)
levels(data.gg$land.cover) <- c("Forest", "Degraded", "Agricultural", "Urban", "Tea")

library(ggplot2)

pdf("Figures/Envenomings.pdf", width = 4, height = 4.5)
ggplot(data.gg) + geom_boxplot(aes(x = land.cover, y = envenomings, colour = land.cover, fill = land.cover), alpha = 0.3) +
    labs(x = "Land cover", y = "Number of envenomings") + theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.5)) 
dev.off()

pdf("Figures/Env-incidence.pdf", width = 4, height = 4.5)
ggplot(data.gg) + geom_boxplot(aes(x = land.cover, y = incidence, colour = land.cover, fill = land.cover), alpha = 0.3) +
    labs(x = "Land cover", y = "Envenoming incidence") + theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()
