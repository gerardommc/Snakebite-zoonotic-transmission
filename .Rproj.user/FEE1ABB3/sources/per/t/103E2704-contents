library(raster); library(ggplot2); library(rgdal); library(SpatialPack);
library(cowplot)

incid.data <- read.csv("Data/Incid-data.csv")[, -1]

incid.data$No.bites <- with(incid.data, bites * hum.pop)
incid.data$No.envs <- with(incid.data, envenomings * hum.pop)

inc.data <- rasterFromXYZ(incid.data[, c("x", "y", "bites")])
bit.data <- rasterFromXYZ(incid.data[, c("x", "y", "No.bites")])


##Snakebite plots
bites <- rasterToPoints(raster("Model-results/Spatial/Snakebites/No-bites.asc"))
incid <- rasterToPoints(raster("Model-results/Spatial/Snakebites/Incidence.asc"))

sl <- readOGR("Data/Sri Lanka boundaries/LKA_adm0.shp")
sl <- spTransform(sl, CRSobj = CRS("+init=epsg:5235"))


bites.df <- data.frame(bites)
bites.df$bites.data <- incid.data$No.bites
names(bites.df) <- c("x", "y", "model", "data")

incid.df <- data.frame(incid)
incid.df$inc.data <- incid.data$bites
names(incid.df) <- c("x", "y", "data", "model")

library(RColorBrewer)
cool.cols <- colorRampPalette(c("midnightblue",
                                "turquoise3",
                                "violetred3",
                                "goldenrod1",
                                "grey95"))
#Bites
bites <- ggplot() + geom_tile(data = bites.df, 
                     aes(x = x, y = y, fill = model)) +
      scale_fill_gradientn(colours = rev(cool.cols(100)),
                           breaks = c(10, 160, 300)) +
      geom_polygon(data = sl,
                   aes(x=long, y=lat, group = group),
                   fill=NA, color="grey50", size=1.25) +
      coord_fixed(ratio = 1) +
      labs(fill = "", title = "Number of snakebites")+
      theme(panel.background =  element_rect(colour = NA, fill = NA),
            legend.position = c(0.95, 0.65),
            legend.background = element_rect("transparent"),
            legend.text = element_text(size = 12),
            legend.key.height = unit(5, units = "mm"),
            legend.key.width = unit(15, units = "mm"),
            plot.title = element_text(face = "bold",
                                 size = 22,
                                 hjust = 0.5),
            axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(face = "plain", size = 12))

cor.bit <- modified.ttest(x = bites.df$data, y = bites.df$model,
                      coords = bites.df[, c("x", "y")])

cor.bites <- ggplot(bites.df) + geom_hex(aes(x = data, y = model)) +
      scale_fill_gradientn(colours = rev(cool.cols(100)), trans = "log",
                           labels = c(1, 7, 50, 400)) + 
   xlim(0, 350) + ylim(0, 350) + coord_fixed(ratio = 1) + 
      geom_abline(slope = 1, intercept = 0, colour = "red", size = 1.5, alpha = 0.5) +
      geom_text(aes(x = 200, y = 100, label= paste0(
            "r = ", round(cor.bit$corr,2),
            "\n d.f. = ", round(cor.bit$dof, 0),
            "\n P = ", round(cor.bit$p.value, 2)
      )), fontface = "plain", size = 6) + 
      labs(x = "Data", y = "Model") +
      theme(
            legend.position = "none",
            panel.background = element_rect(colour = "grey40", fill = NA),
            rect = element_rect(fill = "transparent"),
            axis.title = element_text(size = 18)
      ) 

bites.inset = ggdraw() +
      draw_plot(bites) +
      draw_plot(cor.bites, x = 0.6, y = 0.7, width = 0.3, height = 0.3)

#####
#Incidence
incidence <- ggplot() + geom_tile(data = incid.df, 
                     aes(x = x, y = y, fill = model)) +
      scale_fill_gradientn(colours = rev(cool.cols(100))) +
      geom_polygon(data = sl,
                   aes(x=long, y=lat, group = group),
                   fill=NA, color="grey50", size=1.25) +
      coord_fixed(ratio = 1) + 
      labs(fill = "", title = "Snakebite incidence")+
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.95, 0.65),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(15, units = "mm"),
         plot.title = element_text(face = "bold",
                                   size = 22,
                                   hjust = 0.5),
         axis.text = element_text(size = 12),
         axis.text.x = element_text(angle = 45, hjust = 1),
         axis.title = element_text(face = "plain", size = 12))

cor.inc <- modified.ttest(x = incid.df$data, y = incid.df$model,
                          coords = incid.df[, c("x", "y")])

cor.incidence <- ggplot(incid.df) + geom_hex(aes(x = data, y = model)) +
      scale_fill_gradientn(colours = rev(cool.cols(100)), trans = "log",
                           labels = c(1, 3, 7, 20, 1)) + 
      xlim(0, 0.01) + ylim(0, 0.01) + coord_fixed(ratio = 1) + 
      geom_abline(slope = 1, intercept = 0, colour = "red", size = 1.5, alpha = 0.5) +
      geom_text(aes(x = 0.0035, y = 0.0065, label= paste0(
            "r = ", round(cor.inc$corr,2),
            "\n d.f. = ", round(cor.inc$dof, 0),
            "\n P = ", round(cor.inc$p.value, 2)
      )), fontface = "plain", size = 6) +
      labs(x = "Data", y = "Model") +
      theme(
            legend.position = "none",
            panel.background = element_rect(colour = "grey40", fill = "white"),
            rect = element_rect(fill = "transparent"),
            axis.title = element_text(size = 18)
      ) 

bit.incid.inset = ggdraw() +
      draw_plot(incidence) +
      draw_plot(cor.incidence, x = 0.6, y = 0.7, width = 0.3, height = 0.3)


################################
##Envenoming

envs <- rasterToPoints(raster("Model-results/Spatial/Envenomings/No-envenomings.asc"))
env.inc <- rasterToPoints(raster("Model-results/Spatial/Envenomings/Incidence.asc"))

envs.df <- data.frame(envs)
envs.df$data <- incid.data$No.envs
names(envs.df) <- c("x", "y", "model", "data")

env.inc.df <- data.frame(env.inc)
env.inc.df$data <- incid.data$envenomings
names(env.inc.df) <- c("x", "y", "model", "data")

#Number of envenomings
envs <- ggplot() + geom_tile(data = envs.df, 
                              aes(x = x, y = y, fill = model)) +
   scale_fill_gradientn(colours = rev(cool.cols(100)),
                        breaks = c(10, 75, 150)) +
   geom_polygon(data = sl,
                aes(x=long, y=lat, group = group),
                fill=NA, color="grey50", size=1.25) +
   coord_fixed(ratio = 1) +
   labs(fill = "", title = "Number of envenomings")+
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.95, 0.65),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(15, units = "mm"),
         plot.title = element_text(face = "bold",
                                   size = 22,
                                   hjust = 0.5),
         axis.text = element_text(size = 12),
         axis.text.x = element_text(angle = 45, hjust = 1),
         axis.title = element_text(face = "plain", size = 12))

cor.env <- modified.ttest(x = envs.df$data, y = envs.df$model,
                          coords = envs.df[, c("x", "y")])

cor.envs <- ggplot(envs.df) + geom_hex(aes(x = data, y = model)) +
   scale_fill_gradientn(colours = rev(cool.cols(100)), trans = "log") + 
   xlim(0, 160) + ylim(0, 160) + coord_fixed(ratio = 1) + 
   geom_abline(slope = 1, intercept = 0, colour = "red", size = 1.5, alpha = 0.5) +
   geom_text(aes(x = 45, y = 100, label= paste0(
      "r = ", round(cor.env$corr,2),
      "\n d.f. = ", round(cor.env$dof, 0),
      "\n P = ", round(cor.env$p.value, 2)
   )), fontface = "plain", size = 6) + 
   labs(x = "Data", y = "Model") +
   theme(
      legend.position = "none",
      panel.background = element_rect(colour = "grey40", fill = NA),
      rect = element_rect(fill = "transparent"),
      axis.title = element_text(size = 18)
   ) 

envs.inset = ggdraw() +
   draw_plot(envs) +
   draw_plot(cor.envs, x = 0.6, y = 0.7, width = 0.3, height = 0.3)


#Envenoming incidence

env.inc <- ggplot() + geom_tile(data = env.inc.df, 
                                  aes(x = x, y = y, fill = model)) +
   scale_fill_gradientn(colours = rev(cool.cols(100))) +
   geom_polygon(data = sl,
                aes(x=long, y=lat, group = group),
                fill=NA, color="grey50", size=1.25) +
   coord_fixed(ratio = 1) + 
   labs(fill = "", title = "Envenoming incidence")+
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.95, 0.65),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(15, units = "mm"),
         plot.title = element_text(face = "bold",
                                   size = 22,
                                   hjust = 0.5),
         axis.text = element_text(size = 12),
         axis.text.x = element_text(angle = 45, hjust = 1),
         axis.title = element_text(face = "plain", size = 12))

cor.env.inc <- modified.ttest(x = env.inc.df$data, y = env.inc.df$model,
                          coords = env.inc.df[, c("x", "y")])

cor.inc.env <- ggplot(env.inc.df) + geom_hex(aes(x = data, y = model)) +
   scale_fill_gradientn(colours = rev(cool.cols(100)), trans = "log") + 
   xlim(0, 0.0075) + ylim(0, 0.0075) + coord_fixed(ratio = 1) + 
   geom_abline(slope = 1, intercept = 0, colour = "red", size = 1.5, alpha = 0.5) +
   geom_text(aes(x = 0.0025, y = 0.005, label= paste0(
      "r = ", round(cor.env.inc$corr,2),
      "\n d.f. = ", round(cor.env.inc$dof, 0),
      "\n P = ", round(cor.env.inc$p.value, 2)
   )), fontface = "plain", size = 6) +
   labs(x = "Data", y = "Model") +
   theme(
      legend.position = "none",
      panel.background = element_rect(colour = "grey40", fill = "white"),
      rect = element_rect(fill = "transparent"),
      axis.title = element_text(size = 18)
   ) 

env.incid.inset = ggdraw() +
   draw_plot(env.inc) +
   draw_plot(cor.inc.env, x = 0.6, y = 0.7, width = 0.3, height = 0.3)

##
source("Random functions/multiplot.R")

pdf("Figures/Envenoming-maps.pdf", width = 7.5, height = 21)
multiplot(env.incid.inset, envs.inset, cols = 1)
dev.off()

pdf("Figures/Snakebite-maps.pdf", width = 7.5, height = 21)
multiplot(bit.incid.inset, bites.inset, cols = 1)
dev.off()

################
# Rho

rho.bites <- raster("Model-results/Spatial/Snakebites/Rho.asc")
rho.b.df <- data.frame(rasterToPoints(rho.bites))

rho.envs <- raster("Model-results/Spatial/Envenomings/Rho.asc")
rho.e.df <- data.frame(rasterToPoints(rho.envs))

cool <- colorRampPalette(c("orangered1", "goldenrod1", "darkturquoise"))

rho.b.pl <- ggplot(rho.b.df) + geom_raster(aes(x = x, y = y, fill = Rho)) +
   scale_fill_gradientn(colours = rev(cool(100))) +
   geom_polygon(data = sl,
                aes(x=long, y=lat, group = group),
                fill=NA, color="grey50", size=1.25) +
   coord_fixed(ratio = 1) + 
   labs(fill = "", title = "Random effects")+
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.95, 0.75),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(7.5, units = "mm"),
         plot.title = element_text(size = 20,
                                   hjust = 0.5),
         axis.ticks = element_blank(),
         axis.text = element_text(size = 0),
         axis.title = element_text(size = 0))

rho.e.pl <- ggplot(rho.e.df) + geom_raster(aes(x = x, y = y, fill = Rho)) +
   scale_fill_gradientn(colours = rev(cool(100))) +
   geom_polygon(data = sl,
                aes(x=long, y=lat, group = group),
                fill=NA, color="grey50", size=1.25) +
   coord_fixed(ratio = 1) + 
   labs(fill = "", title = "Random effects")+
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.95, 0.75),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(7.5, units = "mm"),
         plot.title = element_text(size = 20,
                                   hjust = 0.5),
         axis.ticks = element_blank(),
         axis.text = element_text(size = 0),
         axis.title = element_text(size = 0))

pdf("Figures/Rho-snakebites.pdf", width = 7.5/2, height = 21/4)
   rho.b.pl
dev.off()

pdf("Figures/Rho-envenomings.pdf", width = 7.5/2, height = 21/4)
   rho.e.pl
dev.off()

### RMSE

rmse.sb <- readRDS("Model-results/SLA-disctricts-CAR-RMSE.rds")
rmse.env <- readRDS("Model-results/SLA-districts-CAR-RMSE-envs.rds")

# Snakebites

rmse.sb.pl <- ggplot(rmse.sb, aes(x = long, y = lat, group = group, fill = rmse)) + 
   geom_polygon(colour = "grey50") +
   scale_fill_gradientn(colours = rev(cool(100))) +
   coord_equal() +
   labs(x = "x", y = "y", title = "RMSE", fill = "") +
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.9, 0.8),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(7.5, units = "mm"),
         plot.title = element_text(size = 20,
                                   hjust = 0.5),
         axis.ticks = element_blank(),
         axis.text = element_text(size = 0),
         axis.title = element_text(size = 0))

rmse.env.pl <- ggplot(rmse.env, aes(x = long, y = lat, group = group, fill = rmse.env)) + 
   geom_polygon(colour = "grey50") +
   scale_fill_gradientn(colours = rev(cool(100))) +
   coord_equal()+
   labs(x = "x", y = "y", title = "RMSE", fill = "") +
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.9, 0.8),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(7.5, units = "mm"),
         plot.title = element_text(size = 20,
                                   hjust = 0.5),
         axis.ticks = element_blank(),
         axis.text = element_text(size = 0),
         axis.title = element_text(size = 0))

pdf("Figures/RMSE-snakebites.pdf", width = 7.5/2, height = 21/4)
   rmse.sb.pl
dev.off()

pdf("Figures/RMSE-envenomings.pdf", width = 7.5/2, height = 21/4)
   rmse.env.pl
dev.off()

## Residuals

## Incidence estimates

res.bit.inc <- ggplot(incid.df) + geom_raster(aes(x = x, y = y, fill = (data - model))) +
   scale_fill_gradientn(colours = rev(cool(100))) +
   geom_polygon(data = sl,
                aes(x=long, y=lat, group = group),
                fill=NA, color="grey50", size=1.25) +
   coord_fixed(ratio = 1) + 
   labs(fill = "", title = "Estimates residuals")+
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.95, 0.75),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(7.5, units = "mm"),
         plot.title = element_text(size = 20,
                                   hjust = 0.5),
         axis.ticks = element_blank(),
         axis.text = element_text(size = 0),
         axis.title = element_text(size = 0))

res.env.inc <- ggplot(env.inc.df) + geom_raster(aes(x = x, y = y, fill = (data - model))) +
   scale_fill_gradientn(colours = rev(cool(100))) +
   geom_polygon(data = sl,
                aes(x=long, y=lat, group = group),
                fill=NA, color="grey50", size=1.25) +
   coord_fixed(ratio = 1) + 
   labs(fill = "", title = "Estimates residuals")+
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.95, 0.75),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(7.5, units = "mm"),
         plot.title = element_text(size = 20,
                                   hjust = 0.5),
         axis.ticks = element_blank(),
         axis.text = element_text(size = 0),
         axis.title = element_text(size = 0))

res.bit.inc
res.env.inc

pdf("Figures/Bite-incid-resids.pdf", width = 7.5/2, height = 21/4)
res.bit.inc
dev.off()

pdf("Figures/Envenom-incid-resids.pdf", width = 7.5/2, height = 21/4)
res.env.inc
dev.off()
### Survey data

surv.data <- read.csv("Data/data.for.prevmap.csv")

library(sp)

coordinates(surv.data) <- ~ cluster_long + cluster_lat
proj4string(surv.data) <- CRS("+init=epsg:4326")

surv.data.sl <- spTransform(surv.data, CRSobj = CRS("+init=epsg:5235"))

surv.data.df <- data.frame(coordinates(surv.data.sl), surv.data.sl@data)

sb.incid.r <- raster("Model-results/Spatial/Snakebites/Incidence.asc")
env.incid.r <- raster("Model-results/Spatial/Envenomings/Incidence.asc")

surv.data.df$sb.pred <- extract(sb.incid.r, coordinates(surv.data.sl))
surv.data.df$env.pred <- extract(env.incid.r, coordinates(surv.data.sl))

surv.data.df$sb.prev <- with(surv.data.df, n_bitten/n_sampled)
surv.data.df$env.prev <- with(surv.data.df, n_envenom/n_sampled)

surv.data.df$sb.resids <- with(surv.data.df, sb.prev - sb.pred)
surv.data.df$env.resids <- with(surv.data.df, env.prev - env.pred)

sb.res.pl <-ggplot() + geom_polygon(data = sl, aes(x = long, y = lat, group = group), fill = "lightgrey") +
   geom_point(data = surv.data.df, aes(x = cluster_long, y = cluster_lat, colour = sb.resids), shape = 20, size = 2)+
   scale_color_gradientn(colours = rev(cool(100))) +
   coord_equal() +
   labs(x = "x", y = "y", title = "Survey data residuals", colour = "") +
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.9, 0.8),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(7.5, units = "mm"),
         plot.title = element_text(size = 20,
                                   hjust = 0.5),
         axis.ticks = element_blank(),
         axis.text = element_text(size = 0),
         axis.title = element_text(size = 0))


env.res.pl <-ggplot() + geom_polygon(data = sl, aes(x = long, y = lat, group = group), fill = "lightgrey") +
   geom_point(data = surv.data.df, aes(x = cluster_long, y = cluster_lat, colour = env.resids), shape = 20, size = 2)+
   scale_color_gradientn(colours = rev(cool(100))) +
   coord_equal() +
   labs(x = "x", y = "y", title = "Survey data residuals", colour = "") +
   theme(panel.background =  element_rect(colour = NA, fill = NA),
         legend.position = c(0.9, 0.8),
         legend.background = element_rect("transparent"),
         legend.text = element_text(size = 12),
         legend.key.height = unit(5, units = "mm"),
         legend.key.width = unit(7.5, units = "mm"),
         plot.title = element_text(size = 20, 
                                   hjust = 0.5),
         axis.ticks = element_blank(),
         axis.text = element_text(size = 0),
         axis.title = element_text(size = 0))

pdf("Figures/Survey-residuals-snakebites.pdf", width = 7.5/2, height = 21/4)
   sb.res.pl
dev.off()

pdf("Figures/Survey-residuals-envenomings.pdf", width = 7.5/2, height = 21/4)
   env.res.pl
dev.off()



