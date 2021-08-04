jet.3 <- colorRampPalette(c("midnightblue", "cyan", "yellow", "darkgoldenrod1", "red", "darkred", "darkred"))
red.yl.bk  <- colorRampPalette(c("black","black", "yellow","darkgoldenrod1", "red", "red"))

dir.create("Figures")

pdf("Figures/Response-shapes-SI-inverse.pdf", width = 5, height = 7.5)
par(mar = c(2,2,4,3))
par(mfrow = c(3, 2))
for(i in 1:5){
      image2D(x = x[[i]], y = y[[i]], z = z.inc[[i]], ticktype = "detailed",
              col = rev(another.1(100)),
              xlab = "Human density", ylab = "Snakes",
              main = cover.classes[i],
              cex.main = 1, contour = T,
              cex.contour = 2,
              colkey = F
      )
}

pdf("Figures/Response-shapes-SI-jet.pdf", width = 5, height = 7.5)
par(mar = c(2,2,4,3))
par(mfrow = c(3, 2))
for(i in 1:5){
   image2D(x = x[[i]], y = y[[i]], z = z.inc[[i]], ticktype = "detailed",
           col = jet.3(100),
           xlab = "Human density", ylab = "Snakes",
           main = cover.classes[i],
           cex.main = 1, contour = T,
           cex.contour = 2,
           colkey = F
   )
}
dev.off()

#
library(ggplot2)

pdf("Figures/Snakebites-CAR-RMSE.pdf", width = 5, height = 7.5)
ggplot(sla.fort, aes(x = long, y = lat, group = group, fill = rmse)) + geom_polygon(colour = "lightgrey") +
   scale_fill_gradientn(colours = another.1(24)) +
   labs(fill = "RMSE") +
   coord_equal() +
   theme_minimal()
dev.off()


