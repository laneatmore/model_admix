library(ggplot2)
setwd("../Admixture/")

pdf("cv_error_plot.pdf")
error = read.table("../Admixture/cv_error.txt")
cv_plot = ggplot(error, aes(x = V3, y = V4)) +
  geom_point(size = 2) + xlab("K-Value") +
  ylab("Error")

cv_plot

dev.off()

