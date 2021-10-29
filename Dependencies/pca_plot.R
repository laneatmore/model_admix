library(ggplot2)
setwd("../Admixture/")
pdf("pca_plot.pdf")
eigenvec_table <- read.table("../Admixture/model_pca.eigenvec")
populations <- read.table("../Admixture/population_information.txt", header = TRUE)
names(eigenvec_table)[2] = "ID"
Merged <- merge(eigenvec_table, populations, by = "ID")
eigenvec_plot <- ggplot(Merged, aes(x = V3, y = V4, color = factor(POP))) + 
                    geom_point(size = 0.5) +
                    labs(title = "PCA",
                         x = "PC1", y = "PC2") +
                    theme_bw() +
                    theme(plot.title = element_text(hjust=0.5)) +
                    guides(color=guide_legend("Populations"))

eigenvec_plot

dev.off()
