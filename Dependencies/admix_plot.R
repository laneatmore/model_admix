##plot admixture!
setwd("../Admixture/")

pdf("admix_plot2.pdf")
admixtable = read.table("../Admixture/pruned_model.2.Q")
barplot(t(as.matrix(admixtable)), col = rainbow(2), 
        main = "K=2",
        xlab = "Individual", ylab = "Ancestry", 
        border = NA, width = 1, space = 0)

dev.off()

pdf("admix_plot3.pdf")
admixtable = read.table("../Admixture/pruned_model.3.Q")
barplot(t(as.matrix(admixtable)), col = rainbow(3), 
        main = "K=3",
        xlab = "Individual", ylab = "Ancestry", 
        border = NA, width = 1, space = 0)

dev.off()

