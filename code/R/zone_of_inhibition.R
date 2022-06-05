

vanc <- c(11,12,14,10,12)
bhi <- c(0,0,0,0,0)
coop <- c(0,0,0,0,0)
comp <- c(0,0,0,0,0)


pdf(file='~/Desktop/Jenior_Consortia_2022/results/inhibition.pdf', width=4, height=3)
par(mar=c(6,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,15), xlim=c(0,4), 
     ylab='Zone of Inhibition (mm)', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
axis(side=2, at=seq(0,15,3), cex.axis=0.8, lwd=2)
boxplot(vanc, cex=0, lwd=2, at=0.5, col='gray40', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(bhi, cex=0, lwd=2, at=1.5, col='orange', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(coop, cex=0, lwd=2, at=2.5, col='firebrick3', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(comp, cex=0, lwd=2, at=3.5, col='dodgerblue4', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
par(xpd=TRUE)
text(x=c(0.5,1.5,2.5,3.5), y=-5, srt=45, cex=0.8,
     c('Vancomycin\n(ug/ml)','R20291\ndepleted media','Cooperative\ndepleted media','Competative\ndepleted media'))
dev.off()

