
# Read in and subset data
colonies <- read.delim('~/Desktop/Jenior_Consortia_2022/data/JL66_microscopy/colony_metrics.tsv', sep='\t', header=TRUE)
control <- subset(colonies, treatment == 'control')
competitive <- subset(colonies, treatment == 'competitive')
cooperative <- subset(colonies, treatment == 'cooperative')

#---------------------------------------------------------------------------#

# Calculate differences 

# Area
area_pvals <- c(round(wilcox.test(control$area, competitive$area, exact=FALSE)$p.value, 4),
                round(wilcox.test(control$area, cooperative$area, exact=FALSE)$p.value, 4),
                round(wilcox.test(competitive$area, cooperative$area, exact=FALSE)$p.value, 4))
area_pvals <- p.adjust(area_pvals, method='BH')

# Perimeter
perimeter_pvals <- c(round(wilcox.test(control$perimeter, competitive$perimeter, exact=FALSE)$p.value, 4),
                     round(wilcox.test(control$perimeter, cooperative$perimeter, exact=FALSE)$p.value, 4),
                     round(wilcox.test(competitive$perimeter, cooperative$perimeter, exact=FALSE)$p.value, 4))
perimeter_pvals <- p.adjust(perimeter_pvals, method='BH')

#---------------------------------------------------------------------------#

# Generate figure
pdf(file='~/Desktop/Jenior_Consortia_2022/results/colonies.pdf', width=2.5, height=4)

ymax <- signif(as.numeric(max(colonies$area))*1.2, digits=2)
par(mar=c(3.5,4,0.5,0.5), las=1, mgp=c(3,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(50000,800000), xlim=c(0,3), 
     ylab='Colony Area (AU)', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(control$area, cex=0, lwd=2, at=0.5, col='gray40',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(competitive$area, cex=0, lwd=2, at=1.5, col='dodgerblue4',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(cooperative$area, cex=0, lwd=2, at=2.5, col='firebrick3',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(0,800000,100000), cex.axis=0.8, lwd=2)
segments(x0=0.5, y0=ymax*0.85, x1=1.5, lwd=2)
segments(x0=0.5, y0=ymax*0.9, x1=2.5, lwd=2)
segments(x0=1.5, y0=ymax*0.95, x1=2.5, lwd=2)
text(x=1, y=ymax*0.87, 'n.s.', cex=0.7)
text(x=c(1.5,2), y=c(ymax*0.92,ymax*0.97), '**', cex=1.1, font=2)
par(xpd=TRUE)
text(x=c(0.5,1.5,2.5), y=-60500, srt=40, cex=0.75, labels=c('Control\ngavage','Competitive\nconsortia','Cooperative\nconsortia'))
par(xpd=FALSE)

dev.off()

ymax <- signif(as.numeric(max(colonies$perimeter))*1.2, digits=2)
par(mar=c(3.5,3.8,0.5,0.5), las=1, mgp=c(2.6,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,5000), xlim=c(0,3), 
     ylab='Colony Perimeter (AU)', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(control$perimeter, cex=0, lwd=2, at=0.5, col='gray40',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(competitive$perimeter, cex=0, lwd=2, at=1.5, col='dodgerblue4',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(cooperative$perimeter, cex=0, lwd=2, at=2.5, col='firebrick3',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(0,5000,1000), cex.axis=0.8, lwd=2)
segments(x0=0.5, y0=ymax*0.85, x1=1.5, lwd=2)
segments(x0=0.5, y0=ymax*0.9, x1=2.5, lwd=2)
segments(x0=1.5, y0=ymax*0.95, x1=2.5, lwd=2)
text(x=1, y=ymax*0.87, '*', cex=1.1, font=2)
text(x=c(1.5,2), y=c(ymax*0.92,ymax*0.97), '**', cex=1.1, font=2)
par(xpd=TRUE)
text(x=c(0.5,1.5,2.5), y=-750, srt=40, cex=0.7, labels=c('Control\ngavage','Competitive\nconsortia','Cooperative\nconsortia'))
par(xpd=FALSE)





