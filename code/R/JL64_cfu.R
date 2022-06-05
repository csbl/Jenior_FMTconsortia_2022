
cfu <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/infections/jl64_cfu.tsv', sep='\t', header=TRUE)

cfu$day_1 <- log10(cfu$day_1)
cfu$day_2 <- log10(cfu$day_2)

no_treatment <- subset(cfu, treatment == 'none')
fmt <- subset(cfu, treatment == 'FMT')
positive <- subset(cfu, treatment == 'positive_(competition)')
negative <- subset(cfu, treatment == 'negative_(mutualism)')

none_pval <- wilcox.test(no_treatment$day_1, no_treatment$day_2, exact=FALSE)$p.value
fmt_pval <- wilcox.test(fmt$day_1, fmt$day_2, exact=FALSE)$p.value
positive_pval <- wilcox.test(positive$day_1, positive$day_2, exact=FALSE)$p.value
negative_pval <- wilcox.test(negative$day_1, negative$day_2, exact=FALSE)$p.value

ymax <- ceiling(max(c(max(cfu$day_1),max(cfu$day_2))))
day1_col <- 'seagreen2'
day2_col <- 'royalblue1'

png(filename='~/Desktop/cfu.png', units='in', width=6, height=4, res=300)

par(mar=c(3,3,0.5,0.5), las=1, mgp=c(1.5,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,ymax), xlim=c(0.05,4.7), yaxt='n', 
     ylab='CFU/mL (Log10)', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
axis(side=2, at=c(0:ymax), cex.axis=0.8, lwd=2)
abline(v=c(1.125,2.375,3.625), lwd=2)
text(x=c(0.5,1.75,3,4.25), y=8.5, cex=0.8,
     labels=c('No Treatment','High Competition\nConsortia','Low Competition\nConsortia','Conventional\nFecal Slurry'))
par(xpd=TRUE)
text(x=c(0.25,0.75, 1.5,2, 2.75,3.25, 4,4.5), y=-1.25, labels=c('Day 1','Day 2'), cex=1, srt=45)
par(xpd=FALSE)

boxplot(no_treatment$day_1, cex=0, lwd=2, at=0.25, col=day1_col, ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(no_treatment$day_2, cex=0, lwd=2, at=0.75, col=day2_col, ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)

boxplot(fmt$day_1, cex=0, lwd=2, at=1.5, col=day1_col, ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(fmt$day_2, cex=0, lwd=2, at=2, col=day2_col, ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)

boxplot(positive$day_1, cex=0, lwd=2, at=2.75, col=day1_col, ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(positive$day_2, cex=0, lwd=2, at=3.25, col=day2_col, ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)

boxplot(negative$day_1, cex=0, lwd=2, at=4, col=day1_col, ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(negative$day_2, cex=0, lwd=2, at=4.5, col=day2_col, ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)

dev.off()

