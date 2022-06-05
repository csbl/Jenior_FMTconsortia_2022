

# CFU
cfu <- read.delim('~/Desktop/Jenior_Consortia_2022/data/infections/JL66_cecalCFU.tsv', sep='\t', header=TRUE)
for (x in colnames(cfu)) {cfu[,x] <- as.character(cfu[,x])}; rm(x)
cfu <- subset(cfu, replicate == '1')
cfu$neat <- NULL
cfu$minus1 <- NULL
cfu$minus2 <- NULL
cfu[cfu == 'TNTC'] <- 0

# Dilution factor
cfu$minus3 <- as.numeric(cfu$minus3) * 1000
cfu$minus4 <- as.numeric(cfu$minus4) * 10000
cfu$minus5 <- as.numeric(cfu$minus5) * 100000
cfu$minus6 <- as.numeric(cfu$minus6) * 1000000
cfu$minus7 <- as.numeric(cfu$minus7) * 10000000
cfu$minus3 <- cfu$minus3 * 100 # Plating dilution
cfu$minus4 <- cfu$minus4 * 100 # Plating dilution
cfu$minus5 <- cfu$minus5 * 100 # Plating dilution

control <- subset(cfu, group == 'control')
control$mouse <- NULL
control$replicate <- NULL
control$group <- NULL
control <- as.vector(t(control))
control <- control[control>0]
control <- sample(control, 7, replace=FALSE)
competitive <- subset(cfu, group == 'competitive')
competitive$mouse <- NULL
competitive$replicate <- NULL
competitive$group <- NULL
competitive <- as.vector(t(competitive))
competitive <- competitive[competitive>0]
competitive <- sample(competitive, 7, replace=FALSE)
cooperative <- subset(cfu, group == 'cooperative')
cooperative$mouse <- NULL
cooperative$replicate <- NULL
cooperative$group <- NULL
cooperative <- as.vector(t(cooperative))
cooperative <- cooperative[cooperative>0]
cooperative <- sample(cooperative, 7, replace=FALSE)
rm(cfu)
control_cfu <- c(5.0e+07,3.0e+06,1.5e+07,3.4e+07,4.0e+07,1.0e+07,1.3e+07)
competitive_cfu <- c(3.5e+06,1.2e+08,3.0e+07,1.0e+07,3.0e+07,2.9e+07,1.9e+07)
cooperative_cfu <- c(2.0e+07,5.0e+07,3.0e+07,2.2e+07,1.2e+07,2.3e+07,8.0e+07)
round(wilcox.test(control_cfu, competitive_cfu, exact=FALSE)$p.value, 4)
round(wilcox.test(control_cfu, cooperative_cfu, exact=FALSE)$p.value, 4)
round(wilcox.test(cooperative_cfu, competitive_cfu, exact=FALSE)$p.value, 4)

# Toxin
toxin <- read.delim('~/Desktop/Jenior_Consortia_2022/data/toxin_OD.tsv', sep='\t', header=TRUE)
control_tox <- subset(toxin, group == 'control')[,4]
competitive_tox <- subset(toxin, group == 'competitive')[,4]
cooperative_tox <- subset(toxin, group == 'cooperative')[,4]
round(t.test(control_tox, competitive_tox, exact=FALSE)$p.value, 5)
round(t.test(control_tox, cooperative_tox, exact=FALSE)$p.value, 5)
round(t.test(competitive_tox, cooperative_tox, exact=FALSE)$p.value, 5)

# Generate figure
library(plotrix)
pdf(file='~/Desktop/Jenior_Consortia_2022/results/cfu_toxin.pdf', width=2.5, height=5)
layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE))

par(mar=c(3.5,2.5,0.5,0.5), las=1, mgp=c(1.2,0.7,0), xpd=FALSE, lwd=1.8)
plot(0, type='n', ylim=c(4,10), xlim=c(0,3), 
     ylab='CFU/g (Log10)', xlab='', xaxt='n', yaxt='n', cex.lab=0.9)
boxplot(log10(control_cfu), cex=0, lwd=1.8, at=0.5, col='gray40', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
boxplot(log10(competitive_cfu), cex=0, lwd=1.8, at=1.5, col='dodgerblue4',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
boxplot(log10(cooperative_cfu), cex=0, lwd=1.8, at=2.5, col='firebrick3',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(4,10,2), labels=c(0,seq(6,10,2)), cex.axis=0.7, lwd=1.7)
box(lwd=1.8)
axis.break(2, 5, style='zigzag', brw=0.04)
segments(x0=0.5, y0=8.5, x1=1.5, lwd=1.8)
segments(x0=0.5, y0=9, x1=2.5, lwd=1.8)
segments(x0=1.5, y0=9.5, x1=2.5, lwd=1.8)
text(x=c(1,1.5,2), y=c(8.7,9.2,9.7), 'n.s.', cex=0.6)
par(xpd=TRUE)
text(x=c(0.5,1.5,2.5), y=2.5, srt=40, cex=0.75, labels=c('Control\ngavage','Competitive\nconsortia','Cooperative\nconsortia'))
par(xpd=FALSE)

par(mar=c(3.5,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=1.8)
plot(0, type='n', ylim=c(0,75), xlim=c(0,3), 
     ylab='Rel. Toxin Equiv.', xlab='', xaxt='n', yaxt='n', cex.lab=0.9)
boxplot(control_tox, cex=0, lwd=1.8, at=0.5, col='gray40', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
boxplot(competitive_tox, cex=0, lwd=1.8, at=1.5, col='dodgerblue4',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
boxplot(cooperative_tox, cex=0, lwd=1.8, at=2.5, col='firebrick3',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(0,70,10), cex.axis=0.7, lwd=1.8)
segments(x0=0.5, y0=70*0.87, x1=1.5, lwd=1.8)
segments(x0=0.5, y0=70*0.95, x1=2.5, lwd=1.8)
segments(x0=1.5, y0=70*1.03, x1=2.5, lwd=1.8)
text(x=1, y=70*0.9, 'n.s.', cex=0.6)
text(x=1.5, y=70*0.99, '*', font=2, cex=1.2)
text(x=2, y=70*1.07, '*', font=2, cex=1.2)

par(xpd=TRUE)
text(x=c(0.5,1.5,2.5), y=-19, srt=40, cex=0.75, labels=c('Control\ngavage','Competitive\nconsortia','Cooperative\nconsortia'))
par(xpd=FALSE)

dev.off()




