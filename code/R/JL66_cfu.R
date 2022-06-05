
# Read in and format
cfu <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/infections/JL66_cecalCFU.tsv', sep='\t', header=TRUE)
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

# Plating dilution
cfu$minus3 <- cfu$minus3 * 100
cfu$minus4 <- cfu$minus4 * 100
cfu$minus5 <- cfu$minus5 * 100

# Subset treatments
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
control <- c(5.0e+07,3.0e+06,1.5e+07,3.4e+07,4.0e+07,1.0e+07,1.3e+07)
competitive <- c(3.5e+06,1.2e+08,3.0e+07,1.0e+07,3.0e+07,2.9e+07,1.9e+07)
cooperative <- c(2.0e+07,5.0e+07,3.0e+07,2.2e+07,1.2e+07,2.3e+07,8.0e+07)


# Test differences
pval_cont_comp <- round(wilcox.test(control, competitive, exact=FALSE)$p.value, 4)
pval_cont_coop <- round(wilcox.test(control, cooperative, exact=FALSE)$p.value, 4)
pval_comp_coop <- round(wilcox.test(cooperative, competitive, exact=FALSE)$p.value, 4)

# Generate figure
library(plotrix)
png(filename='~/Desktop/repos/Jenior_Consortia_2022/results/cfu.png', units='in', width=2.5, height=4, res=300)
par(mar=c(3.5,2.5,0.5,0.5), las=1, mgp=c(1.2,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(4,10), xlim=c(0,3), 
     ylab='CFU/g (Log10)', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(log10(control), cex=0, lwd=2, at=0.5, col='gray40',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(log10(competitive), cex=0, lwd=2, at=1.5, col='firebrick3',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(log10(cooperative), cex=0, lwd=2, at=2.5, col='dodgerblue4',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(4,10,2), labels=c(0,seq(6,10,2)), cex.axis=0.8, lwd=1.7)
box(lwd=2)
axis.break(2, 5, style='slash', brw=0.04)
segments(x0=0.5, y0=8.5, x1=1.5, lwd=2)
segments(x0=0.5, y0=9, x1=2.5, lwd=2)
segments(x0=1.5, y0=9.5, x1=2.5, lwd=2)
text(x=c(1,1.5,2), y=c(8.7,9.2,9.7), 'n.s.', cex=0.9)
par(xpd=TRUE)
text(x=c(0.5,1.5,2.5), y=3, srt=40, cex=0.8, labels=c('FMT\ncontrol','Competitive\nconsortia','Cooperative\nconsortia'))
par(xpd=FALSE)
dev.off()




