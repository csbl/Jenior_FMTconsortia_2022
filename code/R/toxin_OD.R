


toxin <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/toxin_OD.tsv', sep='\t', header=TRUE)

control <- subset(toxin, group == 'control')[,4]
competitive <- subset(toxin, group == 'competitive')[,4]
cooperative <- subset(toxin, group == 'cooperative')[,4]

round(t.test(control, competitive, exact=FALSE)$p.value, 5)
round(t.test(control, cooperative, exact=FALSE)$p.value, 5)
round(t.test(competitive, cooperative, exact=FALSE)$p.value, 5)


pdf(file='~/Desktop/repos/Jenior_Consortia_2022/results/toxin_OD.pdf', width=2.5, height=4)

par(mar=c(3.5,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,70), xlim=c(0,3), 
     ylab='Rel. Toxin Equiv.', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(control, cex=0, lwd=2, at=0.5, col='gray40',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(competitive, cex=0, lwd=2, at=1.5, col='firebrick3',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(cooperative, cex=0, lwd=2, at=2.5, col='dodgerblue4',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(0,70,10), cex.axis=0.8, lwd=2)
segments(x0=0.5, y0=70*0.87, x1=1.5, lwd=2)
segments(x0=0.5, y0=70*0.92, x1=2.5, lwd=2)
segments(x0=1.5, y0=70*0.97, x1=2.5, lwd=2)
text(x=1, y=70*0.89, 'n.s.', cex=0.7)
text(x=2, y=70*0.99, 'n.s.', cex=0.7)
text(x=1.5, y=70*0.94, '*', cex=1.1, font=2)
par(xpd=TRUE)
text(x=c(0.5,1.5,2.5), y=-11, srt=40, cex=0.7, labels=c('FMT\ncontrol','Competitive\nconsortia','Cooperative\nconsortia'))
par(xpd=FALSE)

dev.off()

