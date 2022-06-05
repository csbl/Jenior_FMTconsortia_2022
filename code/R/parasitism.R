
pdf(file='~/Desktop/Jenior_Consortia_2022/results/parasitism.pdf', width=4, height=4)
par(mar=c(3,3,1,1), las=1, mgp=c(2,0.9,0), lwd=2)
barplot(c(89,578), names.arg=c('iCdR703 to\nother GENRE','Other GENRE\nto iCdR703'), yaxt='n',
        col='gray40', width=0.5, space=0.25,  ylim=c(0,600), ylab='Metabolites shared')
box()
axis(2, at=seq(0,600,100), lwd=2, cex.axis=0.8)
dev.off()
