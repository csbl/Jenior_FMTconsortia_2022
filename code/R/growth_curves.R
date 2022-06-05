
# Read in data
growth1 <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/invitro_growth/spent_media.1.R20291.tsv', sep='\t', header=TRUE)
growth2 <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/invitro_growth/spent_media.2.R20291.tsv', sep='\t', header=TRUE)
growth3 <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/invitro_growth/spent_media.3.R20291.tsv', sep='\t', header=TRUE)
consortia1 <- read.delim('~/Desktop/repos/Jenior_FMT_2021/data/invitro_growth/consortia.1.R20291.tsv', sep='\t', header=TRUE)

#---------------------------------------------------------------------------------------------#

# Define functions
collateNorm <- function(col_num, data) {
  blank <- paste0('A', as.character(col_num))
  wells <- c(paste0('B', as.character(col_num)), paste0('C', as.character(col_num)), 
             paste0('D', as.character(col_num)), paste0('E', as.character(col_num)), 
             paste0('F', as.character(col_num)), paste0('G', as.character(col_num)))
  blank <- as.matrix(data[,blank])
  growth <- as.matrix(data[,wells])
  blank <- apply(blank, 1, quantile, probs=c(0.25,0.5,0.75))
  growth <- apply(growth, 1, quantile, probs=c(0.25,0.5,0.75))
  growth <- abs(growth - blank)
  growth <- t(growth)
  
  return(growth)}

collateNorm2 <- function(col_num, data) {
  blank <- paste0('A', as.character(col_num))
  wells <- c(paste0('B', as.character(col_num)), paste0('C', as.character(col_num)), 
             paste0('D', as.character(col_num)), paste0('E', as.character(col_num)))
  blank <- as.matrix(data[,blank])
  growth <- as.matrix(data[,wells])
  blank <- apply(blank, 1, quantile, probs=c(0.25,0.5,0.75))
  growth <- apply(growth, 1, quantile, probs=c(0.25,0.5,0.75))
  growth <- abs(growth - blank)
  growth <- t(growth)
  
  return(growth)}

library(bayestestR)
analyze_growth <- function(data) {
  max_rate <- round(diff(data)[which.max(diff(data))], digits=3) # Maximum growth rate
  time_max_rate <- round((which.max(diff(data))*6.390534)/60, digits=3) # Time of maximum growth rate
  max_od <- round(max(data), digits=3) # Maximum OD
  time_max_od <- round((which.max(data)*6.390534)/60, digits=3) # Time of max OD
  area_under <- round(area_under_curve(c(1:169), data), digits=3) # Area under curve
  growth_data <- c(max_rate, time_max_rate, max_od, time_max_od, area_under)
  return(growth_data)}

#---------------------------------------------------------------------------------------------#

# Collate and normalize results
blank <- collateNorm2(1, consortia1)[1:169,]
r20291 <- collateNorm(2, growth1) # (self) Clostridioides difficile R20291
dsm2950 <- collateNorm(3, growth2) # Blautia producta
atcc55813 <- collateNorm(7, growth3) # Bifidobacterium longum 
blank_all <- as.matrix(consortia1[,c('B1','C1','D1','E1')])[1:169,]
r20291_all <- as.matrix(growth1[,c('B2','C2','D2','E2','F2','G2')])
dsm2950_all <- as.matrix(growth2[,c('B3','C3','D3','E3','F3','G3')])
atcc55813_all <- as.matrix(growth3[,c('B7','C7','D7','E7','F7','G7')])
rm(collateNorm, collateNorm2, growth1, growth2, growth3, consortia1)

# Calculate curve metrics
blank_summary1 <- analyze_growth(blank_all[,1])
blank_summary2 <- analyze_growth(blank_all[,2])
blank_summary3 <- analyze_growth(blank_all[,3])
blank_summary4 <- analyze_growth(blank_all[,4])
blank_summary <- as.data.frame(rbind(blank_summary1, blank_summary2, blank_summary3, blank_summary4))
rownames(blank_summary) <- c('blank_1','blank_2','blank_3','blank_4')
colnames(blank_summary) <- c('max_growth_rate', 'hours_to_max_rate', 'max_density', 'hours_to_max_density', 'auc')
rm(blank_summary1, blank_summary2, blank_summary3, blank_summary4)
r20291_summary1 <- analyze_growth(r20291_all[,1])
r20291_summary2 <- analyze_growth(r20291_all[,2])
r20291_summary3 <- analyze_growth(r20291_all[,3])
r20291_summary4 <- analyze_growth(r20291_all[,4])
r20291_summary5 <- analyze_growth(r20291_all[,5])
r20291_summary6 <- analyze_growth(r20291_all[,6])
r20291_summary <- as.data.frame(rbind(r20291_summary1, r20291_summary2, r20291_summary3, r20291_summary4, r20291_summary5, r20291_summary6))
rownames(r20291_summary) <- c('r20291_1','r20291_2','r20291_3','r20291_4','r20291_5','r20291_6')
colnames(r20291_summary) <- c('max_growth_rate', 'hours_to_max_rate', 'max_density', 'hours_to_max_density', 'auc')
rm(r20291_summary1, r20291_summary2, r20291_summary3, r20291_summary4, r20291_summary5, r20291_summary6)
dsm2950_summary1 <- analyze_growth(dsm2950_all[,1])
dsm2950_summary2 <- analyze_growth(dsm2950_all[,2])
dsm2950_summary3 <- analyze_growth(dsm2950_all[,3])
dsm2950_summary4 <- analyze_growth(dsm2950_all[,4])
dsm2950_summary5 <- analyze_growth(dsm2950_all[,5])
dsm2950_summary6 <- analyze_growth(dsm2950_all[,6])
dsm2950_summary <- as.data.frame(rbind(dsm2950_summary1, dsm2950_summary2, dsm2950_summary3, dsm2950_summary4, dsm2950_summary5, dsm2950_summary6))
rownames(dsm2950_summary) <- c('dsm2950_1','dsm2950_2','dsm2950_3','dsm2950_4','dsm2950_5','dsm2950_6')
colnames(dsm2950_summary) <- c('max_growth_rate', 'hours_to_max_rate', 'max_density', 'hours_to_max_density', 'auc')
rm(dsm2950_summary1, dsm2950_summary2, dsm2950_summary3, dsm2950_summary4, dsm2950_summary5, dsm2950_summary6)
atcc55813_summary1 <- analyze_growth(atcc55813_all[,1])
atcc55813_summary2 <- analyze_growth(atcc55813_all[,2])
atcc55813_summary3 <- analyze_growth(atcc55813_all[,3])
atcc55813_summary4 <- analyze_growth(atcc55813_all[,4])
atcc55813_summary5 <- analyze_growth(atcc55813_all[,5])
atcc55813_summary6 <- analyze_growth(atcc55813_all[,6])
atcc55813_summary <- as.data.frame(rbind(atcc55813_summary1, atcc55813_summary2, atcc55813_summary3, atcc55813_summary4, atcc55813_summary5, atcc55813_summary6))
rownames(atcc55813_summary) <- c('atcc55813_1','atcc55813_2','atcc55813_3','atcc55813_4','atcc55813_5','atcc55813_6')
colnames(atcc55813_summary) <- c('max_growth_rate', 'hours_to_max_rate', 'max_density', 'hours_to_max_density', 'auc')
rm(atcc55813_summary1, atcc55813_summary2, atcc55813_summary3, atcc55813_summary4, atcc55813_summary5, atcc55813_summary6)
rm(blank_all, r20291_all, dsm2950_all, atcc55813_all)

# Test differences
pval1 <- round(wilcox.test(x=r20291_summary$auc, y=dsm2950_summary$auc, eact=FALE)$p.value, 4)
pval2 <- round(wilcox.test(x=r20291_summary$auc, y=atcc55813_summary$auc, eact=FALE)$p.value, 4)
pval3 <- round(wilcox.test(x=atcc55813_summary$auc, y=dsm2950_summary$auc, eact=FALE)$p.value, 4)
pvals <- p.adjust(c(pval1,pval2,pval3), method='BH')
rm(pval1,pval2,pval3)

#---------------------------------------------------------------------------------------------#

# Generate figure
library(scales)
blank_col <- 'black'
r20291_col <- 'darkorchid3'
dsm2950_col <- 'chocolate2'
atcc55813_col <- 'cyan3'
r20291_lab <- as.expression(bquote(paste(italic('C. difficile'),' str. R20291',sep='')))
dsm2950_lab <- as.expression(bquote(paste(italic('B. producta'),' DSM 2950',sep='')))
atcc55813_lab <- as.expression(bquote(paste(italic('B. longum'),' ATCC 55813',sep='')))

pdf(file='~/Desktop/Jenior_Consortia_2022/results/growth_curves.pdf', width=5, height=3)
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))

par(mar=c(3,3,1.5,0.5), mgp=c(2, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(blank[,2], type='l', ylab='OD600', xaxt='n', yaxt='n', xlab='Hours', ylim=c(0,0.8), xlim=c(1,145), cex.main=1,
     main=substitute(paste(bolditalic('C. difficile'), bold(' R20291 growth in conditioned media'))))
axis(side=1, at=seq(1,176,16), labels=seq(0,20,2), cex.axis=0.8, lwd=2)
axis(side=2, at=seq(0,0.8,0.2), cex.axis=0.8, lwd=2)
abline(h=c(0.2,0.4,0.6), col='gray', lty=5, lwd=2)
polygon(c(1:169, rev(1:169)), c(blank[,3], rev(blank[,1])), col=alpha(blank_col, 0.25), border=NA)
lines(blank[,2], lwd=3, col=blank_col)
polygon(c(1:169, rev(1:169)), c(r20291[,3], rev(r20291[,1])), col=alpha(r20291_col, 0.25), border=NA)
lines(r20291[,2], lwd=3, col=r20291_col)
polygon(c(1:169, rev(1:169)), c(dsm2950[,3], rev(dsm2950[,1])), col=alpha(dsm2950_col, 0.25), border=NA)
lines(dsm2950[,2], lwd=3, col=dsm2950_col)
polygon(c(1:169, rev(1:169)), c(atcc55813[,3], rev(atcc55813[,1])), col=alpha(atcc55813_col, 0.25), border=NA)
lines(atcc55813[,2], lwd=3, col=atcc55813_col)
legend('topleft', legend=c('Fresh media', r20291_lab, dsm2950_lab, atcc55813_lab), 
       col=c(blank_col, r20291_col, dsm2950_col, atcc55813_col), 
       lwd=3, cex=0.9, bg='white')
box()


par(mar=c(10,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,120), xlim=c(0,4), 
     ylab='AUC', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(blank_summary$auc, cex=0, lwd=2, at=0.5, col='gray40', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(r20291_summary$auc, cex=0, lwd=2, at=1.5, col=r20291_col, ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(dsm2950_summary$auc, cex=0, lwd=2, at=2.5, col=dsm2950_col, ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(atcc55813_summary$auc, cex=0, lwd=2, at=3.5, col=atcc55813_col, ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(0,120,40), cex.axis=0.8, lwd=2)
par(xpd=TRUE)
text(x=c(0.5,1.5,2.5,3.5), y=-5, c('Fresh media',r20291_lab,dsm2950_lab,atcc55813_lab), srt=90, adj=1)
segments(x0=1.5, x1=2.5, y0=50)
text(x=2, y=54, 'n.s.', cex=0.7)
segments(x0=1.5, x1=3.5, y0=80)
text(x=2.5, y=84, '**', cex=1.1, font=2)
segments(x0=2.5, x1=3.5, y0=95)
text(x=3, y=99, '*', cex=1.1, font=2)

dev.off()

