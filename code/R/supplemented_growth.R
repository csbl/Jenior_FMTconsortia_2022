
# Read in and format data
growth <- read.delim('~/Desktop/Jenior_Consortia_2022/data/invitro_growth/supplemented_media.tsv', sep='\t', header=TRUE)
growth1 <- read.delim('~/Desktop/Jenior_Consortia_2022/data/invitro_growth/spent_media.1.R20291.tsv', sep='\t', header=TRUE)

collateNorm <- function(col_num, data) {
  blank <- paste0('A', as.character(col_num))
  wells <- c(paste0('B', as.character(col_num)), paste0('C', as.character(col_num)), 
             paste0('D', as.character(col_num)), paste0('E', as.character(col_num)))
  blank <- as.matrix(data[,blank])
  growth <- as.matrix(data[,wells])
  blank <- apply(blank, 1, quantile, probs=c(0.25,0.5,0.75))
  growth <- apply(growth, 1, quantile, probs=c(0.25,0.5,0.75))
  growth <- abs(growth - blank)
  growth <- t(growth)
  growth <- growth[1:132,] # 12 hours
  return(growth)}

library(bayestestR)
analyze_growth <- function(data) {
  data <- data[1:132]
  max_rate <- round(diff(data)[which.max(diff(data))], digits=3) # Maximum growth rate
  time_max_rate <- round((which.max(diff(data))*6.390534)/60, digits=3) # Time of maximum growth rate
  max_od <- round(max(data), digits=3) # Maximum OD
  time_max_od <- round((which.max(data)*6.390534)/60, digits=3) # Time of max OD
  area_under <- round(area_under_curve(c(1:132), data), digits=3) # Area under curve
  growth_data <- c(max_rate, time_max_rate, max_od, time_max_od, area_under)
  return(growth_data)}

#---------------------------------------------------------------------------------------------#

# Calculate curve intervals
bhi <- collateNorm(1, growth) # BHI
dep <- collateNorm(2, growth1) # R20291 depleted
dep_bhi <- collateNorm(2, growth) # R20291 depleted + BHI
dep_comp <- collateNorm(3, growth) # R20291 depleted + Competative depleted
dep_coop <- collateNorm(4, growth) # R20291 depleted + Cooperative depleted

# Calculate curve summaries
dep_bhi_summary1 <- analyze_growth(growth$B2)
dep_bhi_summary2 <- analyze_growth(growth$C2)
dep_bhi_summary3 <- analyze_growth(growth$D2)
dep_bhi_summary4 <- analyze_growth(growth$E2)
dep_bhi_summary <- as.data.frame(rbind(dep_bhi_summary1, dep_bhi_summary2, dep_bhi_summary3, dep_bhi_summary4))
rownames(dep_bhi_summary) <- c('dep_bhi_1','dep_bhi_2','dep_bhi_3','dep_bhi_4')
colnames(dep_bhi_summary) <- c('max_growth_rate', 'hours_to_max_rate', 'max_density', 'hours_to_max_density', 'auc')
rm(dep_bhi_summary1, dep_bhi_summary2, dep_bhi_summary3, dep_bhi_summary4)
dep_summary1 <- analyze_growth(growth1$B2)
dep_summary2 <- analyze_growth(growth1$C2)
dep_summary3 <- analyze_growth(growth1$D2)
dep_summary4 <- analyze_growth(growth1$E2)
dep_summary <- as.data.frame(rbind(dep_summary1, dep_summary2, dep_summary3, dep_summary4))
rownames(dep_summary) <- c('dep_1','dep_2','dep_3','dep_4')
colnames(dep_summary) <- c('max_growth1_rate', 'hours_to_max_rate', 'max_density', 'hours_to_max_density', 'auc')
rm(dep_summary1, dep_summary2, dep_summary3, dep_summary4)
dep_comp_summary1 <- analyze_growth(growth$B3)
dep_comp_summary2 <- analyze_growth(growth$C3)
dep_comp_summary3 <- analyze_growth(growth$D3)
dep_comp_summary4 <- analyze_growth(growth$E3)
dep_comp_summary <- as.data.frame(rbind(dep_comp_summary1, dep_comp_summary2, dep_comp_summary3, dep_comp_summary4))
rownames(dep_comp_summary) <- c('dep_comp_1','dep_comp_2','dep_comp_3','dep_comp_4')
colnames(dep_comp_summary) <- c('max_growth_rate', 'hours_to_max_rate', 'max_density', 'hours_to_max_density', 'auc')
rm(dep_comp_summary1, dep_comp_summary2, dep_comp_summary3, dep_comp_summary4)
dep_coop_summary1 <- analyze_growth(growth$B4)
dep_coop_summary2 <- analyze_growth(growth$C4)
dep_coop_summary3 <- analyze_growth(growth$D4)
dep_coop_summary4 <- analyze_growth(growth$E4)
dep_coop_summary <- as.data.frame(rbind(dep_coop_summary1, dep_coop_summary2, dep_coop_summary3, dep_coop_summary4))
rownames(dep_coop_summary) <- c('dep_coop_1','dep_coop_2','dep_coop_3','dep_coop_4')
colnames(dep_coop_summary) <- c('max_growth_rate', 'hours_to_max_rate', 'max_density', 'hours_to_max_density', 'auc')
rm(dep_coop_summary1, dep_coop_summary2, dep_coop_summary3, dep_coop_summary4)

# Test differences
pval1 <- round(wilcox.test(dep_summary$auc, dep_comp_summary$auc, exact=FALSE)$p.value, 4)
pval2 <- round(wilcox.test(dep_summary$auc, dep_coop_summary$auc, exact=FALSE)$p.value, 4)
pval3 <- round(wilcox.test(dep_coop_summary$auc, dep_comp_summary$auc, exact=FALSE)$p.value, 4)

#---------------------------------------------------------------------------------------------#

library(scales)
draw_growth <- function(data, clr) {
  polygon(c(1:nrow(data), rev(1:nrow(data))), c(data[,3], rev(data[,1])), col=alpha(clr, 0.25), border=NA)
  lines(data[,2], lwd=3, col=clr)}


par(mar=c(3,3,0.5,0.5), mgp=c(2, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(dep_bhi[,2], type='l', ylab='OD600', xaxt='n', yaxt='n', xlab='Hours', ylim=c(0,0.8), xlim=c(1,145))
axis(side=1, at=seq(1,176,16), labels=seq(0,20,2), cex.axis=0.8, lwd=2)
axis(side=2, at=seq(0,0.8,0.2), cex.axis=0.8, lwd=2)
abline(h=c(0.2,0.4,0.6), col='gray', lty=5, lwd=2)
draw_growth(dep, 'purple')
draw_growth(dep_bhi, 'orange')
draw_growth(dep_coop, 'red')
draw_growth(dep_comp, 'blue')



pdf(file='~/Desktop/Jenior_Consortia_2022/results/supplemented.pdf', width=5, height=3)
par(mar=c(6,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,80), xlim=c(0,4), 
     ylab='AUC', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(dep_summary$auc, cex=0, lwd=2, at=0.5, col='gray40', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(dep_bhi_summary$auc, cex=0, lwd=2, at=1.5, col='orange', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(dep_coop_summary$auc, cex=0, lwd=2, at=2.5, col='firebrick3', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
boxplot(dep_comp_summary$auc, cex=0, lwd=2, at=3.5, col='dodgerblue4', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(0,80,20), cex.axis=0.8, lwd=2)
segments(x0=2.5, x1=3.5, y0=50)
text(x=3, y=54, '*', cex=1.5, font=2)
par(xpd=TRUE)
text(x=c(0.5,1.5,2.5,3.5), y=-30, srt=45, cex=0.8,
     c('R20291 depleted','R20291 depleted +\nFresh media','R20291 depleted +\nCooperative depleted','R20291 depleted +\nCompetative depleted'))
dev.off()


