
if (!require(scales)) install.packages('scales')
library(scales)

# C diff growth in conditioned media

# Read in and format data
growth1 <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/spent_media.1.R20291.tsv', sep='\t', header=TRUE)
growth2 <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/spent_media.2.R20291.tsv', sep='\t', header=TRUE)
growth3 <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/spent_media.3.R20291.tsv', sep='\t', header=TRUE)
growth4 <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/spent_media.4.R20291.tsv', sep='\t', header=TRUE)
growth5 <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/K12_spent_media.R20291.tsv', sep='\t', header=TRUE)
growth6 <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/spent_media.5.R20291.tsv', sep='\t', header=TRUE)
consortia1 <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/consortia.1.R20291.tsv', sep='\t', header=TRUE)
consortia2 <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/invitro_growth/consortia.2.R20291.tsv', sep='\t', header=TRUE)


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
  growth <- growth[1:169,] # 18 hours
  
  return(growth)
}

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
  growth <- growth[1:169,]
  return(growth)}

# Collate and normalize results
fresh <- collateNorm(1, growth1) # no depletion
fresh2 <- collateNorm(1, growth2) # no depletion
fresh3 <- collateNorm(1, growth3) # no depletion
fresh4 <- collateNorm(1, growth4) # no depletion
r20291 <- collateNorm(2, growth1) # (self) Clostridioides difficile R20291

dsm933 <- collateNorm(4, growth1) # Clostridium clostridioforme
dsm795 <- collateNorm(6, growth1) # Clostridium sporogenes
dsm109923 <- collateNorm(8, growth1) # Enterococcus faecium
dsm755 <- collateNorm(10, growth1) # Clostridium indolis
dsm2079 <- collateNorm(12, growth1) # Bacteroides thetaiotaomicron
dsm2950 <- collateNorm(3, growth2) # Blautia producta
dsm3979 <- collateNorm(5, growth2) # Collinsella aerofaciens
dsm1 <- collateNorm(7, growth2) # Bacillus coagulens
atcc8503 <- collateNorm(3, growth3) # Parabacteroides distasonis
atcc27919 <- collateNorm(5, growth3) # Bifidobacterium pseudocatenulatum
atcc55813 <- collateNorm(7, growth3) # Bifidobacterium longum 
lmd9 <- collateNorm(9, growth3) # Streptococcus thermophillus
atcc8482 <- collateNorm(3, growth4) # Bacteroides vulgatus
k12 <- collateNorm(11, growth5) # Escherichia coli
atcc33656 <- collateNorm(1, growth6) # Eubacterium rectale
dsm22959 <- collateNorm2(5, consortia2) # Akkermansia muciniphila
dsm14610 <- collateNorm2(3, consortia2) # Roseburia intestinalis

rm(collateNorm, growth1, growth2, growth3, growth4, growth5, growth6)






blank <- collateNorm2(1, consortia1)
consortia1 <- collateNorm2(2, consortia1)

pos_consortia <- collateNorm2(7, consortia2)
neg_consortia <- collateNorm2(9, consortia2)




library(viridis)
pall <- viridis(6)
draw_growth <- function(data, clr) {
  polygon(c(1:nrow(data), rev(1:nrow(data))), c(data[,3], rev(data[,1])), col=alpha(clr, 0.25), border=NA)
  lines(data[,2], lwd=3, col=clr)}



pos_pal <- c('blue','green','orange','deeppink3','aquamarine3','cyan2')
png(filename='~/Desktop/repos/Jenior_Consortia_2022/results/positive_strains.png', units='in', width=6, height=5, res=300)
par(mar=c(3,3,0.5,0.5), mgp=c(2, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(blank[,2], type='l', ylab='OD600', xaxt='n', yaxt='n', xlab='Hours', ylim=c(0,0.8), xlim=c(1,145))
axis(side=1, at=seq(1,176,16), labels=seq(0,20,2), cex.axis=0.8, lwd=2)
axis(side=2, at=seq(0,0.8,0.2), cex.axis=0.8, lwd=2)
abline(h=c(0.2,0.4,0.6), col='gray', lty=5, lwd=2)

draw_growth(blank, 'black')
draw_growth(atcc33656, pall[1])
draw_growth(atcc27919, pall[2])
#draw_growth(atcc8503, pall[3])
draw_growth(atcc8482, pall[3])
draw_growth(dsm2950, pall[4])
#draw_growth(dsm795, pall[6])

legend('topleft', legend=c('Fresh media',
                           'Eubacterium rectale',
                           'Bifidobacterium pseudocatenulatum',
                           'Phocaeicola vulgatus',
                           'Blautia producta'), 
       col=c('black', pall), 
       lwd=3, cex=0.9, bg='white', text.font=c(1,3,3,3,3))
box()
dev.off()



png(filename='~/Desktop/repos/Jenior_Consortia_2022/results/positive_negative_consortia.png', units='in', width=4, height=3, res=300)
par(mar=c(3,3,0.5,0.5), mgp=c(2, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(blank[,2], type='l', ylab='OD600', xaxt='n', yaxt='n', xlab='Hours', ylim=c(0,0.8), xlim=c(1,145))
axis(side=1, at=seq(1,176,16), labels=seq(0,20,2), cex.axis=0.8, lwd=2)
axis(side=2, at=seq(0,0.8,0.2), cex.axis=0.8, lwd=2)
abline(h=c(0.2,0.4,0.6), col='gray', lty=5, lwd=2)

draw_growth(blank, 'black')
draw_growth(consortia1, 'dodgerblue4')
draw_growth(neg_consortia, 'firebrick3')

legend('topleft', legend=c('Fresh media', 'Competitive consortia', 'Cooperative consortia'), 
       col=c('black', 'dodgerblue4', 'firebrick3'), lwd=3, cex=0.8, bg='white')
box()
dev.off()




png(filename='~/Desktop/repos/Jenior_Consortia_2022/results/negative_strains.png', units='in', width=6, height=5, res=300)
par(mar=c(3,3,1.5,0.5), mgp=c(2, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(blank[,2], type='l', ylab='OD600', xaxt='n', yaxt='n', xlab='Hours', 
     ylim=c(0,0.8), xlim=c(1,145))
axis(side=1, at=seq(1,176,16), labels=seq(0,20,2), cex.axis=0.8, lwd=2)
axis(side=2, at=seq(0,0.8,0.2), cex.axis=0.8, lwd=2)
abline(h=c(0.2,0.4,0.6), col='gray', lty=5, lwd=2)

draw_growth(blank, 'black')

draw_growth(k12, pall[1])
#draw_growth(dsm22959, pall[2])
draw_growth(lmd9, pall[3])
#draw_growth(dsm3979, pall[4])
draw_growth(dsm3979, pall[6])
draw_growth(atcc55813, pall[5])

legend('topleft', legend=c('Fresh media',
                           'Escherichia coli',
                           'Streptococcus thermophilus',
                           'Roseburia intestinalis',
                           'Bifidobacterium longum'), 
       col=c('black', pall), 
       lwd=3, cex=0.9, bg='white', text.font=c(1,3,3,3,3,3,3))
box()
dev.off()


#--------------------------------------------------------------------------------------------------__#






pdf(file='~/Desktop/repos/Jenior_FMT_2021/results/r20291_growth_consortium1.pdf', width=7, height=5)
par(mar=c(3,3,0.5,0.5), mgp=c(2, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(blank[,2], type='l', ylab='OD600', xaxt='n', yaxt='n', xlab='Hours', ylim=c(0,0.8), xlim=c(1,145))
axis(side=1, at=seq(1,176,16), labels=seq(0,20,2), cex.axis=0.8, lwd=2)
axis(side=2, at=seq(0,0.8,0.2), cex.axis=0.8, lwd=2)
abline(h=c(0.2,0.4,0.6), col='gray', lty=5, lwd=2)
polygon(c(1:169, rev(1:169)), c(blank[,3], rev(blank[,1])), col=alpha('black', 0.25), border=NA)
lines(blank[,2], lwd=3, col='black')
polygon(c(1:169, rev(1:169)), c(consortia1[,3], rev(consortia1[,1])), col=alpha('red', 0.25), border=NA)
lines(consortia1[,2], lwd=3, col='red')
legend('topleft', legend=c('Fresh media', 'Consortium 1'), 
       col=c('black','red'), lwd=3, cex=0.6, bg='white')

legend('left', legend=c('P. distasonis ATCC 8503',
                        'C. sporogenes DSM 795',
                        'L. indolis DSM 755',
                        'B. coagulans DSM 1',
                        'E. rectale ATCC 33656',
                        'B. producta DSM 2950'), 
       col=c('blue','green','orange','deeppink3','aquamarine3','cyan2'), 
       lwd=3, cex=0.6, bg='white')

polygon(c(1:169, rev(1:169)), c(atcc8503[,3], rev(atcc8503[,1])), col=alpha('blue', 0.25), border=NA)
lines(atcc8503[,2], lwd=3, col='blue')
polygon(c(1:169, rev(1:169)), c(dsm795[,3], rev(dsm795[,1])), col=alpha('green', 0.25), border=NA)
lines(dsm795[,2], lwd=3, col='green')
polygon(c(1:169, rev(1:169)), c(dsm755[,3], rev(dsm755[,1])), col=alpha('orange', 0.25), border=NA)
lines(dsm755[,2], lwd=3, col='orange')
polygon(c(1:169, rev(1:169)), c(dsm1[,3], rev(dsm1[,1])), col=alpha('deeppink3', 0.25), border=NA)
lines(dsm1[,2], lwd=3, col='deeppink3')
polygon(c(1:169, rev(1:169)), c(atcc33656[,3], rev(atcc33656[,1])), col=alpha('aquamarine3', 0.25), border=NA)
lines(atcc33656[,2], lwd=3, col='aquamarine3')
polygon(c(1:169, rev(1:169)), c(dsm2950[,3], rev(dsm2950[,1])), col=alpha('cyan2', 0.25), border=NA)
lines(dsm2950[,2], lwd=3, col='cyan2')
box()
dev.off()



# Generate figure


pdf(file='~/Desktop/repos/Jenior_FMT_2021/results/r20291_growth_spentmedia.pdf', width=5, height=4)
par(mar=c(3,3,0.5,0.5), mgp=c(2, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(fresh3[,2], type='l', ylab='OD600', xaxt='n', yaxt='n', xlab='Hours', ylim=c(0,0.8), xlim=c(1,145))
axis(side=1, at=seq(1,176,16), labels=seq(0,20,2), cex.axis=0.8, lwd=2)
axis(side=2, at=seq(0,0.8,0.2), cex.axis=0.8, lwd=2)
abline(h=c(0.2,0.4,0.6), col='gray', lty=5, lwd=2)

polygon(c(1:169, rev(1:169)), c(fresh3[,3], rev(fresh3[,1])), col=alpha('black', 0.25), border=NA)
lines(fresh3[,2], lwd=3, col='black')
polygon(c(1:169, rev(1:169)), c(r20291[,3], rev(r20291[,1])), col=alpha('purple3', 0.25), border=NA)
lines(r20291[,2], lwd=3, col='purple3')
polygon(c(1:169, rev(1:169)), c(k12[,3], rev(k12[,1])), col=alpha('darkgoldenrod2', 0.25), border=NA)
lines(k12[,2], lwd=3, col='darkgoldenrod2')
polygon(c(1:169, rev(1:169)), c(dsm933[,3], rev(dsm933[,1])), col=alpha('blue', 0.25), border=NA)
lines(dsm933[,2], lwd=3, col='blue')
polygon(c(1:169, rev(1:169)), c(dsm795[,3], rev(dsm795[,1])), col=alpha('green', 0.25), border=NA)
lines(dsm795[,2], lwd=3, col='green')
polygon(c(1:169, rev(1:169)), c(dsm109923[,3], rev(dsm109923[,1])), col=alpha('firebrick', 0.25), border=NA)
lines(dsm109923[,2], lwd=3, col='firebrick')
polygon(c(1:169, rev(1:169)), c(dsm755[,3], rev(dsm755[,1])), col=alpha('orange', 0.25), border=NA)
lines(dsm755[,2], lwd=3, col='orange')
polygon(c(1:169, rev(1:169)), c(dsm2079[,3], rev(dsm2079[,1])), col=alpha('deeppink2', 0.25), border=NA)
lines(dsm2079[,2], lwd=3, col='deeppink2')
polygon(c(1:169, rev(1:169)), c(dsm2950[,3], rev(dsm2950[,1])), col=alpha('cyan2', 0.25), border=NA)
lines(dsm2950[,2], lwd=3, col='cyan2')
polygon(c(1:169, rev(1:169)), c(dsm3979[,3], rev(dsm3979[,1])), col=alpha('darkorange3', 0.25), border=NA)
lines(dsm3979[,2], lwd=3, col='darkorange3')
polygon(c(1:169, rev(1:169)), c(dsm1[,3], rev(dsm1[,1])), col=alpha('red1', 0.25), border=NA)
lines(dsm1[,2], lwd=3, col='red1')
polygon(c(1:169, rev(1:169)), c(atcc8503[,3], rev(atcc8503[,1])), col=alpha('brown1', 0.25), border=NA)
lines(atcc8503[,2], lwd=3, col='brown1')
polygon(c(1:169, rev(1:169)), c(atcc27919[,3], rev(atcc27919[,1])), col=alpha('deepskyblue', 0.25), border=NA)
lines(atcc27919[,2], lwd=3, col='deepskyblue')
polygon(c(1:169, rev(1:169)), c(atcc55813[,3], rev(atcc55813[,1])), col=alpha('darkolivegreen1', 0.25), border=NA)
lines(atcc55813[,2], lwd=3, col='darkolivegreen1')
polygon(c(1:169, rev(1:169)), c(lmd9[,3], rev(lmd9[,1])), col=alpha('khaki2', 0.25), border=NA)
lines(lmd9[,2], lwd=3, col='khaki2')
polygon(c(1:169, rev(1:169)), c(atcc8482[,3], rev(atcc8482[,1])), col=alpha('forestgreen', 0.25), border=NA)
lines(atcc8482[,2], lwd=3, col='forestgreen')
polygon(c(1:169, rev(1:169)), c(atcc33656[,3], rev(atcc33656[,1])), col=alpha('aquamarine3', 0.25), border=NA)
lines(atcc33656[,2], lwd=3, col='aquamarine3')

legend('topleft', legend=c('Fresh media',
                           'C. aerofaciens DSM 3979',
                           'E. coli K12',
                           'B. vulgatus ATCC 8482',
                           'B. theta DSM 2079',
                           'P. distasonis ATCC 8503',
                           'B. pseudocatenulatum ATCC 27919',
                           'B. longum ATCC 55813',
                           'E. faecium DSM 109923',
                           'C. sporogenes DSM 795',
                           'L. indolis DSM 755',
                           'S. thermophilus LMD-9',
                           'E. clostridioformis DSM 933',
                           'B. coagulans DSM 1',
                           'E. rectale ATCC 33656',
                           'B. producta DSM 2950',
                           'C. difficile R20291'), 
       col=c('black','darkorange3','darkgoldenrod2','forestgreen','deeppink2','brown1','deepskyblue','darkolivegreen1','firebrick','green','orange','khaki2','blue','red1','aquamarine3','cyan2','purple3'), 
       lwd=3, cex=0.6, bg='white')
box()

dev.off()



#---------------------------------------------------------------------------------------------#

# Curve analysis

# Return important growth curve summary statistics
library(flux)
analyzeCurve <-function(data) {
  max_od <- round(max(data), digits=3) # Maximum OD
  max_rate <- round(diff(data)[which.max(diff(data))], digits=3) # Maximum growth rate
  area_under <- round(auc(c(1:length(data)), data), digits=3) # Area under curve
  summ_table <- as.data.frame(cbind(c('MaxOD','MaxRate','AUC'), c(max_od,max_rate,area_under)))
  colnames(summ_table) <- c('metric', 'value')
  return(summ_table)
}


# Calculate curve metrics
fresh <- analyzeCurve(fresh[,2])
r20291 <- analyzeCurve(r20291[,2])
dsm933 <- analyzeCurve(dsm933[,2])
dsm795 <- analyzeCurve(dsm795[,2])
dsm109923 <- analyzeCurve(dsm109923[,2])
dsm755 <- analyzeCurve(dsm755[,2])
dsm2079 <- analyzeCurve(dsm2079[,2])
dsm2950 <- analyzeCurve(dsm2950[,2])
dsm3979 <- analyzeCurve(dsm3979[,2])
dsm1 <- analyzeCurve(dsm1[,2])
atcc8503 <- analyzeCurve(atcc8503[,2])
atcc27919 <- analyzeCurve(atcc27919[,2])
atcc55813 <- analyzeCurve(atcc55813[,2])
lmd9 <- analyzeCurve(lmd9[,2])
atcc8482 <- analyzeCurve(atcc8482[,2])
k12 <- analyzeCurve(k12[,2])
atcc33656 <- analyzeCurve(atcc33656[,2])
dsm14610 <- analyzeCurve(dsm14610[,2])

# Assemble the table
growth_summary <- cbind(fresh, r20291, dsm933, dsm795, dsm109923, dsm755, dsm2079, 
                        dsm2950, dsm3979, dsm1, atcc8503, atcc27919, atcc55813, lmd9, 
                        atcc8482, k12, atcc33656, dsm14610)
strains <- c('none','C.difficileR20291','E.clostridioformisDSM933','C.sporogenesDSM795',
             'E.faeciumDSM109923','L.indolisDSM755','B.thetaiotaomicronDSM2079', 'B.productaDSM2950',
             'C.aerofaciensDSM3979','B.coagulansDSM1','P.distasonisATCC8503','B.pseudocatenulatumATCC27919',
             'B.longumATCC55813','S.thermophilusLMD9','B.vulgatusATCC8482','E.coliK12','E.rectaleATCC33656', 'R.intestinalisdsm14610')

#R intestinalis???
  

growth_summary <- as.data.frame(cbind(strains, growth_summary))
colnames(growth_summary) <- c('previous_strain', 'max_growth_rate', 'hours_to_max_rate', 'max_density', 'hours_to_max_density', 'auc')
rownames(growth_summary) <- growth_summary$metric
growth_summary$metric <- NULL
growth_summary <- as.data.frame(t(growth_summary))
colnames(growth_summary) <- c('R20291_AUC','R20291_MaxOD','R20291_MaxRate')

# Write growth summary data to supplementary table
write.table(growth_summary, file='~/Desktop/repos/Jenior_FMT_2021/data/invitro_growth/growth_summary.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE)




# Correlation analysis????










