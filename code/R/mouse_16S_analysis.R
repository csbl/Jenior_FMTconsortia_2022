
# Read in data
path <- '~/Desktop/mouse_16S/data'
ASVs_counts <- read.delim(paste0(path, "/ASVs_counts.tsv"), sep='\t', header=TRUE, row.names=1)
metadata <- read.delim(paste0(path, "/metadata.tsv"), sep='\t', header=TRUE, row.names=1)
taxonomy <- read.delim(paste0(path, "/ASVs_taxonomy.tsv"), sep='\t', header=TRUE, row.names=1)
rm(path)

# Subsample
library(vegan)
sub_size <- round(min(colSums(ASVs_counts))*0.9) # Determine subsample level
for (x in 1:ncol(ASVs_counts)) {ASVs_counts[,x] <- as.vector(rrarefy(ASVs_counts[,x], sample=sub_size))}
rm(sub_size, x)

# Merge
ASVs_counts <- as.data.frame(t(ASVs_counts))
ASVs_counts <- merge(metadata, ASVs_counts, by='row.names')
rownames(ASVs_counts) <- ASVs_counts$Row.names
ASVs_counts$Row.names <- NULL
ASVs_counts$mouse <- NULL
ASVs_counts$cage <- NULL

# Subset
control <- subset(ASVs_counts, group == 'control')
control$group <- NULL
competitive <- subset(ASVs_counts, group == 'competitive')
competitive$group <- NULL
cooperative <- subset(ASVs_counts, group == 'cooperative')
cooperative$group <- NULL
groups <- as.factor(ASVs_counts$group)
ASVs_counts$group <- NULL

# Calculate dissimilarity
distmat <- vegdist(ASVs_counts, method='bray')

# Mean within-group dissimilarity
meandist(distmat, grouping=groups)
#             control  competitive
# control   0.004 0.006
# competitive 0.006 0.001

# Test difference
ASVs_counts$group <- groups
pval <- adonis(distmat ~ group, data=ASVs_counts, perm=999, method='bray')
pval <- pval$aov.tab[[6]][1]
pval <- as.character(round(pval, 4))
rm(ASVs_counts)

# Ordination analysis - center point
nmds_coord <- as.data.frame(metaMDS(distmat, k=2, trymax=50)$points)
x <- (abs(max(nmds_coord$MDS1)) - abs(min(nmds_coord$MDS1))) / 2
y <- (abs(max(nmds_coord$MDS2)) - abs(min(nmds_coord$MDS2))) / 2
nmds_coord$MDS1 <- nmds_coord$MDS1 - x
nmds_coord$MDS2 <- nmds_coord$MDS2 - y
x <- max(abs(max(nmds_coord$MDS1)), abs(min(nmds_coord$MDS1))) + 0.01
y <- max(abs(max(nmds_coord$MDS2)), abs(min(nmds_coord$MDS2))) + 0.01

# Subset points
nmds_coord <- merge(x=metadata, y=nmds_coord, by='row.names')
rownames(nmds_coord) <- nmds_coord$Row.names
nmds_coord$Row.names <- NULL
control_points <- subset(nmds_coord, group == 'control')
competitive_points <- subset(nmds_coord, group == 'competitive')
cooperative_points <- subset(nmds_coord, group == 'cooperative')

# Calculate centroids
control_centroids <- aggregate(cbind(control_points$MDS1, control_points$MDS2)~control_points$group, data=control_points, mean)
competitive_centroids <- aggregate(cbind(competitive_points$MDS1, competitive_points$MDS2)~competitive_points$group, data=competitive_points, mean)
cooperative_centroids <- aggregate(cbind(cooperative_points$MDS1, cooperative_points$MDS2)~cooperative_points$group, data=cooperative_points, mean)

# Generate ordination figure
library(scales)
pdf(file='~/Desktop/Jenior_Consortia_2022/results/nmds_mouse_16S.pdf', width=4, height=4)
par(mar=c(3,3,1,1), las=1, mgp=c(2,0.6,0), lwd=2)
plot(x=nmds_coord$MDS1, y=nmds_coord$MDS2, xlim=c(-0.4,0.4), ylim=c(-0.4, 0.4),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.8)
segments(x0=control_points$MDS1, y0=control_points$MDS2, lty=5,
         x1=control_centroids[1,2], y1=control_centroids[1,3], col='ivory3')
points(x=control_points$MDS1, y=control_points$MDS2, bg=alpha('gray40',0.8), pch=21, cex=1.7)
segments(x0=competitive_points$MDS1, y0=competitive_points$MDS2, lty=5,
         x1=competitive_centroids[1,2], y1=competitive_centroids[1,3], col='ivory3')
points(x=competitive_points$MDS1, y=competitive_points$MDS2, bg=alpha('dodgerblue4',0.8), pch=21, cex=1.7)
segments(x0=cooperative_points$MDS1, y0=cooperative_points$MDS2, lty=5,
         x1=cooperative_centroids[1,2], y1=cooperative_centroids[1,3], col='ivory3')
points(x=cooperative_points$MDS1, y=cooperative_points$MDS2, bg=alpha('firebrick3',0.8), pch=21, cex=1.7)
legend('topright', legend=c('Control gavage','Competitive consortia','Cooperative consortia'), 
       pt.bg=c('gray40','dodgerblue4','firebrick3'), pch=21, pt.cex=1.6, cex=0.9, box.lwd=1.7)
legend('bottomright', legend=as.expression(bquote(paste(italic('p'),'-value = 0.001 ***'))), 
       bty='n', pt.cex=0, cex=0.8)
box()
dev.off()

rm(metadata, nmds_coord, competitive_centroids, competitive_points, cooperative_centroids, cooperative_points,
   control_centroids, control_points)

#---------------------------------------------------------------------------------#

# Combine datasets
control$condition <- 'control'
competitive$condition <- 'consortia'
cooperative$condition <- 'consortia'
gavage <- rbind(control, competitive, cooperative)
gavage$condition <- as.factor(gavage$condition)
competitive$condition <- 'competitive'
cooperative$condition <- 'cooperative'
consortia <- rbind(competitive, cooperative)
consortia$condition <- as.factor(consortia$condition)
control$condition <- 'dead'
competitive$condition <- 'dead'
cooperative$condition <- 'alive'
outcome <- rbind(control, competitive, cooperative)
outcome$condition <- as.factor(outcome$condition)
rm(control, competitive, cooperative)

# Run Random Forest and parse results
library(randomForest)
condition <- gavage$condition
gavage$condition <- NULL
rf_obj <- randomForest(condition ~ ., data=gavage, ntree=5000, mtry=15, importance=TRUE, replace=FALSE, err.rate=TRUE)
print(rf_obj)
# OOB estimate of  error rate: 33.33%
# Confusion matrix:
#         consortia control class.error
# consortia        11       3   0.2142857
# control           4       3   0.5714286
gavage_mda <- importance(rf_obj, type=1, scale=FALSE)
gavage_mda <- subset(gavage_mda, gavage_mda > (abs(min(gavage_mda)))) # significance
gavage_mda <- as.data.frame(gavage_mda)
gavage_mda$feature <- rownames(gavage_mda)
gavage_mda <- gavage_mda[order(-gavage_mda$MeanDecreaseAccuracy),]
condition <- consortia$condition
consortia$condition <- NULL
rf_obj <- randomForest(condition ~ ., data=consortia, ntree=5000, mtry=15, importance=TRUE, replace=FALSE, err.rate=TRUE)
print(rf_obj)
# OOB estimate of  error rate: 0%
# Confusion matrix:
#         competitive cooperative class.error
# competitive           7           0           0
# cooperative           0           7           0
consortia_mda <- importance(rf_obj, type=1, scale=FALSE)
consortia_mda <- subset(consortia_mda, consortia_mda > (abs(min(consortia_mda)))) # significance
consortia_mda <- as.data.frame(consortia_mda)
consortia_mda$feature <- rownames(consortia_mda)
consortia_mda <- consortia_mda[order(-consortia_mda$MeanDecreaseAccuracy),]
condition <- outcome$condition
outcome$condition <- NULL
rf_obj <- randomForest(condition ~ ., data=outcome, ntree=5000, mtry=15, importance=TRUE, replace=FALSE, err.rate=TRUE)
print(rf_obj)
# OOB estimate of  error rate: 0%
# Confusion matrix:
#         alive dead class.error
# alive     7    0           0
# dead      0   14           0
outcome_mda <- importance(rf_obj, type=1, scale=FALSE)
outcome_mda <- subset(outcome_mda, outcome_mda > (abs(min(outcome_mda)))) # significance
outcome_mda <- as.data.frame(outcome_mda)
outcome_mda$feature <- rownames(outcome_mda)
outcome_mda <- outcome_mda[order(-outcome_mda$MeanDecreaseAccuracy),]
rm(rf_obj, condition)

# Cross-reference ASVs against taxonomy
gavage_taxa <- taxonomy[gavage_mda$feature,]
consortia_taxa <- taxonomy[consortia_mda$feature,]
outcome_taxa <- taxonomy[outcome_mda$feature,]
rm(taxonomy)

# Subset ASV tables and 
gavage <- gavage[,rownames(gavage_taxa)]
gavage <- as.data.frame(t(gavage))
consortia <- consortia[,rownames(consortia_taxa)]
consortia <- as.data.frame(t(consortia))
outcome <- outcome[,rownames(outcome_taxa)]
outcome <- as.data.frame(t(outcome))

# Aggregate by genus, calculate total abundances
gavage_taxa[,c('Kingdom','Phylum','Class','Order','Family')] <- NULL
gavage <- merge(x=gavage_taxa, y=gavage, by='row.names')
gavage$Row.names <- NULL
gavage <- aggregate(. ~ Genus, data=gavage, FUN=sum)
rownames(gavage) <- gavage$Genus
gavage$Genus <- NULL
consortia_taxa[,c('Kingdom','Phylum','Class','Order','Family')] <- NULL
consortia <- merge(x=consortia_taxa, y=consortia, by='row.names')
consortia$Row.names <- NULL
consortia <- aggregate(. ~ Genus, data=consortia, FUN=sum)
rownames(consortia) <- consortia$Genus
consortia$Genus <- NULL
outcome_taxa[,c('Kingdom','Phylum','Class','Order','Family')] <- NULL
outcome <- merge(x=outcome_taxa, y=outcome, by='row.names')
outcome$Row.names <- NULL
outcome <- aggregate(. ~ Genus, data=outcome, FUN=sum)
rownames(outcome) <- outcome$Genus
outcome$Genus <- NULL

# Combine with metadata
gavage <- as.data.frame(t(gavage))
gavage <- merge(x=metadata, y=gavage, by='row.names')
rownames(gavage) <- gavage$Row.names
gavage$Row.names <- NULL
gavage$mouse <- NULL
gavage$cage <- NULL
consortia <- as.data.frame(t(consortia))
consortia <- merge(x=metadata, y=consortia, by='row.names')
rownames(consortia) <- consortia$Row.names
consortia$Row.names <- NULL
consortia$mouse <- NULL
consortia$cage <- NULL
outcome <- as.data.frame(t(outcome))
outcome <- merge(x=metadata, y=outcome, by='row.names')
rownames(outcome) <- outcome$Row.names
outcome$Row.names <- NULL
outcome$mouse <- NULL
outcome$cage <- NULL

# Subset groups by metadata
gavage_control <- subset(gavage, group == 'control')
gavage_control$group <- NULL
gavage_consortia <- subset(gavage, group %in% c('competitive','cooperative'))
gavage_consortia$group <- NULL
rm(gavage)
consortia_competitive <- subset(consortia, group == 'competitive')
consortia_competitive$group <- NULL
consortia_cooperative <- subset(consortia, group == 'cooperative')
consortia_cooperative$group <- NULL
rm(consortia)
outcome_alive <- subset(outcome, group == 'cooperative')
outcome_alive$group <- NULL
outcome_dead <- subset(outcome, group %in% c('competitive','control'))
outcome_dead$group <- NULL
rm(outcome)

# Transform abundances
gavage_control <- log10(gavage_control + 1)
gavage_consortia <- log10(gavage_consortia + 1)
consortia_competitive <- log10(consortia_competitive + 1)
consortia_cooperative <- log10(consortia_cooperative + 1)
outcome_alive <- log10(outcome_alive + 1)
outcome_dead <- log10(outcome_dead + 1)

# Rank by median absolute differences
abs_diffs <- c()
for (x in 1:ncol(gavage_control)) {abs_diffs[x] <- abs(median(gavage_control[,x]) - median(gavage_consortia[,x]))}
gavage_control <- gavage_control[,order(abs_diffs)]
gavage_consortia <- gavage_consortia[,order(abs_diffs)]
abs_diffs <- c()
for (x in 1:ncol(consortia_competitive)) {abs_diffs[x] <- abs(median(consortia_competitive[,x]) - median(consortia_cooperative[,x]))}
consortia_competitive <- consortia_competitive[,order(abs_diffs)]
consortia_cooperative <- consortia_cooperative[,order(abs_diffs)]
abs_diffs <- c()
for (x in 1:ncol(outcome_alive)) {abs_diffs[x] <- abs(median(outcome_alive[,x]) - median(outcome_dead[,x]))}
outcome_alive <- outcome_alive[,order(abs_diffs)]
outcome_dead <- outcome_dead[,order(abs_diffs)]
rm(abs_diffs)

# Generate figures
asvPlot <- function(group1, group2, oob) {
        xmax <- max(c(max(group1), max(group2)))
        ymax <- (ncol(group1)*3) + 0.5
        genera <- gsub('_', ' ', colnames(group1))
        
        par(mar=c(3, 0.5, 1.5, 0.5), mgp=c(1.8, 0.5, 0), xpd=FALSE, yaxs='i', lwd=2)
        plot(1, type='n', ylim=c(0.5,ymax), xlim=c(0,xmax), 
             ylab='', xlab='Relative Abundance (Log10)', yaxt='n', cex.lab=0.9)
        index <- 2
        for(i in c(1:ncol(group1))){
                stripchart(at=index+0.4, group1[,i], 
                           pch=21, bg='white', method='jitter', jitter=0.12, cex=1.5, add=TRUE)
                stripchart(at=index-0.8, group2[,i], 
                           pch=21, bg='gray40', method='jitter', jitter=0.12, cex=1.5, add=TRUE)
                if (i != ncol(group1)){abline(h=index+1.5, lty=2)}
                text(x=0, y=index+0.9, labels=genera[i], cex=0.75, font=3, pos=4)
                segments(median(group1[,i]), index+0.8, median(group1[,i]), index, lwd=2)
                segments(median(group2[,i]), index-1.2, median(group2[,i]), index-0.4, lwd=2)
                index <- index + 3}
        mtext(paste0('OOB = ', oob, '%'), side=3, cex=0.75, padj=-1, adj=1)}


#layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow=TRUE))
#asvPlot(gavage_control, gavage_consortia, 33.33)
#asvPlot(consortia_competitive, consortia_cooperative, 0)
pdf(file='~/Desktop/Jenior_Consortia_2022/results/outcome_features_16S.pdf', width=3, height=6)
asvPlot(outcome_alive, outcome_dead, 0.0)
legend('bottomleft', legend=c('Survived','Deceased'), pt.bg=c('white', 'gray40'), 
       pch=21, cex=0.8, pt.cex=1.5, bty='n', horiz=TRUE)
box()
dev.off()

#---------------------------------------------------------------------------------#


# Beta-diversity
control_div <- as.vector(diversity(control, index='invsimpson'))
competitive_div <- as.vector(diversity(competitive, index='invsimpson'))
cooperative_div <- as.vector(diversity(cooperative, index='invsimpson'))
pvals <- p.adjust(c(round(wilcox.test(control_div, competitive_div, exact=FALSE)$p.value, 4),
                    round(wilcox.test(control_div, cooperative_div, exact=FALSE)$p.value, 4),
                    round(wilcox.test(cooperative_div, competitive_div, exact=FALSE)$p.value, 4)), method='BH')
cont_comp_pval <- pvals[1]
cont_coop_pval <- pvals[2]
coop_comp_pval <- pvals[3]
rm(pvals)


y_max <- ceiling(max(c(control_div,competitive_div,cooperative_div)))

library(plotrix)

pdf(file='~/Desktop/Jenior_Consortia_2022/results/div_mouse_16S.pdf', width=2.5, height=4)
par(mar=c(3.5,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(20,95), xlim=c(0,3), 
     ylab='Inv. Simpson Diversity', xlab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(control_div, cex=0, lwd=1.8, at=0.5, col='gray40', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
boxplot(competitive_div, cex=0, lwd=1.8, at=1.5, col='dodgerblue4',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
boxplot(cooperative_div, cex=0, lwd=1.8, at=2.5, col='firebrick3',ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=1.8, xaxt='n', yaxt='n', add=TRUE)
axis(side=2, at=seq(20,95,15), labels=c(0,seq(35,95,15)), cex.axis=0.8)
axis.break(2, 25, style='slash', brw=0.04)
segments(x0=c(0.5,0.5,1.5), y0=c(80,85,90), x1=c(1.5,1.5,2.5))
segments(x0=0.5, y0=85, x1=2.5)
text(x=1, y=82, 'n.s.', cex=0.6)
text(x=c(1.5,2), y=c(87,92), '**', font=2, cex=1.1)

par(xpd=TRUE)
text(x=c(0.5,1.5,2.5), y=8, srt=40, cex=0.8, 
     labels=c('Control\ngavage','Competitive\nconsortia','Cooperative\nconsortia'))
dev.off()




