
# Completeness analysis
file_prefix <- '~/Desktop/repos/Jenior_Consortia_2022/data/pathway_analysis/completeness/completeness_'
s01044A <- read.delim(paste0(file_prefix,'01044A','.tsv'), sep='\t', header=TRUE)
s01090A <- read.delim(paste0(file_prefix,'01090A','.tsv'), sep='\t', header=TRUE)
s01092A <- read.delim(paste0(file_prefix,'01092A','.tsv'), sep='\t', header=TRUE)
s01093A <- read.delim(paste0(file_prefix,'01093A','.tsv'), sep='\t', header=TRUE)
n05042G <- read.delim(paste0(file_prefix,'05042G','.tsv'), sep='\t', header=TRUE)
n05053B <- read.delim(paste0(file_prefix,'05053B','.tsv'), sep='\t', header=TRUE)
n05080M <- read.delim(paste0(file_prefix,'05080M','.tsv'), sep='\t', header=TRUE)
n05098A <- read.delim(paste0(file_prefix,'05098A','.tsv'), sep='\t', header=TRUE)
n28043C <- read.delim(paste0(file_prefix,'28043C','.tsv'), sep='\t', header=TRUE)
n28045A <- read.delim(paste0(file_prefix,'28045A','.tsv'), sep='\t', header=TRUE)
n28045D <- read.delim(paste0(file_prefix,'28045D','.tsv'), sep='\t', header=TRUE)
n28047D <- read.delim(paste0(file_prefix,'28047D','.tsv'), sep='\t', header=TRUE)
super <- merge(s01044A, s01090A, by='pathway')
super <- merge(super, s01092A, by='pathway')
super <- merge(super, s01093A, by='pathway')
normal <- merge(n05042G, n05053B, by='pathway')
normal <- merge(normal, n05080M, by='pathway')
normal <- merge(normal, n05098A, by='pathway')
normal <- merge(normal, n28043C, by='pathway')
normal <- merge(normal, n28045A, by='pathway')
normal <- merge(normal, n28045D, by='pathway')
normal <- merge(normal, n28047D, by='pathway')
success <- merge(s01044A, s01090A, by='pathway')
success <- merge(success, s01092A, by='pathway')
success <- merge(success, s01093A, by='pathway')
success <- merge(success, n05098A, by='pathway')
success <- merge(success, n05080M, by='pathway')
success <- merge(success, n28043C, by='pathway')
success <- merge(success, n28045D, by='pathway')
failure <- merge(n05042G, n05053B, by='pathway')
failure <- merge(failure, n28045A, by='pathway')
failure <- merge(failure, n28047D, by='pathway')
rm(s01044A, s01090A, s01092A, s01093A)
rm(n05042G, n05053B, n05080M, n05098A)
rm(n28043C, n28045A, n28045D, n28047D)

library(Matrix)
completeness <- merge(super, normal, by='pathway')
rownames(completeness) <- completeness$pathway
completeness$pathway <- NULL
metabolic_paths <- c("2-Oxocarboxylic_acid_metabolism","ABC_transporters","Alanine,_aspartate_and_glutamate_metabolism",
                     "Amino_sugar_and_nucleotide_sugar_metabolism","Arginine_and_proline_metabolism","Ascorbate_and_aldarate_metabolism",
                     "beta-Alanine_metabolism","Carbohydrate_digestion_and_absorption",
                     "Carbon_metabolism","Citrate_cycle_(TCA_cycle)","Cyanoamino_acid_metabolism","Cysteine_and_methionine_metabolism",
                     "D-Alanine_metabolism","D-Arginine_and_D-ornithine_metabolism","D-Glutamine_and_D-glutamate_metabolism","Degradation_of_aromatic_compounds",
                     "Fatty_acid_degradation","Fatty_acid_metabolism","Fructose_and_mannose_metabolism","Galactose_metabolism",
                     "Glutathione_metabolism","Glycerolipid_metabolism","Glycerophospholipid_metabolism","Glycine,_serine_and_threonine_metabolism","Glycolysis_Gluconeogenesis",
                     "Glycosaminoglycan_degradation","Histidine_metabolism","Lysine_degradation",
                     "Microbial_metabolism_in_diverse_environments","Mineral_absorption",
                     "Nicotinate_and_nicotinamide_metabolism","Nitrogen_metabolism",
                     "Pentose_and_glucuronate_interconversions","Pentose_phosphate_pathway","Phenylalanine_metabolism",
                     "Phosphonate_and_phosphinate_metabolism","Phosphotransferase_system_(PTS)",
                     "Polycyclic_aromatic_hydrocarbon_degradation","Propanoate_metabolism","Protein_digestion_and_absorption",
                     "Purine_metabolism","Pyrimidine_metabolism","Pyruvate_metabolism","Retinol_metabolism","Riboflavin_metabolism",
                     "Selenocompound_metabolism","Sphingolipid_metabolism","Starch_and_sucrose_metabolism",
                     "Sulfur_metabolism","Synthesis_and_degradation_of_ketone_bodies","Taurine_and_hypotaurine_metabolism",
                     "Thiamine_metabolism","Toluene_degradation","Tryptophan_metabolism","Tyrosine_metabolism",
                     "Valine,_leucine_and_isoleucine_degradation","Vitamin_B6_metabolism",
                     "Vitamin_digestion_and_absorption","Xylene_degradation")
completeness <- completeness[metabolic_paths,]
rownames(super) <- super$pathway
super$pathway <- NULL
super <- super[metabolic_paths,]
rownames(normal) <- normal$pathway
normal$pathway <- NULL
normal <- normal[metabolic_paths,]
rownames(success) <- success$pathway
success$pathway <- NULL
success <- success[metabolic_paths,]
rownames(failure) <- failure$pathway
failure$pathway <- NULL
failure <- failure[metabolic_paths,]
completeness <- as.data.frame(t(completeness))
colnames(completeness) <- make.names(colnames(completeness))
keep <- c()
for (x in colnames(completeness)) {
  if (nnzero(completeness[,x]) >= 4) {
    keep <- c(keep, x)}}
rm(x)
completeness <- completeness[,keep]
super <- as.data.frame(t(super))
colnames(super) <- make.names(colnames(super))
super <- super[,keep]
super <- droplevels(super)
normal <- as.data.frame(t(normal))
colnames(normal) <- make.names(colnames(normal))
normal <- normal[,keep]
normal <- droplevels(normal)
success <- as.data.frame(t(success))
colnames(success) <- make.names(colnames(success))
success <- success[,keep]
success <- droplevels(success)
failure <- as.data.frame(t(failure))
colnames(failure) <- make.names(colnames(failure))
failure <- failure[,keep]
failure <- droplevels(failure)
rm(keep)

# Top pathway differences from abundance analysis
donor <- c("Amino_sugar_and_nucleotide_sugar_metabolism","Arginine_and_proline_metabolism",
           "Ascorbate_and_aldarate_metabolism","Cysteine_and_methionine_metabolism",
           "Fructose_and_mannose_metabolism","Glutathione_metabolism","Glycosaminoglycan_degradation",
           "Nitrogen_metabolism","Protein_digestion_and_absorption","Purine_metabolism",
           "Retinol_metabolism","Riboflavin_metabolism","Starch_and_sucrose_metabolism",
           "Synthesis_and_degradation_of_ketone_bodies","Valine._leucine_and_isoleucine_degradation",
           "Vitamin_B6_metabolism")
outcome <- c("Arginine_and_proline_metabolism","Cyanoamino_acid_metabolism","Glycerolipid_metabolism",
             "Glycosaminoglycan_degradation","Nitrogen_metabolism","Protein_digestion_and_absorption")
super <- super[,donor]
normal <- normal[,donor]
success <- success[,outcome]
failure <- failure[,outcome]
rm(donor, outcome)


wilcox.test(super[,'Arginine_and_proline_metabolism'], normal[,'Arginine_and_proline_metabolism'], exact=FALSE)
wilcox.test(success[,'Arginine_and_proline_metabolism'], failure[,'Arginine_and_proline_metabolism'], exact=FALSE)


wilcox.test(super[,'Glycosaminoglycan_degradation'], normal[,'Glycosaminoglycan_degradation'], exact=FALSE)
wilcox.test(success[,'Glycosaminoglycan_degradation'], failure[,'Glycosaminoglycan_degradation'], exact=FALSE)





library(randomForest)
metadata <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/metadata.tsv', sep='\t', header=TRUE)
donor <- as.factor(metadata$status)
outcome <- as.factor(metadata$outcome)
rm(metadata)

rf_obj <- randomForest(donor ~ ., data=completeness, ntree=10000, mtry=50,
                       importance=TRUE, replace=FALSE, err.rate=TRUE)
# OOB estimate of  error rate: 0%
# Confusion matrix:
#           normal super class.error
# normal      8     0           0
# super       0     4           0
rf_feat <- importance(rf_obj, type=1, scale=FALSE)
super_rf_importance <- subset(rf_feat, rf_feat > abs(min(rf_feat)))
super_rf_importance <- as.data.frame(super_rf_importance)
super_rf_importance$feature <- rownames(super_rf_importance)
colnames(super_rf_importance) <- c('mda', 'pathway')
super_rf_importance <- super_rf_importance[order(-super_rf_importance$mda),]
rm(rf_obj, rf_feat)


super_completeness <- completeness[,super_rf_importance$pathway]
super_completeness <- droplevels(super_completeness)
super_completeness <- as.data.frame(t(super_completeness))
super <- super_completeness[,c("X01044A","X01090A","X01092A","X01093A")]
super <- as.data.frame(t(super))
normal <- super_completeness[,c("X05042G","X05053B","X05080M","X05098A","X28043C","X28045A","X28045D","X28047D")]
normal <- as.data.frame(t(normal))
rm(super_completeness)


super_pvals <- c()
for (x in 1:ncol(super)) {super_pvals[x] <- wilcox.test(super[,x], normal[,x], exact=FALSE)$p.value}
super_pvals <- p.adjust(super_pvals, method='BH')


super <- as.data.frame(apply(super, MARGIN=2, median))
normal <- as.data.frame(apply(normal, MARGIN=2, median))
donor <- merge(super, normal, by='row.names')
rownames(donor) <- donor$Row.names
donor$Row.names <- NULL
colnames(donor) <- c('super', 'normal')

# Calculate percentages
donor$super_perc <- (donor$super / (donor$super + donor$normal)) * 100.0
donor$normal_perc <- (donor$normal / (donor$super + donor$normal)) * 100.0
donor$super <- NULL
donor$normal <- NULL
colnames(donor) <- c('super', 'normal')

# Radar plot tables
donor <- as.data.frame(t(donor))
colnames(donor) <- gsub('_', ' ', colnames(donor))
colnames(donor) <- gsub('\\.', ', ', colnames(donor))
donor <- rbind(rep(ceiling(max(donor)), ncol(donor)) , rep(floor(min(donor)), ncol(donor)) , donor)
write.table(donor, file='~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar.tsv', append = FALSE, quote=FALSE, sep='\t')
donor <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/donor_radar.tsv', sep='\t', header=TRUE)
labs <- gsub('_', ' ', colnames(donor))
labs <- gsub('\\.', '\n', labs)

# Generate figures
library(fmsb)
library(scales)
super_col <- '#8e7cc3'
normal_col <- '#e69138'

#png(file='~/Desktop/repos/Jenior_FMT_2021/results/donor_pathway_radar.png', units='in', width=5, height=5, res=300)
par(mar=c(0,2,0,2), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(donor, pcol=c(super_col,normal_col), pfcol=c(alpha(super_col,0.7),alpha(normal_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.6, cglcol='grey30', cglty=5, cglwd=1, vlabels=labs)
#text(x=0, y=0, 'OOB = 8.33%', font=2, cex=0.8)
legend('bottomright', legend=c('Super Donor', 'Normal Donors'), col=c(super_col,normal_col),
       pt.bg=c(alpha(super_col,0.7),alpha(normal_col,0.7)), pch=22, pt.cex=1.5, cex=0.8, bty='n')
#dev.off()


#------#

rf_obj <- randomForest(outcome ~ ., data=completeness, ntree=10000, mtry=50,
importance=TRUE, replace=FALSE, err.rate=TRUE)
# OOB estimate of  error rate: 41.67%
# Confusion matrix:
#  failure success class.error
# failure       1       3        0.75
# success       2       6        0.25
rf_feat <- importance(rf_obj, type=1, scale=FALSE)
success_rf_importance <- subset(rf_feat, rf_feat > abs(min(rf_feat)))
success_rf_importance <- as.data.frame(success_rf_importance)
success_rf_importance$feature <- rownames(success_rf_importance)
colnames(success_rf_importance) <- c('mda', 'pathway')
success_rf_importance <- success_rf_importance[order(-success_rf_importance$mda),]
rm(rf_obj, rf_feat)


success_completeness <- completeness[,success_rf_importance$pathway]
success_completeness <- droplevels(success_completeness)
success_completeness <- as.data.frame(t(success_completeness))


success <- success_completeness[,c("X01044A","X01090A","X01092A","X01093A","X05080M","X05098A","X28043C","X28045D")]
success <- as.data.frame(t(success))
failure <- success_completeness[,c("X05042G","X05053B","X28045A","X28047D")]
failure <- as.data.frame(t(failure))
rm(success_completeness)


success_pvals <- c()
for (x in 1:ncol(success)) {success_pvals[x] <- wilcox.test(success[,x], failure[,x], exact=FALSE)$p.value}
success_pvals <- p.adjust(success_pvals, method='BH')


success <- as.data.frame(apply(success, MARGIN=2, median))
failure <- as.data.frame(apply(failure, MARGIN=2, median))
outcome <- merge(success, failure, by='row.names')
rownames(outcome) <- outcome$Row.names
outcome$Row.names <- NULL
colnames(outcome) <- c('success', 'failure')

# Calculate percentages
outcome$success_perc <- (outcome$success / (outcome$success + outcome$failure)) * 100.0
outcome$failure_perc <- (outcome$failure / (outcome$success + outcome$failure)) * 100.0
outcome$success <- NULL
outcome$failure <- NULL
colnames(outcome) <- c('success', 'failure')

# Radar plot tables
outcome <- as.data.frame(t(outcome))
colnames(outcome) <- gsub('_', ' ', colnames(outcome))
colnames(outcome) <- gsub('\\.', ', ', colnames(outcome))
outcome <- rbind(rep(ceiling(max(outcome)), ncol(outcome)) , rep(floor(min(outcome)), ncol(outcome)) , outcome)
write.table(outcome, file='~/Desktop/repos/Jenior_Consortia_2022/data/outcome_radar.tsv', append = FALSE, quote=FALSE, sep='\t')
outcome <- read.delim('~/Desktop/repos/Jenior_Consortia_2022/data/outcome_radar.tsv', sep='\t', header=TRUE)
labs <- gsub('_', ' ', colnames(outcome))
labs <- gsub('\\.', '\n', labs)

# Generate figures
library(fmsb)
library(scales)
success_col <- '#6aa84f'
failure_col <- '#990000'

#png(file='~/Desktop/repos/Jenior_FMT_2021/results/outcome_pathway_radar.png', units='in', width=5, height=5, res=300)
par(mar=c(0,2,0,2), xpd=TRUE, lwd=2, las=1, font=2)
radarchart(outcome, pcol=c(success_col,failure_col), pfcol=c(alpha(success_col,0.7),alpha(failure_col,0.7)), 
           plwd=2 , plty=1, vlcex=0.6, cglcol='grey30', cglty=5, cglwd=1, vlabels=labs)
#text(x=0, y=0, 'OOB = 8.33%', font=2, cex=0.8)
legend('bottomright', legend=c('success outcome', 'failure outcomes'), col=c(success_col,failure_col),
       pt.bg=c(alpha(success_col,0.7),alpha(failure_col,0.7)), pch=22, pt.cex=1.5, cex=0.8, bty='n')
#dev.off()


