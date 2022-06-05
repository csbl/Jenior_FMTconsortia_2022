rm(list=ls())
gc()

for (x in c('Hmisc','qgraph','corrplot','xts','igraph','randomForest','vegan')) {
  if (x %in% installed.packages()[,'Package'] == FALSE){
    install.packages(as.character(x), quiet=TRUE)
    }
  library(x, verbose=FALSE, character.only=TRUE)
}
rm(x)

# Sample variance filter for columns
colVarFilter <- function(data) {
  variances <- c()
  data <- as.data.frame(data)
  for (x in 1:ncol(data)) {variances[x] <- var(as.numeric(data[,x]))}
  minimum <- as.vector(quantile(variances))[4]
  final_data <- data[,which(variances > minimum)]
  final_data <- droplevels(final_data)
  return(final_data)
}

# Significant feature filter with random forest
diagRF <- function(training_data){
 
  # Breiman (2001). Random Forests. Machine Learning.
  lvls <- as.vector(unique(training_data$diagnosis))
  x <- round(length(rownames(training_data[which(training_data$diagnosis==lvls[1]),])) * 0.623)
  y <- round(length(rownames(training_data[which(training_data$diagnosis==lvls[2]),])) * 0.623)
  z <- max(c(round(x / y), round(y / x))) * 3
  nTrees <- round(ncol(training_data) - 1) * z
  mTries <- round(sqrt(ncol(training_data) - 1))
  diagnosis <- training_data$diagnosis
  diagnosis <- as.factor(diagnosis)
  training_data$diagnosis <- NULL
  rm(lvls, x, y, z)
  
  # Run random forest and get MDA values
  set.seed(906801)
  names(training_data) <- make.names(names(training_data))
  modelRF <- randomForest(diagnosis ~ ., data=training_data, 
                          importance=TRUE, replace=FALSE, err.rate=TRUE, 
                          ntree=nTrees, mtry=mTries)
  print(modelRF)
  featRF <- importance(modelRF, type=1)
  rm(modelRF)
  
  # Filter to significant features (Strobl 2002)
  sigFeatRF <- as.data.frame(subset(featRF, featRF > (abs(min(featRF)))))
  #sigFeatRF$metabolite <- rownames(sigFeatRF)
  #sigFeatRF <- sigFeatRF[order(-sigFeatRF[,1]),] 
  #topFeat <- round(nrow(sigFeatRF) * 0.10)
  #sigFeatRF <- sigFeatRF[1:topFeat,]
  rm(featRF)
  
  # Subset trainng data to significant features
  finalFeat <- training_data[,which(rownames(sigFeatRF) %in% colnames(training_data))]
  rm(training_data, sigFeatRF)
  
  return(finalFeat)
}

#----------------------------------------------------------------------------------------#

# Read in data
# Metabolomes
metabolome <- read.delim('/home/mjenior/Desktop/current_analyses/hmp2_metabolomics.tsv', sep='\t', header=TRUE)
# Metadata
metadata <- read.delim('/home/mjenior/Desktop/Repositories/iHMP_IBD_multiomics/data/hmp2_metadata.tsv', sep='\t', header=T, na.strings='NA')

#----------------------------------------------------------------------------------------#

# Format data
# Metabolome
metabolome <- metabolome[-which(metabolome$Metabolite == ''), ]
metabolome$Metabolite <- NULL
rownames(metabolome) <- metabolome$Method.Metabolite
metabolome$Method.Metabolite <- NULL
metabolome <- log10(metabolome + 1)

# Metadata
metadata <- subset(metadata, data_type == 'metabolomics')
keep <- c('site_sub_coll','diagnosis','Age_at_diagnosis','consent_age','fecalcal','sccai')
metadata <- metadata[, which(colnames(metadata) %in% keep)]
rownames(metadata) <- metadata$site_sub_coll
metadata$site_sub_coll <- NULL
rm(keep)

# Combine metabolome and metadata
# Subset by diagnosis
UC_metadata <- subset(metadata, diagnosis == 'UC')
nonIBD_metadata <- subset(metadata, diagnosis == 'nonIBD')
rm(metadata)
# Focus analysis on adults
UC_metadata <- subset(UC_metadata, Age_at_diagnosis >= 18)
nonIBD_metadata <- subset(nonIBD_metadata, consent_age >= 18)
# Limit by most accurate clinical readout
UC_metadata <- subset(UC_metadata, sccai >= 3) # 
nonIBD_metadata <- subset(nonIBD_metadata, fecalcal < 110)

# Prune metabolome to only samples of interest
metabolome_UC <- metabolome[, which(colnames(metabolome) %in% rownames(UC_metadata))]
temp_names <- rownames(metabolome_UC)
metabolome_UC <- sapply(metabolome_UC, as.numeric)
rownames(metabolome_UC) <- temp_names
metabolome_UC <- as.data.frame(t(metabolome_UC))
metabolome_UC$diagnosis <- rep('UC', nrow(metabolome_UC))
metabolome_nonIBD <- metabolome[, which(colnames(metabolome) %in% rownames(nonIBD_metadata))]
temp_names <- rownames(metabolome_nonIBD)
metabolome_nonIBD <- sapply(metabolome_nonIBD, as.numeric)
rownames(metabolome_nonIBD) <- temp_names
metabolome_nonIBD <- as.data.frame(t(metabolome_nonIBD))
metabolome_nonIBD$diagnosis <- rep('nonIBD', nrow(metabolome_nonIBD))
rm(UC_metadata, nonIBD_metadata, metabolome, temp_names)

# Combine for machine learning
metab_UC_nonIBD <- rbind(metabolome_UC, metabolome_nonIBD)
groups_UC_nonIBD <- as.data.frame(cbind(rownames(metab_UC_nonIBD), metab_UC_nonIBD$diagnosis))
colnames(groups_UC_nonIBD) <- c('sample', 'diagnosis')
rm(metabolome_UC, metabolome_nonIBD)

# Find most informative metabolites 
metab_UC_nonIBD_filter <- diagRF(metab_UC_nonIBD)
rm(metab_UC_nonIBD)
# OOB estimate of  error rate: 11.29%
# Confusion matrix:
#   nonIBD UC class.error
# nonIBD     37  3   0.0750000
# UC          4 18   0.1818182

# Perform ordination analysis
metabolome_nmds <- as.data.frame(metaMDS(metab_UC_nonIBD_filter, k=2, trymax=100)$points) # default = Bray
ordination <- merge(groups_UC_nonIBD, metabolome_nmds, by.x='sample', by.y='row.names')
uc_ordination <- subset(ordination, diagnosis == 'UC')
nonibd_ordination <- subset(ordination, diagnosis == 'nonIBD')




plot(x=ordination$MDS1, y=ordination$MDS2,
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
points(x=uc_ordination$MDS1, y=uc_ordination$MDS2, bg='red', pch=21, cex=2, lwd=1.2)
points(x=nonibd_ordination$MDS1, y=nonibd_ordination$MDS2, bg='blue', pch=21, cex=2, lwd=1.2)



#----------------------------------------------------------------------------------------#

# Calculate correlation matrix
corr_data <- rcorr(as.matrix(t(metab_UC_nonIBD_filter)), type='spearman')
coeff_UC_nonIBD <- corr_data$r # correlation coefficients
pval_UC_nonIBD <- corr_data$P # p-values
pval_UC_nonIBD[is.na(pval_UC_nonIBD)] <- 0
rm(corr_data)

# Prune correlations based on p-values
coeff_UC_nonIBD[which(pval_UC_nonIBD > 0.05)] <- 0.0

#----------------------------------------------------------------------------------------#

# Generate correlation network plot
pdf(file='corr_net.pdf', width=10, height=10)


par(mar=c(0,0,0,0), xpd=NA)
qgraph(coeff_UC_nonIBD, minimum=0.3, vsize=3, layout='spring', weighted=TRUE,
       groups=groups_UC_nonIBD, borders=TRUE, labels=FALSE, legend=FALSE,
       posCol='royalblue3', negCol='red3', color=c('firebrick3','white'))
text(x=-1.1, y=1.1, labels='D', font=2, cex=2)
legend('topright', legend=c('Ulcerative Colitis','Healthy'), 
       pt.bg=c('firebrick3','white'), pch=21, pt.cex=2, cex=1)
legend('bottomright', legend=c('Metabolites = 10','AUC = 0.973'), 
       pch=1, pt.cex=0, cex=0.9, bty='n', ncol=2)




dev.off()

#----------------------------------------------------------------------------------------#

#rm(list=ls())
#gc()


