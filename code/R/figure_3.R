
# MAG abundance analyses
readAbund <- function(file) {
  abund <- read.delim(file, sep='\t', header=FALSE)
  donor <- gsub('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.','',file)
  donor <- gsub('.binAbund.tsv','',donor)
  colnames(abund) <- c('mag','abund')
  x <- paste0('fmt.metaG.', donor, '.bin.')
  abund$mag <- gsub(x,'',abund$mag)
  abund$mag <- gsub('.genes.pep.60.fasta','',abund$mag)
  abund$sample <- rep(donor, nrow(abund))
  return(abund)}

# CheckM screening data
readCheckMTab <- function(sample_name, screen=TRUE) {
  path <- '/home/mjenior/Desktop/active_projects/CDI_FMT_project/fmt_metaG/checkM_results/'
  sample <- paste0(path, sample_name)
  sample <- read.delim(sample, header=TRUE, sep='\t')
  sample_name <- gsub('.checkM.all.tsv','',sample_name)
  sample$sample <- rep(sample_name, nrow(sample))
  sample$bin <- gsub('bin.','',sample$bin)
  sample$bin <- gsub('.format.pruned','',sample$bin)
  
  if (screen == TRUE) {
    sample <- subset(sample, completeness > 50)
    sample <- subset(sample, contamination < 33)
    sample <- subset(sample, genes >= 1354)
    sample <- subset(sample, genome_size >= 1300000)}
  return(sample)}

# MAG consolidation
consolidateMAGs <- function(mags, sample) {
  mags$species <- rep('unknown', nrow(mags))
  path <- '/home/mjenior/Desktop/Jenior_Consortia_2022/data/kegg_hits/'
  hits <- paste0(path, sample, '.top_kegg_hits.reformat.tsv')
  hits <- read.delim(hits, header=TRUE, sep='\t')
  hits <- subset(hits, percentage >= 40.0)
  
  for (species in unique(hits$species_call)) {
    curr <- subset(hits, species_call == species)
    sub_mags <- subset(mags, mag %in% curr$bin)
    sub_abund <- median(sub_mags$abund)
    if (nrow(sub_mags) > 1) {sub_mags <- sub_mags[1,]}
    sub_mags$species <- rep(species, nrow(sub_mags))
    sub_mags$abund <- rep(sub_abund, nrow(sub_mags))
    mags <- subset(mags, !mag %in% curr$bin)
    mags <- as.data.frame(rbind(mags, sub_mags))}
  
  colnames(mags) <- c('mag','med_abund','sample','species')
  mags <- subset(mags, species != 'unknown')
  return(mags)}

# Species consolidation
consolidateSpecies <- function(mags) {
  
  for (x in unique(mags$species)) {
    sub_mags <- subset(mags, species == x)
    mags <- subset(mags, species != x)
    sub_abund <- median(sub_mags$med_abund)
    if (nrow(sub_mags) > 1) {sub_mags <- sub_mags[1,]}
    sub_mags$med_abund <- rep(sub_abund, nrow(sub_mags))
    mags <- as.data.frame(rbind(mags, sub_mags))}
  
  mags$med_abund <- round(log10(mags$med_abund), 2)
  mags$mag <- NULL
  
  return(mags)}

# Add outcome ta
addta <- function(mags, num=TRUE) {
  mags$outcome <- mags$sample
  mags$outcome <- gsub('01044A', 'success', mags$outcome)
  mags$outcome <- gsub('01090A', 'success', mags$outcome)
  mags$outcome <- gsub('01093A', 'success', mags$outcome)
  mags$outcome <- gsub('01092A', 'success', mags$outcome)
  mags$outcome <- gsub('05098A', 'success', mags$outcome)
  mags$outcome <- gsub('05080M', 'success', mags$outcome)
  mags$outcome <- gsub('28043C', 'success', mags$outcome)
  mags$outcome <- gsub('28045D', 'success', mags$outcome)
  mags$outcome <- gsub('05042G', 'failure', mags$outcome)
  mags$outcome <- gsub('05053B', 'failure', mags$outcome)
  mags$outcome <- gsub('28045A', 'failure', mags$outcome)
  mags$outcome <- gsub('28047D', 'failure', mags$outcome)
  
  mags$donor <- mags$sample
  mags$donor <- gsub('01044A', 'super', mags$donor)
  mags$donor <- gsub('01090A', 'super', mags$donor)
  mags$donor <- gsub('01093A', 'super', mags$donor)
  mags$donor <- gsub('01092A', 'super', mags$donor)
  mags$donor <- gsub('05098A', 'normal', mags$donor)
  mags$donor <- gsub('05080M', 'normal', mags$donor)
  mags$donor <- gsub('28043C', 'normal', mags$donor)
  mags$donor <- gsub('28045D', 'normal', mags$donor)
  mags$donor <- gsub('05042G', 'normal', mags$donor)
  mags$donor <- gsub('05053B', 'normal', mags$donor)
  mags$donor <- gsub('28045A', 'normal', mags$donor)
  mags$donor <- gsub('28047D', 'normal', mags$donor)
  
  mags$sample <- NULL
  
  if (num == TRUE) {
    mags$outcome <- gsub('success', 1, mags$outcome)
    mags$outcome <- gsub('failure', 0, mags$outcome)
    mags$outcome <- as.numeric(mags$outcome)
    mags$donor <- gsub('super', 1, mags$donor)
    mags$donor <- gsub('normal', 0, mags$donor)
    mags$donor <- as.numeric(mags$donor)
  }
  return(mags)}

#---------------------------------#

# Read in data
abund_01044A <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.01044A.binAbund.tsv')
abund_01090A <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.01090A.binAbund.tsv')
abund_01092A <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.01092A.binAbund.tsv')
abund_01093A <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.01093A.binAbund.tsv')
abund_05042G <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.05042G.binAbund.tsv')
abund_05053B <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.05053B.binAbund.tsv')
abund_05080M <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.05080M.binAbund.tsv')
abund_05098A <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.05098A.binAbund.tsv')
abund_28043C <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.28043C.binAbund.tsv')
abund_28045A <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.28045A.binAbund.tsv')
abund_28045D <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.28045D.binAbund.tsv')
abund_28047D <- readAbund('~/Desktop/Jenior_Consortia_2022/data/bin_abundances/fmt.metaG.28047D.binAbund.tsv')
rm(readAbund)

# Prune mags based on checkM screens
checkM_01044A <- readCheckMTab('01044A.checkM.tsv')
abund_01044A <- subset(abund_01044A, mag %in% checkM_01044A$bin)
checkM_01090A <- readCheckMTab('01090A.checkM.tsv')
abund_01090A <- subset(abund_01090A, mag %in% checkM_01090A$bin)
checkM_01092A <- readCheckMTab('01092A.checkM.tsv')
abund_01092A <- subset(abund_01092A, mag %in% checkM_01092A$bin)
checkM_01093A <- readCheckMTab('01093A.checkM.tsv')
abund_01093A <- subset(abund_01093A, mag %in% checkM_01093A$bin)
rm(checkM_01044A, checkM_01090A, checkM_01092A, checkM_01093A)
checkM_05053B <- readCheckMTab('05053B.checkM.tsv')
abund_05053B <- subset(abund_05053B, mag %in% checkM_05053B$bin)
checkM_05080M <- readCheckMTab('05080M.checkM.tsv')
abund_05080M <- subset(abund_05080M, mag %in% checkM_05080M$bin)
checkM_05042G <- readCheckMTab('05042G.checkM.tsv')
abund_05042G <- subset(abund_05042G, mag %in% checkM_05042G$bin)
checkM_05098A <- readCheckMTab('05098A.checkM.tsv')
abund_05098A <- subset(abund_05098A, mag %in% checkM_05098A$bin)
rm(checkM_05053B, checkM_05080M, checkM_05042G, checkM_05098A)
checkM_28045D <- readCheckMTab('28045D.checkM.tsv')
abund_28045D <- subset(abund_28045D, mag %in% checkM_28045D$bin)
checkM_28043C <- readCheckMTab('28043C.checkM.tsv')
abund_28043C <- subset(abund_28043C, mag %in% checkM_28043C$bin)
checkM_28047D <- readCheckMTab('28047D.checkM.tsv')
abund_28047D <- subset(abund_28047D, mag %in% checkM_28047D$bin)
checkM_28045A <- readCheckMTab('28045A.checkM.tsv')
abund_28045A <- subset(abund_28045A, mag %in% checkM_28045A$bin)
rm(checkM_28045D, checkM_28043C, checkM_28047D, checkM_28045A)
rm(readCheckMTab)

# Consolidate duplicate species within samples
abund_01044A <- consolidateMAGs(abund_01044A, '01044A')
abund_01090A <- consolidateMAGs(abund_01090A, '01090A')
abund_01092A <- consolidateMAGs(abund_01092A, '01092A')
abund_01093A <- consolidateMAGs(abund_01093A, '01093A')
abund_05042G <- consolidateMAGs(abund_05042G, '05042G')
abund_05053B <- consolidateMAGs(abund_05053B, '05053B')
abund_05080M <- consolidateMAGs(abund_05080M, '05080M')
abund_05098A <- consolidateMAGs(abund_05098A, '05098A')
abund_28043C <- consolidateMAGs(abund_28043C, '28043C')
abund_28045A <- consolidateMAGs(abund_28045A, '28045A')
abund_28045D <- consolidateMAGs(abund_28045D, '28045D')
abund_28047D <- consolidateMAGs(abund_28047D, '28047D')
rm(consolidateMAGs)

# Add ta
abund_01044A <- addMetadata(abund_01044A)
abund_01090A <- addMetadata(abund_01090A)
abund_01092A <- addMetadata(abund_01092A)
abund_01093A <- addMetadata(abund_01093A)
abund_05042G <- addMetadata(abund_05042G)
abund_05053B <- addMetadata(abund_05053B)
abund_05080M <- addMetadata(abund_05080M)
abund_05098A <- addMetadata(abund_05098A)
abund_28043C <- addMetadata(abund_28043C)
abund_28045A <- addMetadata(abund_28045A)
abund_28045D <- addMetadata(abund_28045D)
abund_28047D <- addMetadata(abund_28047D)
rm(addMetadata)

# Combine groups of interest
super_abund <- as.data.frame(rbind(abund_01044A, abund_01090A))
super_abund <- as.data.frame(rbind(super_abund, abund_01092A))
super_abund <- as.data.frame(rbind(super_abund, abund_01093A))
super_abund$outcome <- NULL
normal_abund <- as.data.frame(rbind(abund_05042G, abund_05053B))
normal_abund <- as.data.frame(rbind(normal_abund, abund_05080M))
normal_abund <- as.data.frame(rbind(normal_abund, abund_05098A))
normal_abund <- as.data.frame(rbind(normal_abund, abund_28043C))
normal_abund <- as.data.frame(rbind(normal_abund, abund_28045A))
normal_abund <- as.data.frame(rbind(normal_abund, abund_28045D))
normal_abund <- as.data.frame(rbind(normal_abund, abund_28047D))
normal_abund$outcome <- NULL
success_abund <- as.data.frame(rbind(abund_01044A, abund_01090A))
success_abund <- as.data.frame(rbind(success_abund, abund_01092A))
success_abund <- as.data.frame(rbind(success_abund, abund_01093A))
success_abund <- as.data.frame(rbind(success_abund, abund_05098A))
success_abund <- as.data.frame(rbind(success_abund, abund_05080M))
success_abund <- as.data.frame(rbind(success_abund, abund_28043C))
success_abund <- as.data.frame(rbind(success_abund, abund_28045D))
success_abund$donor <- NULL
failure_abund <- as.data.frame(rbind(abund_05042G, abund_05053B))
failure_abund <- as.data.frame(rbind(failure_abund, abund_28045A))
failure_abund <- as.data.frame(rbind(failure_abund, abund_28047D))
failure_abund$donor <- NULL
#abund_01 <- as.data.frame(rbind(abund_01044A,abund_01090A,abund_01092A,abund_01093A))
rm(abund_01044A,abund_01090A,abund_01092A,abund_01093A)
#abund_05 <- as.data.frame(rbind(abund_05042G,abund_05053B,abund_05080M,abund_05098A))
rm(abund_05042G,abund_05053B,abund_05080M,abund_05098A)
#abund_28 <- as.data.frame(rbind(abund_28043C,abund_28045A,abund_28045D,abund_28047D))
rm(abund_28043C,abund_28045A,abund_28045D,abund_28047D)

# Consolidate duplicate species within donors
super_abund <- consolidateSpecies(super_abund)
normal_abund <- consolidateSpecies(normal_abund)
success_abund <- consolidateSpecies(success_abund)
failure_abund <- consolidateSpecies(failure_abund)
rm(consolidateSpecies)

# Combine data
donor_shared <- super_abund$species[which(super_abund$species %in% normal_abund$species)]
donor_super_only <- super_abund$species[which(!super_abund$species %in% normal_abund$species)]
donor_normal_only <- normal_abund$species[which(!normal_abund$species %in% super_abund$species)]
outcome_success_only <- success_abund$species[which(!success_abund$species %in% failure_abund$species)]
outcome_failure_only <- failure_abund$species[which(!failure_abund$species %in% success_abund$species)]
success_normal <- outcome_success_only[which(outcome_success_only %in% normal_abund$species)]
success_super <- outcome_success_only[which(outcome_success_only %in% super_abund$species)]
temp1 <- subset(super_abund, species %in% donor_shared)
temp2 <- subset(normal_abund, species %in% donor_shared)
donor_abund <- merge(temp1, temp2, by='species')
donor_abund[,c('donor.x','donor.y')] <- NULL
colnames(donor_abund) <- c('species','med_abund_super','med_abund_normal')
temp1 <- subset(success_abund, species %in% outcome_shared)
temp2 <- subset(failure_abund, species %in% outcome_shared)
outcome_abund <- merge(temp1, temp2, by='species')
outcome_abund[,c('outcome.x','outcome.y')] <- NULL
colnames(donor_abund) <- c('species','med_abund_success','med_abund_failure')
rm(temp1, temp2)

# Get just genus names
getGenera <- function(data) {
  for (x in 1:length(data)) {data[x] <- unlist(strsplit(data[x], '_'))[1]}
  data <- gsub('Anaerobutyricum', 'Eubacterium', data)
  data <- gsub('Bacteroides', 'Phocaeicola', data)
  data <- data[data != 'Methanobrevibacter']
  data <- unique(data)
  return(data)}

super_genera <- getGenera(super_abund$species)
success_genera <- getGenera(success_abund$species)
success_super_genera <- getGenera(success_super)



# Generate figure
library(VennDiagram)
# group1, group2, shared
draw.pairwise.venn(length(super_genera), length(success_genera), length(success_super_genera), 
                   category=c('Super\nDonor', 'Successful\nOutcome'), lty=rep('blank', 2), 
                   fill=c('darkorange', 'dodgerblue4'), alpha=rep(0.5, 2), cat.pos=c(0, 0), 
                   cat.dist=rep(0.025, 2), scaled=TRUE)

success_super
#rm(super_abund, normal_abund, success_abund, failure_abund)


# Pathway abundances
path <- '~/Desktop/Jenior_Consortia_2022/data/bin_gene_abundances/'
medabund_01044A <- read.delim(paste0(path,'01044A.gene_abundances/median_abundance.tsv'), sep='\t', header=TRUE)
medabund_05098A <- read.delim(paste0(path,'05098A.gene_abundances/median_abundance.tsv'), sep='\t', header=TRUE)
medabund_05080M <- read.delim(paste0(path,'05080M.gene_abundances/median_abundance.tsv'), sep='\t', header=TRUE)
medabund_28045A <- read.delim(paste0(path,'28045A.gene_abundances/median_abundance.tsv'), sep='\t', header=TRUE)
medabund_28043C <- read.delim(paste0(path,'28043C.gene_abundances/median_abundance.tsv'), sep='\t', header=TRUE)
medabund_28045D <- read.delim(paste0(path,'28045D.gene_abundances/median_abundance.tsv'), sep='\t', header=TRUE)
medabund_228047D <- read.delim(paste0(path,'28047D.gene_abundances/median_abundance.tsv'), sep='\t', header=TRUE)
rm(path)

# Merge by outcome
success_medabund <- merge(medabund_01044A, medabund_05098A, by='pathway')
success_medabund <- merge(success_medabund, medabund_05080M, by='pathway')
success_medabund <- merge(success_medabund, medabund_28043C, by='pathway')
success_medabund <- merge(success_medabund, medabund_28045D, by='pathway')
colnames(success_medabund) <- c('pathway','01044A','05098A','05080M','28043C','28045D')
failure_medabund <- merge(medabund_28045A, medabund_228047D, by='pathway')
colnames(failure_medabund) <- c('pathway','28045A','228047D')
rm(medabund_01044A, medabund_05098A, medabund_05080M, medabund_28043C, medabund_28045D, 
   medabund_28045A, medabund_228047D)

# Calculate medians
success_medabund$medabund <- apply(success_medabund[,2:ncol(success_medabund)], 1, median)
success_medabund$`01044A` <- NULL
success_medabund$`05098A` <- NULL
success_medabund$`05080M` <- NULL
success_medabund$`28043C` <- NULL
success_medabund$`28045D` <- NULL
failure_medabund$medabund <- apply(failure_medabund[,2:ncol(failure_medabund)], 1, median)
failure_medabund$`28045A` <- NULL
failure_medabund$`228047D` <- NULL

# Combine datasets
medabunds <- merge(success_medabund, failure_medabund, by='pathway')
colnames(medabunds) <- c('pathway','success_medabund','failure_medabund')
medabunds$abs_diff <- abs(medabunds$success_medabund - medabunds$failure_medabund)
medabunds <- subset(medabunds, abs_diff > 1.0)
medabunds$success_medabund_log <- log10(medabunds$success_medabund)
medabunds$failure_medabund_log <- log10(medabunds$failure_medabund)
rm(success_medabund, failure_medabund)

# Format names and remove host pathways
medabunds$pathway <- gsub('_', ' ', medabunds$pathway)
medabunds <- subset(medabunds, pathway != 'Type I diabetes mellitus')
medabunds <- subset(medabunds, pathway != 'Prion disease')
medabunds <- subset(medabunds, pathway != 'Pathways in cancer')
medabunds <- subset(medabunds, pathway != 'Peroxisome')
medabunds <- subset(medabunds, pathway != 'NOD-like receptor signaling pathway')
medabunds <- subset(medabunds, pathway != 'Pancreatic secretion')
medabunds <- subset(medabunds, pathway != 'Parkinson disease')
medabunds <- subset(medabunds, pathway != 'Pathways in cancer')
medabunds <- subset(medabunds, pathway != 'Pathways of neurodegeneration - multiple diseases')
medabunds <- subset(medabunds, pathway != 'Necroptosis')
medabunds <- subset(medabunds, pathway != 'MAPK signaling pathway - fly')
medabunds <- subset(medabunds, pathway != 'MAPK signaling pathway - plant')
medabunds <- subset(medabunds, pathway != 'MAPK signaling pathway - yeast')
medabunds <- subset(medabunds, pathway != 'Longevity regulating pathway - worm')
medabunds <- subset(medabunds, pathway != 'Longevity regulating pathway - multiple species')
medabunds <- subset(medabunds, pathway != 'Longevity regulating pathway')
medabunds <- subset(medabunds, pathway != 'Insulin signaling pathway')
medabunds <- subset(medabunds, pathway != 'Insulin resistance')
medabunds <- subset(medabunds, pathway != 'Insect hormone biosynthesis')
medabunds <- subset(medabunds, pathway != 'IL-17 signaling pathway')
medabunds <- subset(medabunds, pathway != 'Huntington disease')
medabunds <- subset(medabunds, pathway != 'Human T-cell leukemia virus 1 infection')
medabunds <- subset(medabunds, pathway != 'Human papillomavirus infection')
medabunds <- subset(medabunds, pathway != 'Hepatocellular carcinoma')
medabunds <- subset(medabunds, pathway != 'Glycosphingolipid biosynthesis - globo and isoglobo series')
medabunds <- subset(medabunds, pathway != 'Glycosphingolipid biosynthesis - ganglio series')
medabunds <- subset(medabunds, pathway != 'Fluid shear stress and atherosclerosis')
medabunds <- subset(medabunds, pathway != 'Cushing syndrome')
medabunds <- subset(medabunds, pathway != 'Axon regeneration')
medabunds <- subset(medabunds, pathway != 'Autophagy - yeast')
medabunds <- subset(medabunds, pathway != 'Apoptosis - fly')
medabunds <- subset(medabunds, pathway != 'Apoptosis')
medabunds <- subset(medabunds, pathway != 'Antigen processing and presentation')
medabunds <- subset(medabunds, pathway != 'Amoebiasis')
medabunds <- subset(medabunds, pathway != 'Alzheimer disease')
medabunds <- subset(medabunds, pathway != 'Nitrotoluene degradation')
medabunds <- subset(medabunds, pathway != 'Glutamatergic synapse')
medabunds <- subset(medabunds, pathway != 'Central carbon metabolism in cancer')
medabunds <- subset(medabunds, pathway != 'GABAergic synapse')
medabunds <- subset(medabunds, pathway != 'Thermogenesis')
medabunds <- subset(medabunds, pathway != 'Choline metabolism in cancer')
medabunds <- subset(medabunds, pathway != 'Carbon fixation in photosynthetic organisms')
medabunds <- subset(medabunds, pathway != 'AMPK signaling pathway')
medabunds <- subset(medabunds, pathway != 'African trypanosomiasis')
medabunds <- subset(medabunds, pathway != 'Adipocytokine signaling pathway')
medabunds <- subset(medabunds, pathway != 'Drug metabolism - cytochrome P450')
medabunds <- subset(medabunds, pathway != 'Drug metabolism - other enzymes')
medabunds <- subset(medabunds, pathway != 'Ferroptosis')
medabunds <- subset(medabunds, pathway != 'HIF-1 signaling pathway')
medabunds <- subset(medabunds, pathway != 'Glucagon signaling pathway')
medabunds <- subset(medabunds, pathway != 'Lysosome')
medabunds <- subset(medabunds, pathway != 'Meiosis - yeast')
medabunds <- subset(medabunds, pathway != 'MicroRNAs in cancer')
medabunds <- subset(medabunds, pathway != 'mRNA surveillance pathway')
medabunds <- subset(medabunds, pathway != 'Neuroactive ligand-receptor interaction')
medabunds <- subset(medabunds, pathway != 'PPAR signaling pathway')
medabunds <- subset(medabunds, pathway != 'Primary immunodeficiency')
medabunds <- subset(medabunds, pathway != 'Prolactin signaling pathway')
medabunds <- subset(medabunds, pathway != 'Prostate cancer')
medabunds <- subset(medabunds, pathway != 'Protein processing in endoplasmic reticulum')
medabunds <- subset(medabunds, pathway != 'Proteoglycans in cancer')
medabunds <- subset(medabunds, pathway != 'Proximal tubule bicarbonate reclamation')
medabunds <- subset(medabunds, pathway != 'Renal cell carcinoma')
medabunds <- subset(medabunds, pathway != 'Renin-angiotensin system')
medabunds <- subset(medabunds, pathway != 'Ribosome biogenesis in eukaryotes')
medabunds <- subset(medabunds, pathway != 'Salivary secretion')
medabunds <- subset(medabunds, pathway != 'Primary bile acid biosynthesis')
medabunds <- subset(medabunds, pathway != 'Spinocerebellar ataxia')
medabunds <- subset(medabunds, pathway != 'Th17 cell differentiation')
medabunds <- subset(medabunds, pathway != 'Type II diabetes mellitus')
medabunds <- subset(medabunds, pathway != 'Viral carcinogenesis')
medabunds <- subset(medabunds, pathway != 'African trypanosomiasis')
medabunds$pathway <- gsub('Arabinogalactan biosynthesis - Mycobacterium', 'Arabinogalactan biosynthesis', medabunds$pathway)
biofilm <- c('Biofilm formation - Escherichia coli','Biofilm formation - Pseudomonas aeruginosa','Biofilm formation - Vibrio cholerae')
biofilmabund <- subset(medabunds, pathway %in% biofilm)
medabunds <- subset(medabunds, !pathway %in% biofilm)
biofilm <- c('Biofilm formation',sum(biofilmabund$success_medabund),sum(biofilmabund$failure_medabund),
             abs(sum(biofilmabund$success_medabund)-sum(biofilmabund$failure_medabund)),
             log10(sum(biofilmabund$success_medabund)),log10(sum(biofilmabund$failure_medabund)))
medabunds <- as.data.frame(rbind(medabunds, biofilm))
rm(biofilm,biofilmabund)

# Remove antibiotic biosynthesis pathawys for future analysis
antibiotics <- c('Acarbose and validamycin biosynthesis','Novobiocin biosynthesis','Biosynthesis of vancomycin group antibiotics')
antibiotics <- subset(medabunds, pathway %in% antibiotics)
medabunds <- subset(medabunds, !pathway %in% antibiotics$pathway)


# Compiled pathway abundances
path <- '~/Desktop/Jenior_Consortia_2022/data/bin_gene_abundances/'
compiled_01044A <- read.delim(paste0(path,'01044A.gene_abundances/compiled_abundances.tsv'), sep='\t', header=FALSE)
colnames(compiled_01044A) <- c('pathway', paste('01044A', seq(1,ncol(compiled_01044A)-1), sep='_'))
compiled_05098A <- read.delim(paste0(path,'05098A.gene_abundances/compiled_abundances.tsv'), sep='\t', header=FALSE)
colnames(compiled_05098A) <- c('pathway', paste('05098A', seq(1,ncol(compiled_05098A)-1), sep='_'))
compiled_05080M <- read.delim(paste0(path,'05080M.gene_abundances/compiled_abundances.tsv'), sep='\t', header=FALSE)
colnames(compiled_05080M) <- c('pathway', paste('05080M', seq(1,ncol(compiled_05080M)-1), sep='_'))
compiled_28045A <- read.delim(paste0(path,'28045A.gene_abundances/compiled_abundances.tsv'), sep='\t', header=FALSE)
colnames(compiled_28045A) <- c('pathway', paste('28045A', seq(1,ncol(compiled_28045A)-1), sep='_'))
compiled_28043C <- read.delim(paste0(path,'28043C.gene_abundances/compiled_abundances.tsv'), sep='\t', header=FALSE)
colnames(compiled_28043C) <- c('pathway', paste('28043C', seq(1,ncol(compiled_28043C)-1), sep='_'))
compiled_28045D <- read.delim(paste0(path,'28045D.gene_abundances/compiled_abundances.tsv'), sep='\t', header=FALSE)
colnames(compiled_28045D) <- c('pathway', paste('28045D', seq(1,ncol(compiled_28045D)-1), sep='_'))
compiled_228047D <- read.delim(paste0(path,'28047D.gene_abundances/compiled_abundances.tsv'), sep='\t', header=FALSE)
colnames(compiled_228047D) <- c('pathway', paste('228047D', seq(1,ncol(compiled_228047D)-1), sep='_'))
rm(path)

# Combine by outcome
compiled_success <- merge(compiled_01044A, compiled_05098A, by='pathway')
#compiled_success <- merge(compiled_success, compiled_05080M, by='pathway')
#compiled_success <- merge(compiled_success, compiled_28043C, by='pathway')
compiled_success <- merge(compiled_success, compiled_28045D, by='pathway')
compiled_failed <- merge(compiled_28045A, compiled_228047D, by='pathway')
rm(compiled_01044A, compiled_05098A, compiled_05080M, compiled_28043C, 
   compiled_28045D, compiled_28045A, compiled_228047D)

# Remove NAs
compiled_success[is.na(compiled_success)] <- 0
compiled_failed[is.na(compiled_failed)] <- 0

# Subsample MAGS
#subSample <- sample(c(1:ncol(compiled_failed)), round(ncol(compiled_failed) * 0.9))
#compiled_success <- compiled_success[,c(1,subSample)]
#compiled_failed <- compiled_failed[,c(1,subSample)]
#rm(subSample)

# Remove host associated pathways
compiled_success$pathway <- gsub('_', ' ', compiled_success$pathway)
compiled_success <- subset(compiled_success, pathway != 'Carbon metabolism') # Too general
compiled_success <- subset(compiled_success, pathway != 'Metabolic pathways') # Too general
compiled_success <- subset(compiled_success, pathway != 'Type I diabetes mellitus')
compiled_success <- subset(compiled_success, pathway != 'Prion disease')
compiled_success <- subset(compiled_success, pathway != 'Pathways in cancer')
compiled_success <- subset(compiled_success, pathway != 'Peroxisome')
compiled_success <- subset(compiled_success, pathway != 'NOD-like receptor signaling pathway')
compiled_success <- subset(compiled_success, pathway != 'Pancreatic secretion')
compiled_success <- subset(compiled_success, pathway != 'Parkinson disease')
compiled_success <- subset(compiled_success, pathway != 'Pathways in cancer')
compiled_success <- subset(compiled_success, pathway != 'Pathways of neurodegeneration - multiple diseases')
compiled_success <- subset(compiled_success, pathway != 'Necroptosis')
compiled_success <- subset(compiled_success, pathway != 'MAPK signaling pathway - fly')
compiled_success <- subset(compiled_success, pathway != 'MAPK signaling pathway - plant')
compiled_success <- subset(compiled_success, pathway != 'MAPK signaling pathway - yeast')
compiled_success <- subset(compiled_success, pathway != 'Longevity regulating pathway - worm')
compiled_success <- subset(compiled_success, pathway != 'Longevity regulating pathway - multiple species')
compiled_success <- subset(compiled_success, pathway != 'Longevity regulating pathway')
compiled_success <- subset(compiled_success, pathway != 'Insulin signaling pathway')
compiled_success <- subset(compiled_success, pathway != 'Insulin resistance')
compiled_success <- subset(compiled_success, pathway != 'Insect hormone biosynthesis')
compiled_success <- subset(compiled_success, pathway != 'IL-17 signaling pathway')
compiled_success <- subset(compiled_success, pathway != 'Huntington disease')
compiled_success <- subset(compiled_success, pathway != 'Human T-cell leukemia virus 1 infection')
compiled_success <- subset(compiled_success, pathway != 'Human papillomavirus infection')
compiled_success <- subset(compiled_success, pathway != 'Hepatocellular carcinoma')
compiled_success <- subset(compiled_success, pathway != 'Glycosphingolipid biosynthesis - globo and isoglobo series')
compiled_success <- subset(compiled_success, pathway != 'Glycosphingolipid biosynthesis - ganglio series')
compiled_success <- subset(compiled_success, pathway != 'Fluid shear stress and atherosclerosis')
compiled_success <- subset(compiled_success, pathway != 'Cushing syndrome')
compiled_success <- subset(compiled_success, pathway != 'Axon regeneration')
compiled_success <- subset(compiled_success, pathway != 'Autophagy - yeast')
compiled_success <- subset(compiled_success, pathway != 'Apoptosis - fly')
compiled_success <- subset(compiled_success, pathway != 'Apoptosis')
compiled_success <- subset(compiled_success, pathway != 'Antigen processing and presentation')
compiled_success <- subset(compiled_success, pathway != 'Amoebiasis')
compiled_success <- subset(compiled_success, pathway != 'Alzheimer disease')
compiled_success <- subset(compiled_success, pathway != 'Nitrotoluene degradation')
compiled_success <- subset(compiled_success, pathway != 'Glutamatergic synapse')
compiled_success <- subset(compiled_success, pathway != 'Central carbon metabolism in cancer')
compiled_success <- subset(compiled_success, pathway != 'GABAergic synapse')
compiled_success <- subset(compiled_success, pathway != 'Thermogenesis')
compiled_success <- subset(compiled_success, pathway != 'Choline metabolism in cancer')
compiled_success <- subset(compiled_success, pathway != 'Carbon fixation in photosynthetic organisms')
compiled_success <- subset(compiled_success, pathway != 'AMPK signaling pathway')
compiled_success <- subset(compiled_success, pathway != 'African trypanosomiasis')
compiled_success <- subset(compiled_success, pathway != 'Adipocytokine signaling pathway')
compiled_success <- subset(compiled_success, pathway != 'Drug metabolism - cytochrome P450')
compiled_success <- subset(compiled_success, pathway != 'Drug metabolism - other enzymes')
compiled_success <- subset(compiled_success, pathway != 'Ferroptosis')
compiled_success <- subset(compiled_success, pathway != 'HIF-1 signaling pathway')
compiled_success <- subset(compiled_success, pathway != 'Glucagon signaling pathway')
compiled_success <- subset(compiled_success, pathway != 'Lysosome')
compiled_success <- subset(compiled_success, pathway != 'Meiosis - yeast')
compiled_success <- subset(compiled_success, pathway != 'MicroRNAs in cancer')
compiled_success <- subset(compiled_success, pathway != 'mRNA surveillance pathway')
compiled_success <- subset(compiled_success, pathway != 'Neuroactive ligand-receptor interaction')
compiled_success <- subset(compiled_success, pathway != 'PPAR signaling pathway')
compiled_success <- subset(compiled_success, pathway != 'Primary immunodeficiency')
compiled_success <- subset(compiled_success, pathway != 'Prolactin signaling pathway')
compiled_success <- subset(compiled_success, pathway != 'Prostate cancer')
compiled_success <- subset(compiled_success, pathway != 'Protein processing in endoplasmic reticulum')
compiled_success <- subset(compiled_success, pathway != 'Proteoglycans in cancer')
compiled_success <- subset(compiled_success, pathway != 'Proximal tubule bicarbonate reclamation')
compiled_success <- subset(compiled_success, pathway != 'Renal cell carcinoma')
compiled_success <- subset(compiled_success, pathway != 'Renin-angiotensin system')
compiled_success <- subset(compiled_success, pathway != 'Ribosome biogenesis in eukaryotes')
compiled_success <- subset(compiled_success, pathway != 'Salivary secretion')
compiled_success <- subset(compiled_success, pathway != 'Primary bile acid biosynthesis')
compiled_success <- subset(compiled_success, pathway != 'Spinocerebellar ataxia')
compiled_success <- subset(compiled_success, pathway != 'Th17 cell differentiation')
compiled_success <- subset(compiled_success, pathway != 'Type II diabetes mellitus')
compiled_success <- subset(compiled_success, pathway != 'Viral carcinogenesis')
compiled_success <- subset(compiled_success, pathway != 'African trypanosomiasis')
compiled_success <- subset(compiled_success, pathway != 'Protein digestion and absorption')
compiled_success <- subset(compiled_success, pathway != 'FoxO signaling pathway')
compiled_success$pathway <- gsub('Arabinogalactan biosynthesis - Mycobacterium', 'Arabinogalactan biosynthesis', compiled_success$pathway)
biofilm <- c('Biofilm formation - Escherichia coli','Biofilm formation - Pseudomonas aeruginosa','Biofilm formation - Vibrio cholerae')
biofilmabund <- subset(compiled_success, pathway %in% biofilm)
compiled_success <- subset(compiled_success, !pathway %in% biofilm)
biofilm <- c('Biofilm formation',sum(biofilmabund$success_medabund),sum(biofilmabund$failure_medabund),
             abs(sum(biofilmabund$success_medabund)-sum(biofilmabund$failure_medabund)),
             log10(sum(biofilmabund$success_medabund)),log10(sum(biofilmabund$failure_medabund)))
compiled_success <- as.data.frame(rbind(compiled_success, biofilm))
rm(biofilm,biofilmabund)

compiled_failed$pathway <- gsub('_', ' ', compiled_failed$pathway)
compiled_failed <- subset(compiled_failed, pathway != 'Carbon metabolism') # Too general
compiled_failed <- subset(compiled_failed, pathway != 'Metabolic pathways') # Too general
compiled_failed <- subset(compiled_failed, pathway != 'Type I diabetes mellitus')
compiled_failed <- subset(compiled_failed, pathway != 'Prion disease')
compiled_failed <- subset(compiled_failed, pathway != 'Pathways in cancer')
compiled_failed <- subset(compiled_failed, pathway != 'Peroxisome')
compiled_failed <- subset(compiled_failed, pathway != 'NOD-like receptor signaling pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Pancreatic secretion')
compiled_failed <- subset(compiled_failed, pathway != 'Parkinson disease')
compiled_failed <- subset(compiled_failed, pathway != 'Pathways in cancer')
compiled_failed <- subset(compiled_failed, pathway != 'Pathways of neurodegeneration - multiple diseases')
compiled_failed <- subset(compiled_failed, pathway != 'Necroptosis')
compiled_failed <- subset(compiled_failed, pathway != 'MAPK signaling pathway - fly')
compiled_failed <- subset(compiled_failed, pathway != 'MAPK signaling pathway - plant')
compiled_failed <- subset(compiled_failed, pathway != 'MAPK signaling pathway - yeast')
compiled_failed <- subset(compiled_failed, pathway != 'Longevity regulating pathway - worm')
compiled_failed <- subset(compiled_failed, pathway != 'Longevity regulating pathway - multiple species')
compiled_failed <- subset(compiled_failed, pathway != 'Longevity regulating pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Insulin signaling pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Insulin resistance')
compiled_failed <- subset(compiled_failed, pathway != 'Insect hormone biosynthesis')
compiled_failed <- subset(compiled_failed, pathway != 'IL-17 signaling pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Huntington disease')
compiled_failed <- subset(compiled_failed, pathway != 'Human T-cell leukemia virus 1 infection')
compiled_failed <- subset(compiled_failed, pathway != 'Human papillomavirus infection')
compiled_failed <- subset(compiled_failed, pathway != 'Hepatocellular carcinoma')
compiled_failed <- subset(compiled_failed, pathway != 'Glycosphingolipid biosynthesis - globo and isoglobo series')
compiled_failed <- subset(compiled_failed, pathway != 'Glycosphingolipid biosynthesis - ganglio series')
compiled_failed <- subset(compiled_failed, pathway != 'Fluid shear stress and atherosclerosis')
compiled_failed <- subset(compiled_failed, pathway != 'Cushing syndrome')
compiled_failed <- subset(compiled_failed, pathway != 'Axon regeneration')
compiled_failed <- subset(compiled_failed, pathway != 'Autophagy - yeast')
compiled_failed <- subset(compiled_failed, pathway != 'Apoptosis - fly')
compiled_failed <- subset(compiled_failed, pathway != 'Apoptosis')
compiled_failed <- subset(compiled_failed, pathway != 'Antigen processing and presentation')
compiled_failed <- subset(compiled_failed, pathway != 'Amoebiasis')
compiled_failed <- subset(compiled_failed, pathway != 'Alzheimer disease')
compiled_failed <- subset(compiled_failed, pathway != 'Nitrotoluene degradation')
compiled_failed <- subset(compiled_failed, pathway != 'Glutamatergic synapse')
compiled_failed <- subset(compiled_failed, pathway != 'Central carbon metabolism in cancer')
compiled_failed <- subset(compiled_failed, pathway != 'GABAergic synapse')
compiled_failed <- subset(compiled_failed, pathway != 'Thermogenesis')
compiled_failed <- subset(compiled_failed, pathway != 'Choline metabolism in cancer')
compiled_failed <- subset(compiled_failed, pathway != 'Carbon fixation in photosynthetic organisms')
compiled_failed <- subset(compiled_failed, pathway != 'AMPK signaling pathway')
compiled_failed <- subset(compiled_failed, pathway != 'African trypanosomiasis')
compiled_failed <- subset(compiled_failed, pathway != 'Adipocytokine signaling pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Drug metabolism - cytochrome P450')
compiled_failed <- subset(compiled_failed, pathway != 'Drug metabolism - other enzymes')
compiled_failed <- subset(compiled_failed, pathway != 'Ferroptosis')
compiled_failed <- subset(compiled_failed, pathway != 'HIF-1 signaling pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Glucagon signaling pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Lysosome')
compiled_failed <- subset(compiled_failed, pathway != 'Meiosis - yeast')
compiled_failed <- subset(compiled_failed, pathway != 'MicroRNAs in cancer')
compiled_failed <- subset(compiled_failed, pathway != 'mRNA surveillance pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Neuroactive ligand-receptor interaction')
compiled_failed <- subset(compiled_failed, pathway != 'PPAR signaling pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Primary immunodeficiency')
compiled_failed <- subset(compiled_failed, pathway != 'Prolactin signaling pathway')
compiled_failed <- subset(compiled_failed, pathway != 'Prostate cancer')
compiled_failed <- subset(compiled_failed, pathway != 'Protein processing in endoplasmic reticulum')
compiled_failed <- subset(compiled_failed, pathway != 'Proteoglycans in cancer')
compiled_failed <- subset(compiled_failed, pathway != 'Proximal tubule bicarbonate reclamation')
compiled_failed <- subset(compiled_failed, pathway != 'Renal cell carcinoma')
compiled_failed <- subset(compiled_failed, pathway != 'Renin-angiotensin system')
compiled_failed <- subset(compiled_failed, pathway != 'Ribosome biogenesis in eukaryotes')
compiled_failed <- subset(compiled_failed, pathway != 'Salivary secretion')
compiled_failed <- subset(compiled_failed, pathway != 'Primary bile acid biosynthesis')
compiled_failed <- subset(compiled_failed, pathway != 'Spinocerebellar ataxia')
compiled_failed <- subset(compiled_failed, pathway != 'Th17 cell differentiation')
compiled_failed <- subset(compiled_failed, pathway != 'Type II diabetes mellitus')
compiled_failed <- subset(compiled_failed, pathway != 'Viral carcinogenesis')
compiled_failed <- subset(compiled_failed, pathway != 'African trypanosomiasis')
compiled_failed <- subset(compiled_failed, pathway != 'Protein digestion and absorption')
compiled_failed <- subset(compiled_failed, pathway != 'FoxO signaling pathway')
compiled_failed$pathway <- gsub('Arabinogalactan biosynthesis - Mycobacterium', 'Arabinogalactan biosynthesis', compiled_failed$pathway)
biofilm <- c('Biofilm formation - Escherichia coli','Biofilm formation - Pseudomonas aeruginosa','Biofilm formation - Vibrio cholerae')
biofilmabund <- subset(compiled_failed, pathway %in% biofilm)
compiled_failed <- subset(compiled_failed, !pathway %in% biofilm)
biofilm <- c('Biofilm formation',sum(biofilmabund$failed_medabund),sum(biofilmabund$failure_medabund),
             abs(sum(biofilmabund$failed_medabund)-sum(biofilmabund$failure_medabund)),
             log10(sum(biofilmabund$failed_medabund)),log10(sum(biofilmabund$failure_medabund)))
compiled_failed <- as.data.frame(rbind(compiled_failed, biofilm))
rm(biofilm,biofilmabund)

# Fix pathway names
compiled_success$pathway <- make.names(compiled_success$pathway)
rownames(compiled_success) <- compiled_success$pathway
compiled_success$pathway <- NULL
compiled_failed$pathway <- make.names(compiled_failed$pathway)
rownames(compiled_failed) <- compiled_failed$pathway
compiled_failed$pathway <- NULL

# Add group labels
compiled_success <- as.data.frame(t(compiled_success))
compiled_success$outcome <- rep('success', nrow(compiled_success))
compiled_failed <- as.data.frame(t(compiled_failed))
compiled_failed$outcome <- rep('failed', nrow(compiled_failed))

# Merge for feature selection
compiled_success <- as.data.frame(t(compiled_success))
compiled_failed <- as.data.frame(t(compiled_failed))
compiled_abund <- merge(compiled_success, compiled_failed, by='row.names')
rownames(compiled_abund) <- compiled_abund$Row.names
compiled_abund$Row.names <- NULL
compiled_abund <- as.data.frame(t(compiled_abund))

# Run random forest
library(randomForest)
set.seed(9861)
condition <- as.factor(compiled_abund$outcome)
compiled_abund$outcome <- NULL
compiled_abund <- droplevels(compiled_abund)
for (x in 1:ncol(compiled_abund)) {compiled_abund[,x] <- as.numeric(as.character(compiled_abund[,x]))}
compiled_abund[compiled_abund == -Inf] <- 0
mTries <- round(sqrt(ncol(compiled_abund)))
nTrees <- as.numeric(length(as.vector(levels(condition))) * ncol(compiled_abund))
rf_obj <- randomForest(condition ~ ., data=compiled_abund, importance=TRUE, err.rate=TRUE, na.action=na.roughfix, 
                       ntree=nTrees, mtry=mTries)
print(rf_obj)
# Number of trees: 368
# No. of variables tried at each split: 14
# OOB estimate of  error rate: 25.88%
# Confusion matrix:
#   failed success class.error
# failed      93      21   0.1842105
# success     38      76   0.3333333
rf_obj <- importance(rf_obj, type=1, scale=TRUE)
rf_mda <- subset(rf_obj, rf_obj > (abs(min(rf_obj)))) # significance
rf_mda <- as.data.frame(rf_mda)
rf_mda$feature <- rownames(rf_mda)
rf_mda <- rf_mda[order(rf_mda$MeanDecreaseAccuracy),]
rf_mda <- subset(rf_mda, MeanDecreaseAccuracy >= 3)
rf_mda$feature <- gsub('\\.', ' ', rf_mda$feature)
rf_mda$feature <- gsub('Glycolysis   Gluconeogenesis', 'Glycolysis/Gluconeogenesis', rf_mda$feature)
rm(x, nTrees, mTries, rf_obj, compiled_abund, condition)

rf_mda <- read.delim('~/Desktop/Jenior_Consortia_2022/data/mag_abund_rf.tsv', sep='\t', header=TRUE)
success_mda <- subset(rf_mda, color == 'gray60')
failure_mda <- subset(rf_mda, color == 'white')


# MDA plotting
pdf(file='~/Desktop/Jenior_Consortia_2022/results/success_mag_abund_rf.pdf', width=4, height=3)
par(mar=c(3, 0.5, 0.5, 0.5), mgp=c(1.4, 0.7, 0), xpd=FALSE, lwd=2)
dotchart(success_mda$MeanDecreaseAccuracy, bg=success_mda$color, xlim=c(0,6),  
           pch=21, lwd=1.7, pt.cex=1.7, cex=0.8)
text(x=-0.2, y=seq(1.4,nrow(rf_mda)+0.4,1), labels=success_mda$feature, cex=0.95, pos=4)
mtext('Mean Decrease Accuracy', side=1, padj=2.5)
#legend('bottomright', legend=c('Greater in Success','Greater in Failure'), cex=0.8, pt.cex=1.4, 
#       pch=21, pt.bg=c(success_col,failure_col), bg='white')
box()
dev.off()

pdf(file='~/Desktop/Jenior_Consortia_2022/results/failure_mag_abund_rf.pdf', width=4, height=2)
par(mar=c(3, 0.5, 0.5, 0.5), mgp=c(1.4, 0.7, 0), xpd=FALSE, lwd=2)
dotchart(failure_mda$MeanDecreaseAccuracy, bg=failure_mda$color, xlim=c(0,6),  
         pch=21, lwd=1.7, pt.cex=1.7, cex=0.8)
text(x=-0.2, y=seq(1.4,nrow(rf_mda)+0.4,1), labels=failure_mda$feature, cex=0.95, pos=4)
mtext('Mean Decrease Accuracy', side=1, padj=2.5)
#legend('bottomright', legend=c('Greater in Success','Greater in Failure'), cex=0.8, pt.cex=1.4, 
#       pch=21, pt.bg=c(success_col,failure_col), bg='white')
box()
dev.off()

write.table(rf_mda, file='~/Desktop/Jenior_Consortia_2022/data/mag_abund_rf.tsv', quote=FALSE, sep='\t',
            row.names=FALSE, col.names=TRUE)

# Check medians
compiled_success <- as.data.frame(t(compiled_success))
compiled_failed <- as.data.frame(t(compiled_failed))

greaterMedian <- function(pathway) {
  success_med <- median(as.numeric(compiled_success[,pathway]))
  failed_med <- median(as.numeric(compiled_failed[,pathway]))
  if (success_med >= failed_med) {print('Success')} else {print('Failure')}}

greaterMedian('Homologous.recombination')
greaterMedian('Biosynthesis.of.secondary.metabolites')
greaterMedian('Mismatch.repair')
greaterMedian('Phenylpropanoid.biosynthesis')
greaterMedian('ABC.transporters')
greaterMedian('Starch.and.sucrose.metabolism')
greaterMedian('Valine..leucine.and.isoleucine.biosynthesis')
greaterMedian('Amino.sugar.and.nucleotide.sugar.metabolism')
greaterMedian('DNA.replication')
greaterMedian('Nicotinate.and.nicotinamide.metabolism')
greaterMedian('Glycine..serine.and.threonine.metabolism')
greaterMedian('Cysteine.and.methionine.metabolism')
greaterMedian('Aminoacyl.tRNA.biosynthesis')
greaterMedian('Quorum.sensing')
greaterMedian('Staphylococcus.aureus.infection')
greaterMedian('Phenylalanine..tyrosine.and.tryptophan.biosynthesis')
greaterMedian('Galactose.metabolism')
greaterMedian('Purine.metabolism')
greaterMedian('Biosynthesis.of.amino.acids')

success_col <- 'chocolate2'
failure_col <- 'blue3'
mda_cols <- c(failure_col,success_col,success_col,success_col,success_col,success_col,success_col,
              success_col,failure_col,success_col,failure_col,failure_col,success_col,success_col,
              success_col,failure_col,success_col,failure_col,success_col)




compiled_abund$Biofilm.formation <-NULL
compiled_abund$outcome <- NULL
compiled_abund <- as.data.frame(t(compiled_abund))
for (x in 1:ncol(compiled_abund)) {
  compiled_abund[,x] <- as.numeric(compiled_abund[,x])
  compiled_abund[,x] <- (compiled_abund[,x] / sum(compiled_abund[,x])) * 100}
compiled_abund <- as.data.frame(t(compiled_abund))

focused_abund <- compiled_abund[,c('Biosynthesis.of.amino.acids','Galactose.metabolism','Amino.sugar.and.nucleotide.sugar.metabolism','Starch.and.sucrose.metabolism')]
focused_abund$combined_carbohydrate <- as.numeric(focused_abund$Galactose.metabolism) + as.numeric(focused_abund$Amino.sugar.and.nucleotide.sugar.metabolism) + as.numeric(focused_abund$Starch.and.sucrose.metabolism)
focused_abund$Galactose.metabolism <- NULL
focused_abund$Amino.sugar.and.nucleotide.sugar.metabolism <- NULL
focused_abund$Starch.and.sucrose.metabolism <- NULL
focused_abund$Biosynthesis.of.amino.acids <- as.numeric(focused_abund$Biosynthesis.of.amino.acids)

aa_cutoff <- as.numeric(quantile(focused_abund$Biosynthesis.of.amino.acids, 0.9))
special_interest1 <- subset(focused_abund, Biosynthesis.of.amino.acids > aa_cutoff)
carb_cutoff <- as.numeric(quantile(focused_abund$combined_carbohydrate, 0.9))
special_interest2 <- subset(focused_abund, combined_carbohydrate > carb_cutoff)
special_interest <- rbind(special_interest1, special_interest2)

special_interest$mag <- rownames(special_interest)
write.table(special_interest, file='~/Desktop/special_interest.tsv', quote=FALSE, sep='\t',
            row.names=FALSE, col.names=TRUE)

special_interest <- read.delim('~/Desktop/special_interest.tsv', sep='\t', header=TRUE)

xmin <- min(focused_abund$Biosynthesis.of.amino.acids, na.rm=TRUE) * 0.9
xmax <- max(focused_abund$Biosynthesis.of.amino.acids, na.rm=TRUE) * 1.1
ymin <- min(focused_abund$combined_carbohydrate, na.rm=TRUE) * 0.9
ymax <- max(focused_abund$combined_carbohydrate, na.rm=TRUE) * 1.1

library(scales)
pdf(file='~/Desktop/Jenior_Consortia_2022/results/mag_sctterplot.pdf', width=4, height=4)
par(mar=c(3,3,1,1), mgp=c(2, 0.6, 0), new=FALSE, xpd=FALSE, lwd=2, las=1)
plot(x=focused_abund$Biosynthesis.of.amino.acids, y=focused_abund$combined_carbohydrate, 
     xlim=c(xmin,xmax), ylim=c(ymin,ymax), pch=16, col=alpha('gray40',0.5), cex=1.2, cex.axis=0.8,
     xlab='Amino acid biosynthesis genes (%)', ylab='Combined carbohydrate metabolism genes (%)')
abline(v=aa_cutoff, lty=5, col='gray60')
abline(h=carb_cutoff, lty=5, col='gray60')
points(x=special_interest$Biosynthesis.of.amino.acids, y=special_interest$combined_carbohydrate, 
       pch=21, bg=special_interest$color, cex=1.5)
legend('topright', legend=c('Donor A','Donor B','Donor C'), 
       pt.bg=c('chocolate1','cornflowerblue','chartreuse3'), pch=21, pt.cex=2, cex=1.1)
box()
dev.off()




#---------------------------------#

# Correlation analysis
# Plotting function
corr_plot <- function(data, x_label='X', y_label='Y', x_cutoff=0.25, y_cutoff=0.75) {
  test <- cor.test(x=as.numeric(data[,2]), y=as.numeric(data[,3]), method='spearman', exact=FALSE)
  pval <- round(test$p.value, 3)
  rho <- round(test$estimate, 3)
  reg <- lm(as.numeric(data[,3]) ~ as.numeric(data[,2])) 
  xmin <- min(as.numeric(data[,2]), na.rm=TRUE) * 0.9
  xmax <- max(as.numeric(data[,2]), na.rm=TRUE) * 1.1
  x_cutoff <- as.numeric(quantile(as.numeric(data[,2]), x_cutoff))
  ymin <- min(as.numeric(data[,3]), na.rm=TRUE) * 0.9
  ymax <- max(as.numeric(data[,3]), na.rm=TRUE) * 1.1
  y_cutoff <- as.numeric(quantile(as.numeric(data[,3]), y_cutoff))
  
  sub_data <- subset(data, data[,2] >= x_cutoff)
  sub_data <- subset(sub_data, sub_data[,3] <= y_cutoff)
  sub_data <- subset(sub_data, sub_data[,2] > sub_data[,3])
  
  print(nrow(sub_data))
  sub_genera <- getGenera(sub_data$species)
  print(length(sub_genera))

  par(mar=c(2.5,2.5,0.5,0.5), mgp=c(1.5, 0.6, 0), new=FALSE, xpd=FALSE, lwd=2, las=1)
  plot(x=as.numeric(data[,2]), y=as.numeric(data[,3]), xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       pch=21, bg='dodgerblue3', cex=1.5, cex.axis=0.7, xlab=x_label, ylab=y_label)
  abline(v=x_cutoff, col='gray70', lty=5)
  abline(h=y_cutoff, col='gray70', lty=5)
  abline(reg, col='firebrick4', lwd=3)
  points(x=sub_data[,2], y=sub_data[,3], pch=21, bg='green3', cex=1.5)
  box()
  
  if (pval <= 0.05) {pval_str <- as.expression(bquote(paste(italic('p'),' = ', .(as.character(pval)))))} 
  if (pval == 0) {pval_str <- as.expression(bquote(paste(italic('p'),' < 0.001')))
  } else {pval_str <- as.expression(bquote(paste(italic('p'),' = ', .(as.character(pval)))))}
  rho_str <- paste('R =', as.character(rho))
  legend('topright', legend=c(pval_str, rho_str), cex=0.7, pt.cex=0, bty='n')
  
  text_ycoord <- seq(ymin, as.numeric(quantile(as.numeric(data[,3]), 0.6)), ymax*0.022)[1:length(sub_genera)]
  text(x=xmax*0.65, y=text_ycoord, labels=sub_genera, cex=0.65, font=3, pos=4)}

# Generate figures

# MAG abundances
png(filename='~/Desktop/Jenior_Consortia_2022/results/mag_abund_corr.png', units='in', width=6, height=3, res=300)
layout(matrix(c(1,2), nrow=2, ncol=2, byrow=TRUE))
corr_plot(donor_abund, x_label='Super Donor Abundance (Log10)', y_label='Normal Donor Abundance (Log10)')
corr_plot(outcome_abund, x_label='Success Abundance (Log10)', y_label='Failure Abundance (Log10)')
dev.off()

png(filename='~/Desktop/Jenior_Consortia_2022/results/mag_abund_venn.png', units='in', width=6, height=3, res=300)
draw.pairwise.venn(length(super_genera), length(success_genera), length(success_super_genera), 
                   category=c('', ''), lty=rep('blank', 2), label.col='white', cex=2, 
                   fill=c('darkorange', 'dodgerblue4'), alpha=rep(0.5, 2), cat.pos=c(0, 0), 
                   cat.dist=rep(0.025, 2), scaled=TRUE)
print(success_super)
dev.off()

# Pathway abundances
test <- cor.test(x=as.numeric(medabunds$success_medabund_log), y=as.numeric(medabunds$failure_medabund_log), method='spearman', exact=FALSE)
pval <- round(test$p.value, 3)
rho <- round(test$estimate, 3)
reg <- lm(as.numeric(medabunds$failure_medabund_log) ~ as.numeric(medabunds$success_medabund_log)) 
xmin <- min(as.numeric(medabunds$success_medabund_log), na.rm=TRUE) * 0.9
xmax <- max(as.numeric(medabunds$success_medabund_log), na.rm=TRUE) * 1.1
x_cutoff <- as.numeric(quantile(as.numeric(medabunds$success_medabund_log), 0.4))
ymin <- min(as.numeric(medabunds$failure_medabund_log), na.rm=TRUE) * 0.9
ymax <- max(as.numeric(medabunds$failure_medabund_log), na.rm=TRUE) * 1.1
y_cutoff <- as.numeric(quantile(as.numeric(medabunds$failure_medabund_log), 0.6))
sub_medabunds <- subset(medabunds, medabunds$success_medabund_log >= x_cutoff)
sub_medabunds <- subset(sub_medabunds, sub_medabunds$failure_medabund_log <= y_cutoff)
sub_medabunds <- subset(sub_medabunds, sub_medabunds$success_medabund_log > sub_medabunds$failure_medabund_log)
sub_medabunds <- sub_medabunds[order(as.numeric(sub_medabunds$abs_diff)),] 

png(filename='~/Desktop/Jenior_Consortia_2022/results/pathway_abund_corr.png', units='in', width=5, height=4, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
par(mar=c(2.5,2.5,0.5,0.5), mgp=c(1.5, 0.6, 0), new=FALSE, xpd=FALSE, lwd=2, las=1)
plot(x=as.numeric(medabunds$success_medabund_log), y=as.numeric(medabunds$failure_medabund_log), xlim=c(1.5,3), ylim=c(1,3), cex.lab=0.8,
     pch=21, bg='dodgerblue3', cex=1.5, cex.axis=0.7, xlab='Success Med. Abund. (Log10)', ylab='Failure Med. Abund. (Log10)')
abline(v=x_cutoff, col='gray70', lty=5)
abline(h=y_cutoff, col='gray70', lty=5)
abline(reg, col='firebrick4', lwd=3)
points(x=sub_medabunds$success_medabund_log, y=sub_medabunds$failure_medabund_log, pch=21, bg='firebrick2', cex=1.5)
box()

if (pval <= 0.05) {pval_str <- as.expression(bquote(paste(italic('p'),' = ', .(as.character(pval)))))} 
if (pval == 0) {pval_str <- as.expression(bquote(paste(italic('p'),' < 0.001')))
} else {pval_str <- as.expression(bquote(paste(italic('p'),' = ', .(as.character(pval)))))}
rho_str <- paste('R =', as.character(rho))
legend('bottomright', legend=c(pval_str, rho_str), cex=0.9, pt.cex=0, bty='n')

par(mar=c(1,1,1,1), xpd=TRUE)
plot(0, type='n', ylim=c(0,27), xlim=c(0,20), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
text(x=c(-3.8,3), y=27, labels=c('Abs Diff','KEGG Pathway'), cex=0.6, pos=4, font=2)
text(x=0.4, y=c(1:26), labels=sub_medabunds$pathway, cex=0.55, pos=4)
text(x=-2.5, y=c(1:26), labels=round(as.numeric(sub_medabunds$abs_diff)), cex=0.6, pos=4, font=2)

dev.off()


#rm(corr_plot)

#------------------------------------------------------------------------------------------------------------#

# Carbohydrate metabolism
# Amino acid biosynthesis

paths <- read.delim('~/Desktop/Jenior_Consortia_2022/data/mag_pathway_counts.tsv', sep='\t', header=TRUE, row.names=1)
paths$total <- paths$total * 0.3
paths$aminoacid <- (paths$aminoacid / paths$total) * 100
paths$carbohydrate <- (paths$carbohydrate / paths$total) * 100
paths$total <- NULL
paths$aminoacid_perc <- NULL
paths$carbohydrate_perc <- NULL
paths$aminoacid_log <- NULL
paths$carbohydrate_log <- NULL

paths <- t(paths)
paths <- paths[,order(ncol(paths):1)]
pdf(file='~/Desktop/Jenior_Consortia_2022/results/mag_enriched_paths.pdf', width=2, height=8)
par(mar=c(3,1,0.5,1), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2)
barplot(height=paths, beside=TRUE, horiz=TRUE, xlim=c(0,30), xaxt='n', cex.lab=0.9,
        xlab='Percent metabolic genes', col=c('white','gray40'))
box()
axis(side=1, at=c(0,15,30), cex.axis=0.8, lwd=2)



dev.off()

#------------------------------------------------------------------------------------------------------------#

# Data
consistency <- c(51.2, 48.1, 48.1, 48.3, 45.9, 50.3, 45, 42.8, 42.8, 47.9, 46.9, 49.6, 53.3, 48.4, 47.6, 48, 43.6, 47.7, 44.9, 43.8, 48, 47.4)
balance <- c(88.8, 87.9, 86, 88.1, 87.4, 88, 79.6, 86.2, 86, 88.4, 89.3, 87.9, 87, 88, 87.6, 87.3, 85.8, 85.3, 86.2, 85.8, 86.7, 86.9)
memote <- c(75, 74, 75, 73, 74, 72, 75, 74, 74, 72, 74, 72, 75, 73, 73, 74, 75, 75, 73, 74, 75, 75)


pdf(file='~/Desktop/Jenior_Consortia_2022/results/supplemented.pdf', width=5, height=3)
par(mar=c(6,3,0.5,0.5), las=1, mgp=c(1.7,0.7,0), xpd=FALSE, lwd=2.5)
plot(0, type='n', ylim=c(0,3), xlim=c(0,100), 
     xlab='Average score (%)', ylab='', xaxt='n', yaxt='n', cex.lab=1.1)
boxplot(consistency, cex=0, lwd=2, at=0.5, col='antiquewhite1', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', horizontal=TRUE, add=TRUE)
boxplot(balance, cex=0, lwd=2, at=1.5, col='antiquewhite1', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', horizontal=TRUE, add=TRUE)
boxplot(memote, cex=0, lwd=2, at=2.5, col='antiquewhite1', ylab='', staplewex=0.6, 
        boxwex=1.5, lty=1, medlwd=2, xaxt='n', yaxt='n', horizontal=TRUE, add=TRUE)


axis(side=2, at=seq(0,80,20), cex.axis=0.8, lwd=2.5)
segments(x0=2.5, x1=3.5, y0=50)
text(x=3, y=54, '*', cex=1.5, font=2)
par(xpd=TRUE)
text(x=c(0.5,1.5,2.5,3.5), y=-30, srt=45, cex=0.8,
     c('R20291 depleted','R20291 depleted +\nFresh media','R20291 depleted +\nCooperative depleted','R20291 depleted +\nCompetative depleted'))
dev.off()



# Calculate summary stats
consistency <- quantile(consistency)
balance <- quantile(balance)
memote <- quantile(memote)

# Generate figure
#pdf(file='~/Desktop/Jenior_Consortia_2022/results/figure_3B.pdf', width=1.5, height=4)
png(filename='~/Desktop/Jenior_Consortia_2022/results/figure_3B.png', units='in', width=1.5, height=4, res=300)

par(mar=c(9,0.5,1,3), mgp=c(1.3,0.4,0), lwd=2, xaxs='i', yaxs='i', xpd=FALSE)
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,3.25), ylim=c(0,100), ylab='Average score (%)', xlab='')
axis(4, at=c(0,20,40,60,80,100), cex.axis=0.7, lwd=2)

rect(xleft=0.25, ybottom=0, xright=1, ytop=consistency[3], col='gray75')
segments(x0=0.35, x1=0.9, y0=consistency[2])
segments(x0=0.35, x1=0.9, y0=consistency[4])
segments(x0=0.625, y0=consistency[2], y1=consistency[4])

rect(xleft=1.25, ybottom=0, xright=2, ytop=balance[3], col='gray75')
segments(x0=1.35, x1=1.9, y0=balance[2])
segments(x0=1.35, x1=1.9, y0=balance[4])
segments(x0=1.625, y0=balance[2], y1=balance[4])

rect(xleft=2.25, ybottom=0, xright=3, ytop=memote[3], col='gray75')
segments(x0=2.35, x1=2.9, y0=memote[2])
segments(x0=2.35, x1=2.9, y0=memote[4])
segments(x0=2.625, y0=memote[2], y1=memote[4])

par(xpd=TRUE)
text(x=4.75, y=80, 'Average score (%)', srt=90, adj=1, cex=0.8)
text(x=c(0.625,1.625,2.625), y=-6, c('Stoichiometric consistency','Mass balance','MEMOTE'), 
     srt=90, adj=1, cex=0.8)

dev.off()

