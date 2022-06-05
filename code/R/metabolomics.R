
# Read in and format data
metabolome <- read.delim('~/Desktop/Jenior_Consortia_2022/data/combined_metabolome.tsv', sep='\t', header=TRUE)
exclude <- c('ZDXPYRJPNDTMRX-VKHMYHEASA-N','CQIUKKVOEOPUDV-IYSWYEEDSA-N',
             'DAXYUDFNWXHGBE-VKCGGMIFSA-N','LELBFTMXCIIKKX-UHFFFAOYNA-N',
             'SPLUHYDIQKURQW-UHFFFAOYSA-N','HJBWJAPEBGSQPR-GQCTYLIASA-N')
metabolome <- subset(metabolome, !INCHIKEY %in% exclude)
metabolome$MetaboliteName <- make.names(metabolome$MetaboliteName)
#inchikey <- metabolome[,c(1:2)]
metabolome$INCHIKEY <- NULL
rownames(metabolome) <- metabolome$MetaboliteName
metabolome$MetaboliteName <- NULL
metabolome$blank_1 <- NULL
metabolome$blank_2 <- NULL
metabolome$blank_3 <- NULL
metabolome$blank_4 <- NULL
metabolome <- as.data.frame(t(metabolome))
metadata <- read.delim('~/Desktop/Jenior_Consortia_2022/data/metabolome_metadata.tsv', sep='\t', header=TRUE)
metabolome <- as.data.frame(merge(metadata, metabolome, by.x='Label', by.y='row.names'))
rownames(metabolome) <- metabolome$Label
metabolome$Label <- NULL
rm(exclude)

# Focus on validation set
keep <- c('C_difficile_R20291','B_producta_DSM2950','B_longum_ATCC55813','Fresh_media')
metabolome <- subset(metabolome, Depleting_species %in% keep)
fresh <- subset(metabolome, Depleting_species == 'Fresh_media')
fresh$R20291_secondary_growth <- NULL
fresh$Depleting_species <- NULL

#-------------------------------------------------------------#

# Define functions

# Create subset of decreased growth only
decreasedSubset <- function(data, species) {
  keep <- subset(data, R20291_secondary_growth == 'no')
  keep <- subset(keep, Depleting_species %in% c(species,'Fresh_media'))
  keep$R20291_secondary_growth <- NULL
  keep <- droplevels(keep)
  return(keep)}

# Test for significant differences - decreased
sigDecreased <- function(data) {
  group1 <- subset(data, Depleting_species == 'Fresh_media')
  group1$Depleting_species <- NULL
  group2 <- subset(data, Depleting_species != 'Fresh_media')
  group2$Depleting_species <- NULL
  keep <- c()
  for (x in colnames(group1)) {
    test1 <- as.vector(group1[,x])
    test2 <- as.vector(group2[,x])
    pval <- as.numeric(wilcox.test(test1, test2, exact=FALSE)$p.value)
    if (is.na(pval)) {pval = 1}
    if (pval <= 0.05) {if (median(test1) > median(test2)) {keep <- c(keep, x)}}
  }
  group2 <- group2[, keep]
  return(group2)}

findDecrease <- function(data, species='default') {
  decreased <- decreasedSubset(data, species)
  sig_decreased <- sigDecreased(decreased)
  return(sig_decreased)}


# Test for edges of competition
identifyCompetition <- function(data1, data2) {
  group1 <- colnames(data1)
  group2 <- colnames(data2)
  keep <- intersect(colnames(data1), colnames(data2))
  group1 <- data1[,keep]
  group2 <- data2[,keep]
  competition <- rbind(group1, group2)
  return(competition)}

# Create subset of primary growth only
primarySubset <- function(data, species) {
  keep <- subset(data, R20291_secondary_growth == 'no')
  keep <- subset(keep, Depleting_species %in% c(species,'Fresh_media'))
  keep$R20291_secondary_growth <- NULL
  keep <- droplevels(keep)
  return(keep)}
# Create subset of secondary growth only
secondarySubset <- function(data, species) {
  keep <- subset(data, Depleting_species == species)
  keep$Depleting_species <- NULL
  keep <- droplevels(keep)
  return(keep)}

# Test for significant differences - Primary
sigIncreased <- function(data) {
  group1 <- subset(data, Depleting_species == 'Fresh_media')
  group1$Depleting_species <- NULL
  group2 <- subset(data, Depleting_species != 'Fresh_media')
  group2$Depleting_species <- NULL
  keep <- c()
  for (x in colnames(group1)) {
    test1 <- as.vector(group1[,x])
    test2 <- as.vector(group2[,x])
    pval <- as.numeric(wilcox.test(test1, test2, exact=FALSE)$p.value)
    if (is.na(pval)) {pval = 1}
    if (pval <= 0.05) {if (median(test1) < median(test2)) {keep <- c(keep, x)}}
  }
  group2 <- group2[, keep]
  return(group2)}

findIncrease <- function(data, species='default') {
  increased <- primarySubset(data, species)
  sig_increased <- sigIncreased(increased)
  return(sig_increased)}

# Test for significant differences - Secondary
testCrossfeeding <- function(data1, data2, species) {
  group1 <- subset(data1, R20291_secondary_growth == 'yes')
  group1$R20291_secondary_growth <- NULL
  group1 <- subset(group1, Depleting_species == species)
  group1$Depleting_species <- NULL
  keep <- c()
  for (x in colnames(data2)) {
    test1 <- as.vector(group1[,x])
    test2 <- as.vector(data2[,x])
    pval <- as.numeric(wilcox.test(test1, test2, exact=FALSE)$p.value)
    if (is.na(pval)) {pval = 1}
    if (pval <= 0.05) {if (median(test1) < median(test2)) {keep <- c(keep, x)}}
  }
  group1 <- group1[, keep]
  group2 <- data2[, keep]
  increased <- rbind(group1, group2)
  return(increased)}

#----------------------------------------------------------#

# Find consumed subtrates
Bproducta_decreased <- findDecrease(metabolome, 'B_producta_DSM2950')
Blongum_decreased <- findDecrease(metabolome, 'B_longum_ATCC55813')
Cdifficile_decreased <- findDecrease(metabolome, 'C_difficile_R20291')

# Check edges of competition
Bproducta_competition <- identifyCompetition(Bproducta_decreased, Cdifficile_decreased)
Blongum_competition <- identifyCompetition(Blongum_decreased, Cdifficile_decreased)
rm(Bproducta_decreased, Cdifficile_decreased, Blongum_decreased)

# Identify those unique to each organism
Bproducta_competition_only <- setdiff(colnames(Bproducta_competition), colnames(Blongum_competition))
Bproducta_competition_only <- Bproducta_competition[,Bproducta_competition_only]
Blongum_competition_only <- setdiff(colnames(Blongum_competition), colnames(Bproducta_competition))
Blongum_competition_only <- Blongum_competition[,Blongum_competition_only]
rm(Bproducta_competition, Blongum_competition)

#-------------------------------------------------------------#

# find produced substrates
Blongum_increased <- findIncrease(metabolome, 'B_longum_ATCC55813')
Bproducta_increased <- findIncrease(metabolome, 'B_producta_DSM2950')

# check for potential cross-feeding
Blongum_cooperation <- testCrossfeeding(metabolome, Blongum_increased, 'B_longum_ATCC55813')
Bproducta_cooperation <- testCrossfeeding(metabolome, Bproducta_increased, 'B_producta_DSM2950')
rm(Blongum_increased, Bproducta_increased)

# Identify those unique to each organism
Blongum_cooperation_only <- setdiff(colnames(Blongum_cooperation), colnames(Bproducta_cooperation))
Blongum_cooperation_only <- Blongum_cooperation[,Blongum_cooperation_only]
Bproducta_cooperation_only <- setdiff(colnames(Bproducta_cooperation), colnames(Blongum_cooperation))
Bproducta_cooperation_only <- Bproducta_cooperation[,Bproducta_cooperation_only]
rm(Blongum_cooperation, Bproducta_cooperation)

#-------------------------------------------------------#

# Assemble final tables
Blongum_competition_only <- Blongum_competition_only[1:6,]
Blongum_competition_fresh <- fresh[,colnames(Blongum_competition_only)]
Blongum_competition_only <- Blongum_competition_only - Blongum_competition_fresh
Blongum_competition_only <- log10(abs(Blongum_competition_only) + 1) * -1.0
Blongum_cooperation_only <- Blongum_cooperation_only[7:12,]
Blongum_cooperation_fresh <- fresh[,colnames(Blongum_cooperation_only)]
Blongum_cooperation_only <- Blongum_cooperation_only - Blongum_cooperation_fresh
Blongum_cooperation_only <- log10(abs(Blongum_cooperation_only) + 1)
Blongum <- as.data.frame(cbind(Blongum_competition_only, Blongum_cooperation_only))
rm(Blongum_competition_only, Blongum_competition_fresh, 
   Blongum_cooperation_only, Blongum_cooperation_fresh)

Bproducta_competition_only <- Bproducta_competition_only[1:6,]
Bproducta_competition_fresh <- fresh[,colnames(Bproducta_competition_only)]
Bproducta_competition_only <- Bproducta_competition_only - Bproducta_competition_fresh
Bproducta_competition_only <- log10(abs(Bproducta_competition_only) + 1) * -1.0
Bproducta_cooperation_only <- Bproducta_cooperation_only[7:12,]
Bproducta_cooperation_fresh <- fresh[,colnames(Bproducta_cooperation_only)]
Bproducta_cooperation_only <- Bproducta_cooperation_only - Bproducta_cooperation_fresh
Bproducta_cooperation_only <- log10(abs(Bproducta_cooperation_only) + 1)
Bproducta <- as.data.frame(cbind(Bproducta_competition_only, Bproducta_cooperation_only))
rm(Bproducta_competition_only, Bproducta_competition_fresh, 
   Bproducta_cooperation_only, Bproducta_cooperation_fresh)

# Reformat metabolite names
library(stringr)
Blongum_metabolites <- gsub('_', ' ', colnames(Blongum))
Blongum_metabolites <- gsub('\\.', '-', Blongum_metabolites)
Blongum_metabolites <- gsub('--', '-', Blongum_metabolites)
Blongum_metabolites <- gsub(' -IORA-', '', Blongum_metabolites)
Blongum_metabolites <- gsub(' -TENTATIVE-', '', Blongum_metabolites)
Blongum_metabolites <- gsub('^X', '', Blongum_metabolites)
Blongum_metabolites <- gsub('-2R-7-hydroxy-8-2-hydroxyethyl-5-methoxy-2-methyl-2-3-dihydrochromen-4-one', 'LL-D-253alpha', Blongum_metabolites)
Blongum_metabolites <- gsub('5-Aminoimidazole-4-carboxamide-1-beta-D-ribofuranosyl 5-monophosphate', 'AICAR', Blongum_metabolites)
Blongum_metabolites <- gsub('NCGC00380353-01 C15H18O6 Pentaleno-1-6a-c-pyran-9-carboxylic acid- 1-3-4-7-7a-9a-hexahydro-4-hydroxy-4-hydroxymethyl-6-7-dimethyl-3-oxo- -4R-7S-7aR-9aR-', 'Pentalenolactone O', Blongum_metabolites)
colnames(Blongum) <- Blongum_metabolites
Bproducta_metabolites <- gsub('_', ' ', colnames(Bproducta))
Bproducta_metabolites <- gsub('\\.', '-', Bproducta_metabolites)
Bproducta_metabolites <- gsub('--', '-', Bproducta_metabolites)
Bproducta_metabolites <- gsub(' -IORA-', '', Bproducta_metabolites)
Bproducta_metabolites <- gsub(' -TENTATIVE-', '', Bproducta_metabolites)
Bproducta_metabolites <- gsub('^X', '', Bproducta_metabolites)
Bproducta_metabolites <- gsub('NCGC00380373-01-2--2-amino-3-methylbutanoyl-amino-3-phenylpropanoic acid', 'MCULE-3661074999', Bproducta_metabolites)
Bproducta_metabolites <- gsub('NCGC00380353-01 C15H18O6 Pentaleno-1-6a-c-pyran-9-carboxylic acid- 1-3-4-7-7a-9a-hexahydro-4-hydroxy-4-hydroxymethyl-6-7-dimethyl-3-oxo- -4R-7S-7aR-9aR-', 'Pentalenolactone O', Bproducta_metabolites)
Bproducta_metabolites <- gsub('-2R-7-hydroxy-8-2-hydroxyethyl-5-methoxy-2-methyl-2-3-dihydrochromen-4-one', 'LL-D-253alpha', Bproducta_metabolites)
Bproducta_metabolites <- gsub('3-1-2-dihydroxypropyl-1-6-8-trihydroxyanthracene-9-10-dione', 'MCULE-9148758046', Bproducta_metabolites)
Bproducta_metabolites <- gsub('5-Aminoimidazole-4-carboxamide-1-beta-D-ribofuranosyl 5-monophosphate', 'AICAR', Bproducta_metabolites)
colnames(Bproducta) <- Bproducta_metabolites
rm(Blongum_metabolites, Bproducta_metabolites)

# Eliminate remaining duplicates
Blongum[,'LL-D-253alpha'] <- NULL
Blongum[,'AICAR'] <- NULL
Blongum[,'Pentalenolactone O'] <- NULL
Bproducta[,'LL-D-253alpha'] <- NULL
Bproducta[,'AICAR'] <- NULL
Bproducta[,'Pentalenolactone O'] <- NULL
colnames(Blongum) <- str_to_title(colnames(Blongum))
colnames(Bproducta) <- str_to_title(colnames(Bproducta))

#-------------------------------------------------------#

# Generate figure
set.seed(9)
library(pheatmap)
library(viridis)

pdf(file='~/Desktop/Jenior_Consortia_2022/results/Blongum_metabolomics.pdf', width=10, height=4)
pheatmap(as.matrix(Blongum), col=viridis(100), cluster_rows=FALSE, labels_row=c('','','','','',''))
dev.off()

pdf(file='~/Desktop/Jenior_Consortia_2022/results/Bproducta_metabolomics.pdf', width=10, height=4)
pheatmap(as.matrix(Bproducta), col=viridis(100), cluster_rows=FALSE, labels_row=c('','','','','',''))
dev.off()










