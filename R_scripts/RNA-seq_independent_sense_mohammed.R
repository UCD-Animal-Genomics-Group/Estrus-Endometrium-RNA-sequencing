###############################################################
# RNA-seq data analysis of sense counts (independent samples) #
###############################################################

# Last updated: 20/04/2016

#############################
# List of required packages #
#############################

# If not already present, install packages
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("edgeR")
biocLite("org.Bt.eg.db")
install.packages("magrittr")
install.packages("dplyr")
install.packages("MASS")

# Load the required packages
library(limma)
library(edgeR)
library(org.Bt.eg.db)
library(magrittr)
library(dplyr)
library(MASS)

######################################################
# Use featureCounts output files as input files in R #
######################################################

###############
# Preparation #
###############

# Move to the appropriate folder
setwd("C:/Users/carol/Dropbox/CSF/Animal_Genomics/RNA-seq_Mohammed/edgeR/sense_genes")
getwd()
workDir <- getwd()
workDir
load("RNA-seq_independent_sense_mohammed.RData")

################################################
# Read in and concatenate input files within R #
################################################

# Create vector of all files name
fileDir <- "C:/Users/carol/Dropbox/CSF/Animal_Genomics/RNA-seq_Mohammed/edgeR/sense_genes/Counts"
files <- list.files(path = fileDir,
         pattern         = "^endo", 
		 all.files       = TRUE,
		 full.names      = FALSE,
		 recursive       = FALSE,
		 ignore.case     = FALSE)
files

# Reads and merges a set of files containing counts
Count <- readDGE(path         = fileDir,
                 files        = files,
                 header       = TRUE,
                 comment.char = "#",
                 columns      = c(1, 7))
names(Count)
head(Count$samples)
head(Count$counts)

# Ouptut data
write.table(x         = Count$samples,
            file      = "RNA-seq_Mohammed_samples.txt", 
			sep       = "\t",
			quote     = FALSE,
			row.names = TRUE,
			col.names = NA)
write.table(x         = Count$counts,
            file      = "RNA-seq_Mohammed_rawcounts.txt", 
			sep       = "\t",
			quote     = FALSE,
			row.names = TRUE,
			col.names = NA)

# After cleaning gene names from rawcounts file, read it in
raw.counts <- read.table(file   = "RNA-seq_Mohammed_rawcounts-clean.txt",
                         header = TRUE)
names(raw.counts)
head(raw.counts)

# After changing group numbers to letters and adding time point info,
# import file and prepare group information
group.table <- read.table(file   = "RNA-seq_Mohammed_samples-clean.txt",
                          header = TRUE)
names(group.table)
head(group.table)

group.info <- factor(paste0(group.table$group, "", group.table$time.point))
group.info

# Ouptut data
write.table(x         = group.table,
            file      = "RNA-seq_Mohammed_samples_groupInfo.txt", 
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

##################################################################
# Get gene information using the org.Bt.eg.db annotation package #
##################################################################

# Create annotation table with counts information
annotated.counts <- read.table(file   = "RNA-seq_Mohammed_rawcounts-clean.txt",
                               header = TRUE)
head(annotated.counts)
dim(annotated.counts)
columns(org.Bt.eg.db)

# Get gene names from NCBI gene identifiers
annotated.counts$gene.name <- mapIds(org.Bt.eg.db,
                              keys      = rownames(annotated.counts),
                              column    = "GENENAME",
                              keytype   = "ENTREZID",
                              multiVals = "first")
							  
# Get gene symbols from NCBI gene identifiers
annotated.counts$gene.symbol <- mapIds(org.Bt.eg.db,
                              keys      = rownames(annotated.counts),
                              column    = "SYMBOL",
                              keytype   = "ENTREZID",
                              multiVals = "first")

# Get ENSEMBL gene ids from NCBI gene identifiers
annotated.counts$ENSEMBL.tag <- mapIds(org.Bt.eg.db,
                                keys      = rownames(annotated.counts),
                                column    = "ENSEMBL",
                                keytype   = "ENTREZID",
                                multiVals = "first")

# Ouptut data
write.table(x         = annotated.counts,
            file      = "RNA-seq_Mohammed_annotated-counts.txt", 
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

##########################################
# Create a DGElist of all samples counts #
##########################################

# Select annotation columns from annotated.counts data frame
names(annotated.counts)
annotation <- select(annotated.counts, gene.name:ENSEMBL.tag)
head(annotation)

# Create DGElist containing the group and annotation
endo_dgelist <- DGEList(counts       = raw.counts,
                        group        = group.info,
                        genes        = annotation,
                        lib.size     = NULL,
                        norm.factors = NULL,
                        remove.zeros = FALSE)

names(endo_dgelist)
dim(endo_dgelist)
head(endo_dgelist$counts)
head(endo_dgelist$samples)
head(endo_dgelist$genes)

###########################################################
# Quality check of libraries by plotting density of count #
###########################################################

# Log10 transform the count data for better visualization
count_log10 <- log10(x = (endo_dgelist$counts
                          [, 1 : ncol(endo_dgelist$counts)] + 1))

# Plot density of raw counts for all libraries
png(filename = "Density_raw_Endo.png", width = 1366, height = 768, units = "px")

plot(x = density(count_log10[, 1]),
     main = "Density plot of raw counts per gene",
     lty  = 1, 
     xlab = "Log10 of raw counts per gene",
     ylab = "Density",
     col  = "black",
     ylim = c(0.0,1.25))


for (i in 2 : ncol(count_log10)) {
    lines(density(count_log10[, i]),
          lty = 1,
          col = "black")
}
dev.off()

#####################################
# Filtering of lowly expressed tags #
#####################################

# Filter non expressed tags (all genes that have zero counts in all samples)
Endo_no_zeros <- endo_dgelist[rowSums(endo_dgelist$counts) > 0, ]
dim(Endo_no_zeros$counts)
head(Endo_no_zeros$counts)
colnames(Endo_no_zeros$counts)

# Filter lowly expressed tags, retaining only tags with at least
# 1 count per million in 4 or more libraries
# (4 libraries correspond to the smallest treatment group)
Endo_filt <- Endo_no_zeros[rowSums(cpm(Endo_no_zeros) > 1) >= 4, ]
dim(Endo_filt$counts)
head(Endo_filt$counts)

# Output file of filtered counts
write.table(x         = Endo_filt$counts,
            file      = "Endo_filt_rawcounts.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

# Recompute the library size
Endo_filt$samples$lib.size <- colSums(Endo_filt$counts)
head(Endo_filt$samples)
head(endo_dgelist$samples)

##############################################################
# Normalization of data using trimmed mean of M-values       #
# (normalization based on RNA composition between libraries) #
##############################################################

# Calculate normalisation factor for our DGElist.
# With edgeR, counts are not transformed in any way after normalization,
# instead normalization will modify library size.
Endo_norm <- calcNormFactors(Endo_filt)
head(Endo_norm$samples)

#####################################################################
# Quality check of filtered libraries by plotting density of counts #
#####################################################################

# Log10 transform the filtered count data for better visualization
count_filt_log10 <- log10(Endo_norm$counts[, 1 : ncol(Endo_norm$counts)] + 1)

# Plot density of filtered count for all libraries
png(filename = "Density_Endo_post_filter.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plot(density(count_filt_log10[, 1]),
     main = "Density plot of count per gene",
     lty  = 1,
     xlab = "Log10 of count per gene",
     ylab = "Density",
     col  = "black",
     ylim = c(0.0,0.6))

for (i in 2 : ncol(count_filt_log10)) {
    lines(density(count_filt_log10[, i]),
          lty = 1,
          col = "black")
}
dev.off()

##################################
# Multidimensional scaling plots #
##################################

# set colors and shapes
colours = c('#8c510a', '#d8b365', '#f6e8c3', '#c7eae5', '#5ab4ac', '#01665e')
shapes = c(2, 4, 15, 16, 17, 18)

# Plot MDS of all samples
png(filename = "MDS_all_samples.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x      = Endo_norm,
        top    = 1000000,
        labels = NULL,
        pch    = shapes,
        col    = colours,
        xlab   = "Dimension 1",
        ylab   = "Dimension 2",
        cex    = 2)
legend("topright",
        Endo_norm,
        col    = colours,
        legend = levels(Endo_norm$samples$group),
        pch    = shapes,
        cex    = 1)

dev.off()

# Plot MDS of 12hrPostCIDR time point
png(filename = "MDS_12hrPostCIDR.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x      = Endo_norm[, grep(pattern = "^endo_12hrPostCIDR",
                                  x       = colnames(Endo_norm))],
        labels = NULL,
        pch    = 2,
        col    = '#8c510a',
        xlab   = "Dimension 1",
        ylab   = "Dimension 2",
        cex    = 2.0)
 
dev.off()

# Plot MDS of 24hrPostCIDR time point
png(filename = "MDS_24hrPostCIDR.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x      = Endo_norm
        [, grep(pattern = "^endo_24hrPostCIDR", x = colnames(Endo_norm))],
        labels = NULL,
        pch    = 4,
        col    = '#d8b365',
        xlab   = "Dimension 1",
        ylab   = "Dimension 2",
        cex    = 2.0)

dev.off()

# Plot MDS of Oestrus time point
png(filename = "MDS_Oestrus.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x      = Endo_norm
        [, grep(pattern = "^endo_Oestrus", x = colnames(Endo_norm))],
        labels = NULL,
        pch    = 15,
        col    = '#f6e8c3',
        xlab   = "Dimension 1",
        ylab   = "Dimension 2",
        cex    = 2.0)

dev.off()

# Plot MDS of 12hrPostOestrus time point
png(filename = "MDS_12hrPostOestrus.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x      = Endo_norm[, grep(pattern = "^endo_12hrPostOestrus",
                                  x       = colnames(Endo_norm))],
        labels = NULL,
        pch    = 16,
        col    = '#c7eae5',
        xlab   = "Dimension 1",
        ylab   = "Dimension 2",
        cex    = 2.0)

dev.off()

# Plot MDS of 48hrPostOestrus time point
png(filename = "MDS_48hrPostOestrus.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x      = Endo_norm[, grep(pattern = "^endo_48hrPostOestrus",
                                  x       = colnames(Endo_norm))],
        labels = NULL,
        pch    = 17,
        col    = '#5ab4ac',
        xlab   = "Dimension 1",
        ylab   = "Dimension 2",
        cex    = 2.0)

dev.off()

# Plot MDS of Luteal time point
png(filename = "MDS_Luteal.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x      = Endo_norm[, grep(pattern = "^endo_Luteal",
                                  x       = colnames(Endo_norm))],
        labels = NULL,
        pch    = 18,
        col    = '#01665e',
        xlab   = "Dimension 1",
        ylab   = "Dimension 2",
        cex    = 2.0)

dev.off()

###################################################
# Create a design matrix for independent analysis #
###################################################

# Create a design matrix
design <- model.matrix(~0 + group.info)
rownames(design) <- rownames(Endo_norm$samples)
colnames(design) <- levels(Endo_norm$samples$group)
design

# Output the design matrix info
write.table(x         = design,
            file      = "Endo_design.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)


########################################################################
# Estimate the dispersion parameter for each tag using Cox-Reid method #
# (for multi-factor data)                                              #
########################################################################

Endo_disp <- estimateGLMCommonDisp(y       = Endo_norm,
                                   design  = design,
                                   verbose = TRUE)
Endo_disp <- estimateGLMTrendedDisp(y      = Endo_disp,
                                    design = design)
Endo_disp <- estimateGLMTagwiseDisp(y      = Endo_disp,
                                    design = design)
names(Endo_disp)

# Plot the dispersion
png(filename = "BCV_Endo.png",
    width    = 1366,
    height   = 768,
    units    = "px")
plotBCV(Endo_disp)
dev.off()

# Show the calculated dispersion
Endo_disp$common.dispersion

# And show its square root, the coefficient of biological variation
sqrt(Endo_disp$common.dispersion)

# Create a matrix of the tagwise dispersion associated to ensembl gene
Tagwisedisp <- cbind(Endo_disp$genes, Endo_disp$tagwise.dispersion)
head(Tagwisedisp)
dim(Tagwisedisp)

# Write into a table the calculated tagwise dispersion
write.matrix(x    = Tagwisedisp,
             file = "Tagwise_dispersion.txt",
             sep  = "\t")

##################################################################
# Determine differential expression using negative binomial GLMs #
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
Endo_fit <- glmFit(y = Endo_disp, design = design)
colnames(Endo_fit$design)

# Make a contrast vector for pairwise comparisons between the groups
my.contrasts <- makeContrasts(PostCIDR12vsO = A12hrPostCIDR-COestrus,
                              PostCIDR24vsO = B24hrPostCIDR-COestrus,
                              OvsO          = COestrus-COestrus,
                              PostOes12vsO  = D12hrPostOestrus-COestrus,
                              PostOes48vsO  = E48hrPostOestrus-COestrus,
                              LutealvsO     = FLuteal-COestrus,
                              levels        = Endo_fit)

# Test for differential expression between the different time points/treatments
PostCIDR12vsO.lrt <- glmLRT(Endo_fit,
                            contrast = my.contrasts[,"PostCIDR12vsO"])
DE.12hrPostCIDR <- topTags(object        = PostCIDR12vsO.lrt,
                           n             = "inf",
                           adjust.method = "BH")

PostCIDR24vsO.lrt <- glmLRT(Endo_fit,
                            contrast = my.contrasts[,"PostCIDR24vsO"])
DE.24hrPostCIDR <- topTags(object        = PostCIDR24vsO.lrt,
                           n             = "inf",
                           adjust.method = "BH")

OvsO.lrt <- glmLRT(Endo_fit, contrast = my.contrasts[,"OvsO"])

PostOes12vsO.lrt <- glmLRT(Endo_fit,
                           contrast = my.contrasts[,"PostOes12vsO"])
DE.12hrPostOestrus <- topTags(object        = PostOes12vsO.lrt,
                              n             = "inf",
                              adjust.method = "BH")

PostOes48vsO.lrt <- glmLRT(Endo_fit,
                           contrast = my.contrasts[,"PostOes24vsO"])
DE.48hrPostOestrus <- topTags(object        = PostOes48vsO.lrt,
                              n             = "inf",
                              adjust.method = "BH")

LutealvsO.lrt <- glmLRT(Endo_fit,
                        contrast = my.contrasts[,"LutealvsO"])
DE.Luteal <- topTags(object        = LutealvsO.lrt,
                     n             = "inf",
                     adjust.method = "BH")

# Plot DE genes that passed FDR correction (p-value cut-off 0.05)
png(filename = "12hPostCIDR_vs_Oestrus.png",
    width    = 1366,
    height   = 768,
    units    = "px")
plotSmear(object  = PostCIDR12vsO.lrt,
          de.tags = (rownames(PostCIDR12vsO.lrt$table)
                     [as.logical(decideTestsDGE
                                 (PostCIDR12vsO.lrt,
                                  adjust.method = "BH",
                                  p.value       = 0.05))]))
abline(h = c(-1, 1), col = "blue")
dev.off()


png(filename = "24hPostCIDR_vs_Oestrus.png",
    width    = 1366,
    height   = 768,
    units    = "px")
plotSmear(object  = PostCIDR24vsO.lrt,
          de.tags = (rownames(PostCIDR24vsO.lrt$table)
                     [as.logical(decideTestsDGE
                                 (PostCIDR24vsO.lrt,
                                  adjust.method = "BH",
                                  p.value       = 0.05))]))
abline(h = c(-1, 1), col = "blue")
dev.off()


png(filename = "12hPostOestrus_vs_Oestrus.png",
    width    = 1366,
    height   = 768,
    units    = "px")
plotSmear(object  = PostOes12vsO.lrt,
          de.tags = (rownames(PostOes12vsO.lrt$table)
                     [as.logical(decideTestsDGE
                                 (PostOes12vsO.lrt,
                                  adjust.method = "BH",
                                  p.value       = 0.05))]))
abline(h = c(-1, 1), col = "blue")
dev.off()


png(filename = "48hPostOestrus_vs_Oestrus.png",
    width    = 1366,
    height   = 768,
    units    = "px")
plotSmear(object  = PostOes48vsO.lrt,
          de.tags = (rownames(PostOes48vsO.lrt$table)
                     [as.logical(decideTestsDGE
                                 (PostOes48vsO.lrt,
                                  adjust.method = "BH",
                                  p.value       = 0.05))]))
abline(h = c(-1, 1), col = "blue")
dev.off()


png(filename = "Luteal_vs_Oestrus.png",
    width    = 1366,
    height   = 768,
    units    = "px")
plotSmear(object  = LutealvsO.lrt,
          de.tags = (rownames(LutealvsO.lrt$table)
                     [as.logical(decideTestsDGE
                                 (LutealvsO.lrt,
                                  adjust.method = "BH",
                                  p.value       = 0.05))]))
abline(h = c(-1, 1), col = "blue")
dev.off()

#########################################################
# Merge all DE call data from the different time points #
#########################################################

# Output DE results for each comparison
write.table(x         = DE.12hrPostCIDR,
            file      = "DE_12hrPostCIDR.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE.24hrPostCIDR,
            file      = "DE_24hrPostCIDR.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE.12hrPostOestrus,
            file      = "DE_12hrPostOestrus.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE.48hrPostOestrus,
            file      = "DE_48hrPostOestrus.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE.Luteal,
            file      = "DE_Luteal.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

# Merge DE results from all comparisons into a single dataframe
Full_DE_endo <- merge(x  = DE.12hrPostCIDR$table,
                      y  = DE.24hrPostCIDR$table
                      [,(ncol(DE.24hrPostCIDR$table) - 4) :
                           ncol(DE.24hrPostCIDR$table)],
                      by = "row.names")
head(Full_DE_endo)

rownames(Full_DE_endo) <- Full_DE_endo[, 1] # To correct row names before
Full_DE_endo <- Full_DE_endo[, -1]          # merging the other data frames.

colnames(Full_DE_endo) <- gsub(pattern     = ".x$",
                               replacement = "_12hrPostCIDR",
                               x           = colnames(Full_DE_endo),
                               perl        = TRUE)
colnames(Full_DE_endo) <- gsub(pattern     = ".y$",
                               replacement = "_24hrPostCIDR",
                               x           = colnames(Full_DE_endo),
                               perl        = TRUE)
head(Full_DE_endo)

Full_DE_endo <- merge(x  = Full_DE_endo,
                      y  = DE.12hrPostOestrus$table
                      [, (ncol(DE.12hrPostOestrus$table) - 4) :
                           ncol(DE.12hrPostOestrus$table)],
                      by = "row.names")
head(Full_DE_endo)

rownames(Full_DE_endo) <- Full_DE_endo[, 1]
Full_DE_endo <- Full_DE_endo[, -1]

Full_DE_endo <- merge(x  = Full_DE_endo,
                      y  = DE.48hrPostOestrus$table
                      [, (ncol(DE.48hrPostOestrus$table) - 4) :
                           ncol(DE.48hrPostOestrus$table)],
                      by = "row.names")
head(Full_DE_endo)

rownames(Full_DE_endo) <- Full_DE_endo[, 1]
Full_DE_endo <- Full_DE_endo[, -1]

colnames(Full_DE_endo) <- gsub(pattern     = ".x$",
                               replacement = "_12hrPostOestrus",
                               x           = colnames(Full_DE_endo),
                               perl        = TRUE)
colnames(Full_DE_endo) <- gsub(pattern     = ".y$",
                               replacement = "_48hrPostOestrus",
                               x           = colnames(Full_DE_endo),
                               perl        = TRUE)
head(Full_DE_endo)

Full_DE_endo <- merge(x  = Full_DE_endo,
                      y  = DE.Luteal$table[, (ncol(DE.Luteal$table) - 4) :
                                               ncol(DE.Luteal$table)],
                      by = "row.names")
head(Full_DE_endo)

rownames(Full_DE_endo) <- Full_DE_endo[, 1]
head(Full_DE_endo)

colnames(Full_DE_endo) <- gsub(pattern     = "^logFC$",
                               replacement = "logFC_Luteal",
                               x           = colnames(Full_DE_endo),
                               perl        = TRUE)
colnames(Full_DE_endo) <- gsub(pattern     = "^logCPM$",
                               replacement = "logCPM_Luteal",
                               x           = colnames(Full_DE_endo),
                               perl        = TRUE)
colnames(Full_DE_endo) <- gsub(pattern     = "^LR$",
                               replacement = "LR_Luteal",
                               x           = colnames(Full_DE_endo),
                               perl        = TRUE)
colnames(Full_DE_endo) <- gsub(pattern     = "^PValue$",
                               replacement = "PValue_Luteal",
                               x           = colnames(Full_DE_endo),
                               perl        = TRUE)
colnames(Full_DE_endo) <- gsub(pattern     = "^FDR$",
                               replacement = "FDR_Luteal",
                               x           = colnames(Full_DE_endo),
                               perl        = TRUE)
colnames(Full_DE_endo)[1] <- "RefSeq_gene_id"
head(Full_DE_endo)

write.table(x         = Full_DE_endo,
            file      = "Full_DE_endo.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

####################
# Save .RData file #
####################

save.image(file = "RNA-seq_independent_sense_mohammed.RData")

