## Input data (gene x sample) needs to be TPM normalized.

library(immunedeconv)
library(data.table)
library(readr)
library(ggplot2)

## Data ------------------------------------------------
## CountMatrix: 24 pediatric renal transplant recipients
## 15 acute rejection of transplan
## 9 stable post-transplant course

CountMatrix <- fread("../../extdata/data/CountMatrix_WholeBloodSamples.txt")
GeneAnnotation <- fread("../../extdata/data/GeneAnnotation.txt")
colnames(CountMatrix)[1] <- "HGNC_Symbol"
colnames(GeneAnnotation)[c(1,2)] <- c("HGNC_ID", "HGNC_Symbol")

CountMatrix$HGNC_Symbol <- GeneAnnotation[match(CountMatrix$HGNC_Symbol, CountMatrix$HGNC_Symbol), external_gene_name]
CountMatrix <- as.matrix(CountMatrix, rownames = "HGNC_Symbol")

# wanted to get the sample metadata from here, but didnt manage manage at first
series <- read_file("../../extdata/data/GSE20300_series_matrix.txt")

# I think the data is log2 transformed -> we need to reverse the transformation
CountMatrix <- 2^CountMatrix

# --------------------------------------------------------
# quantiseq
res_quantiseq <- immunedeconv::deconvolute(CountMatrix, method = "quantiseq", arrays = T, tumor = F)
res_quantiseq <- melt(as.data.table(res_quantiseq), id.vars = "cell_type", variable.name = "sample", value.name = "proportion")
res_quantiseq$cell_type <- as.factor(res_quantiseq$cell_type)
ggplot(res_quantiseq, aes(x = sample, y = proportion, fill = cell_type)) + geom_col(position = "fill")
#--------------------------------------------------------
# timer
res_timer <- immunedeconv::deconvolute(CountMatrix, method = "timer", arrays = T, tumor = F)


#--------------------------------------------------------
# cibersort


#-------------------------------------------------------
# mcp_counter
res_mcp_counter <- immunedeconv::deconvolute(CountMatrix, method = "mcp_counter", arrays = T, tumor = F)
res_mcp_counter <- melt(as.data.table(res_mcp_counter), id.vars = "cell_type", variable.name = "sample", value.name = "score")
ggplot(res_mcp_counter, aes(x = sample, y = score, color = cell_type)) + geom_point() + facet_wrap(~cell_type, scales = "free")


#--------------------------------------------------------
# xcell

res_xcell <- immunedeconv::deconvolute(CountMatrix, method = "xcell", arrays = T, tumor = F)
res_xcell <- melt(as.data.table(res_xcell), id.vars = "cell_type", variable.name = "sample", value.name = "proportion")
ggplot(res_xcell, aes(x = sample, y = proportion, fill = cell_type)) + geom_col(position = "fill")


#----------------------------------------------------------
# epic
res_epic <- immunedeconv::deconvolute(CountMatrix, method = "epic", arrays = T, tumor = F)
res_epic <- melt(as.data.table(res_epic), id.vars = "cell_type", variable.name = "sample", value.name = "proportion")
ggplot(res_epic, aes(x = sample, y = proportion, fill = cell_type)) + geom_col(position = "fill")


###-------------------------------------------------------------------------------------------------------
# load real proportions
load("../../extdata/data/hw3/kidneyTransplant.RData")
cellFreq <- data.frame(cellFreq)
cellFreq["patientGroup"] <- as.factor(patientGroups)
names(cellFreq) <- c(cellNames, "patientGroup")
cellFreq <- melt(as.data.table(cellFreq), id.vars = "patientGroup", variable.name = "cell_type", value.name = "proportion")
levels(cellFreq$patientGroup) <- c("normal", "AR")

ggplot(cellFreq, aes(x = patientGroup, y = proportion, fill = cell_type)) + geom_col(position = "fill")
ggplot(cellFreq, aes(x = cell_type, y = proportion, fill = patientGroup)) + geom_col(position = "fill") 


## put all the data in one table
## quantiseq: combine cell types to lymphocytes
levels(res_quantiseq$cell_type)[res_quantiseq$cell_type %in% c("B cell", "NK cell", 
                                                               "T cell regulatory (Tregs)",
                                                               "T cell CD4+ (non-regulatory)",
                                                               "T cell CD8+")] <- "Lymphocytes"
levels(res_quantiseq$cell_type)[res_quantiseq$cell_type %in% c("Macrophage M1", "Macrophage M2", "Monocyte")] <- "Monocytes"
sample_desc = c("normal", "AR", "AR", "AR", "AR", "AR", "AR", "normal", "normal", "AR", "normal", "AR", "AR", "AR", "AR",
                "normal", "normal", "normal", "normal", "AR", "AR", "AR", "normal","AR")
tmp <- data.frame(cbind(levels(res_quantiseq$sample),sample_desc, deparse.level = 1))
names(tmp) <- c("sample", "group")
res_quantiseq <- merge(res_quantiseq, tmp, by = "sample", all.y = T)

ggplot(res_quantiseq, aes(group, proportion, fill = cell_type)) + geom_col(position = "fill") 
