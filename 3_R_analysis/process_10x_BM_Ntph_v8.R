rm(list = ls())

setwd('/Users/minhookim/Data/scRNA-seq/Neutrophil/Code/R_MKim')
options(stringsAsFactors = F)

library('Seurat')
library(SeuratWrappers)
library(SingleR) 
library(singleCellNet)
library(Augur)
library(viridis)
library('muscat')
library(sctransform)
library("singleCellTK")
library(celldex)
library(monocle3)

########################################################################################################################################################
# 0. Read raw data

# read HTO libraries
M_HTO <- Read10X('/Users/minhookim/Data/scRNA-seq/Neutrophil/Code/2_CITEseq_tools/Output/3m_M_HTO_parsed/umi_count/', gene.column=1)
F_HTO <- Read10X('/Users/minhookim/Data/scRNA-seq/Neutrophil/Code/2_CITEseq_tools/Output/3m_F_HTO_parsed/umi_count/', gene.column=1)

M_HTO <- M_HTO[1:4,] # remove unmapped
F_HTO <- F_HTO[1:4,] # remove unmapped

# need to clean up cells with no HTOs
m.zero <- which(colSums(as.matrix(M_HTO)) < 100)
f.zero <- which(colSums(as.matrix(F_HTO)) < 100)

M_HTO <- M_HTO[,-m.zero] # remove cells with no hashtag
F_HTO <- F_HTO[,-f.zero] # remove cells with no hashtag

# read 10x library
ntph_F <- Read10X('/Users/minhookim/Data/scRNA-seq/Neutrophil/Code/1_Cell_Ranger/3m_F_NTPH_v2/outs/raw_feature_bc_matrix/')
ntph_M <- Read10X('/Users/minhookim/Data/scRNA-seq/Neutrophil/Code/1_Cell_Ranger/3m_M_NTPH_v2/outs/raw_feature_bc_matrix/')

dim(ntph_F)
# 32285 433747
dim(ntph_M)
# 32285 467210

########################################################################################################################################################
# 1. Filter and demultiplex cells

###################################################################################################################
# a. Select cell barcodes detected by both RNA and HTO 
joint.bcs.F <- intersect(paste0(colnames(F_HTO),"-1"), colnames(ntph_F))
joint.bcs.M <- intersect(paste0(colnames(M_HTO),"-1"), colnames(ntph_M))

# Clean names to extract HTO from each lane
joint.bcs.F.cl <- unlist(strsplit(joint.bcs.F,"-1"))
joint.bcs.M.cl <- unlist(strsplit(joint.bcs.M,"-1"))

# Subset RNA and HTO counts by joint cell barcodes
ntph.F.umis <- ntph_F[, joint.bcs.F]
ntph.M.umis <- ntph_M[, joint.bcs.M]

# Select joint HTOs for future merging
ntph.htos.F           <- as.matrix(F_HTO[, joint.bcs.F.cl])
ntph.htos.M           <- as.matrix(M_HTO[, joint.bcs.M.cl])
colnames(ntph.htos.F) <- paste0("Female_NTPH_",joint.bcs.F) # to compare with RNA object
colnames(ntph.htos.M) <- paste0("Male_NTPH_"  ,joint.bcs.M) # to compare with RNA object

# Make combination for addition to Seurat
ntph.htos <- cbind(ntph.htos.F,ntph.htos.M)

# Confirm that the HTO have the correct names
rownames(ntph.htos)
#  "HTO1_3m_F_1-ACCCACCAGTAAGAC" "HTO2_3m_M_1-GGTCGAGAGCATTCA" "HTO3_3m_F_2-CTTGCCGCATGTCAT" "HTO4_3m_M_2-AAAGCATTCTTCACG"

###########################
# b. Setup Seurat object

seurat.F <- CreateSeuratObject(counts = ntph.F.umis, project = "10X_ntph_F")
seurat.M <- CreateSeuratObject(counts = ntph.M.umis, project = "10X_ntph_M")

ntph.combined <- merge(seurat.F,
                       y =  seurat.M,
                       add.cell.ids = c("Female_NTPH" ,
                                        "Male_NTPH"),
                       project = "10X_ntph_Young")
ntph.combined
# An object of class Seurat 
# 32285 features across 6264 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ntph.combined@assays$RNA)
num.cells <- Matrix::rowSums(ntph.combined@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ntph.filt <- subset(ntph.combined, features = genes.use)
ntph.filt
# An object of class Seurat
# 11199 features across 6264 samples within 1 assay 
# Active assay: RNA (11199 features, 0 variable features)

# Normalize RNA data with log normalization
ntph.filt <- NormalizeData(ntph.filt, normalization.method = "LogNormalize",  scale.factor = 10000)

# Find and scale variable features
ntph.filt <- FindVariableFeatures(ntph.filt, selection.method = "mean.var.plot")
ntph.filt <- ScaleData(ntph.filt, features = VariableFeatures(ntph.filt))

#########################
# b. Add HTO info to Seurat object

# Add HTO data as an independent assay from RNA
ntph.filt[["HTO"]] <- CreateAssayObject(counts = ntph.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ntph.filt <- NormalizeData(ntph.filt, assay = "HTO", normalization.method = "CLR")

# Demultiplex cells based on HTO enrichment
# https://github.com/satijalab/seurat/issues/2549
ntph.filt <- MULTIseqDemux(ntph.filt, assay = "HTO") 

# We can visualize how many cells are classified as singlets, doublets and negative/ambiguous cells.
table(ntph.filt$MULTI_ID)
# Doublet HTO1-3m-F-1-ACCCACCAGTAAGAC HTO2-3m-M-1-GGTCGAGAGCATTCA HTO3-3m-F-2-CTTGCCGCATGTCAT HTO4-3m-M-2-AAAGCATTCTTCACG 
# 66                        1460                        1843                        1237                        1656 
# Negative 
# 2 

# Visualize enrichment for selected HTOs with ridge plots, and group cells based on the max HTO signal
Idents(ntph.filt) <- "MULTI_ID"

pdf(paste0(Sys.Date(),"_HTO_id_ridge_plot.pdf"), height = 10, width = 15)
RidgePlot(ntph.filt, assay = "HTO", features = rownames(ntph.filt[["HTO"]]), ncol = 2)
dev.off()

pdf(paste0(Sys.Date(),"_HTO1_2_feature_Scatter_plot.pdf"), height = 5, width = 10)
FeatureScatter(ntph.filt, feature1 = "hto_HTO1-3m-F-1-ACCCACCAGTAAGAC", feature2 = "hto_HTO2-3m-M-1-GGTCGAGAGCATTCA")
dev.off()

# Compare number of UMIs for singlets, doublets and negative cells
pdf(paste0(Sys.Date(),"_nCount_RNA_violin_based_on_hashtag.pdf"), height = 10, width = 15)
VlnPlot(ntph.filt, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()


# Generate a two dimensional embedding for HTOs. Here we are grouping cells by singlets and doublets for simplicity.
# First, we will remove negative cells from the object
ntph.filt.subset <- subset(ntph.filt, idents = "Negative", invert = TRUE)


# Calculate a UMAP embedding of the HTO data
DefaultAssay(ntph.filt.subset) <- "HTO"
ntph.filt.subset <- ScaleData(ntph.filt.subset, features = rownames(ntph.filt.subset), verbose = FALSE)
ntph.filt.subset <- RunPCA(ntph.filt.subset   , features = rownames(ntph.filt.subset), approx = FALSE)
ntph.filt.subset <- RunUMAP(ntph.filt.subset  , dims = 1:4, perplexity = 100, check_duplicates = FALSE) # only 4 dimension possible

pdf(paste0(Sys.Date(),"_UMAP_for_hashtag_levels.pdf"), height = 10, width = 15)
DimPlot(ntph.filt.subset)
dev.off()

# create metadata column to plot HTO heatmap
ntph.filt@meta.data$HTO.global <- rep("Singlet",dim(ntph.filt@meta.data)[1])
ntph.filt@meta.data$HTO.global[ntph.filt@meta.data$MULTI_ID %in% "Doublet"] <- "Doublet"
ntph.filt@meta.data$HTO.global[ntph.filt@meta.data$MULTI_ID %in% "Negative"] <- "Negative"

# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.
pdf(paste0(Sys.Date(),"_HTO_levels_heatmap.pdf"), height = 10, width = 15)
HTOHeatmap(ntph.filt, assay = "HTO",   
           classification = "MULTI_ID",
           global.classification = "HTO.global")
dev.off()

########################################################################################################################################################
# 2. Extract singlets and run QC

# Extract the singlets
Idents(ntph.filt) <- "HTO.global"
ntph.singlet <- subset(ntph.filt, idents = "Singlet")
ntph.singlet
# An object of class Seurat 
# 11203 features across 6196 samples within 2 assays 
# Active assay: RNA (11199 features, 730 variable features)
# 1 other assay present: HTO

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.  
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
ntph.singlet[["percent.mito"]] <- PercentageFeatureSet(ntph.singlet, pattern = "^mt-")
head(ntph.singlet@meta.data)
#                                 orig.ident nCount_RNA nFeature_RNA nCount_HTO nFeature_HTO                    MULTI_ID        MULTI_classification HTO.global
# Female_NTPH_GCTTTCGGTTCTTCAT-1 10X_ntph_F       5481         1576        127            2 HTO3-3m-F-2-CTTGCCGCATGTCAT HTO3-3m-F-2-CTTGCCGCATGTCAT    Singlet
# Female_NTPH_GATGGAGGTCGCTTGG-1 10X_ntph_F       2035          851        164            1 HTO1-3m-F-1-ACCCACCAGTAAGAC HTO1-3m-F-1-ACCCACCAGTAAGAC    Singlet
# Female_NTPH_AGTCTCCAGTCTTCCC-1 10X_ntph_F       8253         2045        199            4 HTO3-3m-F-2-CTTGCCGCATGTCAT HTO3-3m-F-2-CTTGCCGCATGTCAT    Singlet
# Female_NTPH_GGAGCAACAGGTCTCG-1 10X_ntph_F      10749         2347        749            2 HTO1-3m-F-1-ACCCACCAGTAAGAC HTO1-3m-F-1-ACCCACCAGTAAGAC    Singlet
# Female_NTPH_ATGGATCCAGACACCC-1 10X_ntph_F       3697         1418        387            2 HTO1-3m-F-1-ACCCACCAGTAAGAC HTO1-3m-F-1-ACCCACCAGTAAGAC    Singlet
# Female_NTPH_GGGAGTACAGATACCT-1 10X_ntph_F       4062         1426        195            3 HTO3-3m-F-2-CTTGCCGCATGTCAT HTO3-3m-F-2-CTTGCCGCATGTCAT    Singlet
#                                percent.mito
# Female_NTPH_GCTTTCGGTTCTTCAT-1   0.23718300
# Female_NTPH_GATGGAGGTCGCTTGG-1   0.04914005
# Female_NTPH_AGTCTCCAGTCTTCCC-1   0.19386890
# Female_NTPH_GGAGCAACAGGTCTCG-1   0.18606382
# Female_NTPH_ATGGATCCAGACACCC-1   0.13524479
# Female_NTPH_GGGAGTACAGATACCT-1   0.09847366

pdf(paste(Sys.Date(),"NTPH_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ntph.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

pdf(paste(Sys.Date(),"NTPH_violinPlots_QC_gene_UMI_mito_by_sample.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ntph.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "MULTI_ID", col = c("deeppink","deepskyblue","deeppink","deepskyblue"))
dev.off()

# plot important QC metrics
plot1 <- FeatureScatter(ntph.singlet, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ntph.singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"NTPH_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2))
dev.off()

# get clean cells
ntph.singlet <- subset(ntph.singlet, subset = nFeature_RNA > 100 & percent.mito < 25)
ntph.singlet
# An object of class Seurat 
# 11203 features across 6073 samples within 2 assays 
# Active assay: RNA (11199 features, 730 variable features)
# 1 other assay present: HTO


pdf(paste(Sys.Date(),"NTPH_violinPlots_QC_gene_UMI_mito_by_sample_POST_FILTER.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ntph.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "MULTI_ID", col = c("deeppink","deepskyblue","deeppink","deepskyblue"))
dev.off()

# Apply sctransform to dataset
ntph.singlet <- SCTransform(object = ntph.singlet, vars.to.regress = c("nFeature_RNA", "percent.mito"))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ntph.singlet), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ntph.singlet)
plot2 <- LabelPoints(plot = plot1, points = top10)

pdf(paste(Sys.Date(),"ntph_Singlets_Variable genes.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

# Run PCA
ntph.singlet <- RunPCA(ntph.singlet)
ElbowPlot(ntph.singlet)  #15

# We select the top 15 PCs for clustering and UMAP based on PCElbowPlot
ntph.singlet <- RunUMAP(ntph.singlet, reduction = "pca", dims = 1:15)
ntph.singlet <- FindNeighbors(ntph.singlet, reduction = "pca", dims = 1:15)
ntph.singlet <- FindClusters(ntph.singlet, resolution = 0.3, verbose = FALSE)

# Projecting singlet identities on UMAP visualization
pdf(paste(Sys.Date(),"ntph_10xRNA_UMAP_by_HTO.pdf", sep = "_"), height = 5, width = 10)
DimPlot(ntph.singlet, reduction = "umap", group.by = "MULTI_ID", cols = c("pink","dodgerblue","deeppink","skyblue"))
dev.off()

pdf(paste(Sys.Date(),"ntph_UMAP_clusters_0.3res.pdf", sep = "_"), height = 5, width = 10)
DimPlot(ntph.singlet, reduction = "umap", group.by = "seurat_clusters")
dev.off()

# plot gene expression of key neutrophil markers and sex differential genes
pdf(paste(Sys.Date(),"ntph_UMAP_Ly6g_Cd11b.pdf", sep = "_"), height = 5, width = 10)
FeaturePlot(object = ntph.singlet,features = c("Ly6g","Itgam"), pt.size = 0.5)
dev.off()

pdf(paste(Sys.Date(),"ntph_UMAP_Sex_Gene_markers.pdf", sep = "_"), height = 5, width = 10)
FeaturePlot(object = ntph.singlet,features = c("Xist","Ddx3y"), pt.size = 0.5)
dev.off()

pdf(paste(Sys.Date(),"ntph_Ridgeplot_Sex_Gene_markers.pdf", sep = "_"), height = 5, width = 12)
RidgePlot(object = ntph.singlet,features = c("Xist","Ddx3y"), group.by = "MULTI_ID", cols = c("pink","dodgerblue","deeppink","skyblue"))
dev.off()

# save data
save(ntph.singlet, file = paste(Sys.Date(),"Seurat_10x_BM_Ntph_withHTO_Singlets.RData",sep = "_"))

########################################################################################################################################################
# 3. SingleR (Immgen) cell annotation

ntph.singlet@meta.data$Sample_ID <- rep('NA')
ntph.singlet@meta.data$Sample_ID[grep("3m-F-1",ntph.singlet@meta.data$MULTI_ID)] <- "3m_F_1"
ntph.singlet@meta.data$Sample_ID[grep("3m-F-2",ntph.singlet@meta.data$MULTI_ID)] <- "3m_F_2"
ntph.singlet@meta.data$Sample_ID[grep("3m-M-1",ntph.singlet@meta.data$MULTI_ID)] <- "3m_M_1"
ntph.singlet@meta.data$Sample_ID[grep("3m-M-2",ntph.singlet@meta.data$MULTI_ID)] <- "3m_M_2"

table(ntph.singlet@meta.data$Sample_ID)
# 3m_F_1 3m_F_2 3m_M_1 3m_M_2 
# 1406   1219   1814   1634 

###########################
# a. Run singleR and save object

my.SingleCellExperiment.object <- as.SingleCellExperiment(ntph.singlet)

immgen <- ImmGenData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))
immgen
# class: SummarizedExperiment 
# dim: 22134 830 
# metadata(0):
#        assays(1): logcounts
# rownames(22134): Zglp1 Vmn2r65 ... Tiparp Kdm1a
# rowData names(0):
#        colnames(830): GSM1136119_EA07068_260297_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_1.CEL GSM1136120_EA07068_260298_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_2.CEL
# ... GSM920654_EA07068_201214_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_1.CEL GSM920655_EA07068_201215_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_2.CEL
# colData names(3): label.main label.fine label.ont

my.singler.immgen = SingleR(test = my.SingleCellExperiment.object,
                            ref  = immgen,
                            assay.type.test = 1,
                            labels = immgen$label.main)

table(my.singler.immgen$labels)
# B cells B cells, pro    Basophils           DC          ILC    Monocytes  Neutrophils   Stem cells      T cells 
# 4            6            3            4            2           16         6025           10            3 

# Add Sample IDs to metadata
my.singler.immgen$meta.data$Sample       =    ntph.singlet@meta.data$Sample_ID                       # Sample IDs

save(my.singler.immgen, file = paste(Sys.Date(),"10x_BM_Ntph_SingleR_object_SCT_Immgen.RData",sep = "_"))


###########################
# b. Transfer cell annotations to Seurat object

my.SingleCellExperiment.object <- as.SingleCellExperiment(ntph.singlet)

# Plot annotation heatmap
annotation_col = data.frame(Labels = factor(my.singler.immgen$labels),
        					Sample_ID = as.data.frame(colData(my.SingleCellExperiment.object)[,"Sample_ID",drop=FALSE]))

ann_colors = list(Labels = c("B cells" = "red",
							 "B cells, pro" = "red4",
               			     "Basophils" = "darkorange",
            		         "DC" = "gold",
    		                 "ILC" = "lightgoldenrod",
             			     "Monocytes" = "palegreen",
        		             "Neutrophils" = "seagreen",
           			         "Stem cells" = "midnightblue",
       			             "T cells" = "purple"),
   			     Sample_ID = c("3m_F_1" = "#EFA8CA", "3m_F_2" = "#ED4C8A", "3m_M_1" = "#73CFEC", "3m_M_2" = "#414DA0"))

pdf(paste(Sys.Date(),"10X_NTPH_SingleR_object_heatmap_Immgen.cell.types.pdf",sep = "_"), height = 5, width = 7)
plotScoreHeatmap(my.singler.immgen, annotation_col=annotation_col, cluster_cols = TRUE, annotation_colors = ann_colors)
dev.off()

# Transfer SingleR annotations to Seurat Object
ntph.singlet[["SingleR.labels"]] <- my.singler.immgen$labels

table(ntph.singlet@meta.data$SingleR.labels )
# B cells B cells, pro    Basophils           DC          ILC    Monocytes  Neutrophils   Stem cells      T cells 
# 4            6            3            4            2           16         6025           10            3 

# Save annotation
write.table(my.singler.immgen$labels, file = paste(Sys.Date(),"10X_NTPH_SingleR_annotated_cells_Immgen.txt",sep = "_"), sep = "\t", quote = F)

# Plot tSNE colored by SingleR cell annotation
pdf(paste(Sys.Date(),"10X_NTPH_SingleR_tSNE_annotated_IMMGEN.pdf", sep = "_"), height = 7, width = 7)
DimPlot(ntph.singlet, reduction = "umap", group.by = "SingleR.labels", cols = c("blue3",
                                                                                "dodgerblue",
                                                                                "cyan3",
                                                                                "chartreuse3",
                                                                                "deeppink2",
                                                                                "coral1",
                                                                                "darkgoldenrod1",
                                                                                "aquamarine2",
                                                                                "firebrick1"))
dev.off()

###########################
# c. Extract neutrophils to obtain clean neutrophils Seurat object

ntph.singlet.cl <- subset(ntph.singlet, subset = SingleR.labels %in% "Neutrophils")    # Neutrophils, 6025 cells 
ntph.singlet.cl
# An object of class Seurat 
# 22402 features across 6025 samples within 3 assays 
# Active assay: SCT (11199 features, 3000 variable features)
#  2 other assays present: RNA, HTO
#  2 dimensional reductions calculated: pca, umap

save(ntph.singlet.cl, file = paste(Sys.Date(),"10x_BM_Ntph_Singlets_Seurat_SingleR_QC.object.RData",sep = "_"))

########################################################################################################################################################
# 4. Neutrophil subpopulation annotation using singleCellNet - Xie et al. dataset

# Load query data
ntph.singlet.cl
# An object of class Seurat 
# 22402 features across 6025 samples within 3 assays 
# Active assay: SCT (11199 features, 3000 variable features)
#  2 other assays present: RNA, HTO
#  2 dimensional reductions calculated: pca, umap

# Load reference data
# Import seurat object from GEO on droplets
# Get Xie NTPH object
load("./2020-08-31_NTPH_Xie_paper_Seurat_obj.RData")

Xie.combined.filt
# An object of class Seurat 
# 27998 features across 12285 samples within 1 assay 
# Active assay: RNA (27998 features, 0 variable features)

# Remake object to have clean metadata slots
xie.clean <- CreateSeuratObject(counts = Xie.combined.filt@assays$RNA@counts, assay = "RNA", meta.data = Xie.combined.filt@meta.data)

# Subset to remove "non QC" cells (Cont and GM) based on paper
table(xie.clean@meta.data$cluster)
# Cont   G0   G1   G2   G3   G4  G5a  G5b  G5c   GM 
#  130  509  436  699 2496 1837 3545  823 1647  163 

xie.clean.2 <- subset(xie.clean, 
                      subset = cluster %in% c("G0",
                                              "G1",
                                              "G2",
                                              "G3",
                                              "G4",
                                              "G5a",
                                              "G5b",
                                              "G5c"))   # QC cells, 11992 cells

########################## TRAINING ##############################
# Extract info from reference object for query data
# exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
XIE.seuratfile    <- extractSeurat(xie.clean.2, exp_slot_name = "counts")
XIE.sampTab       <- XIE.seuratfile$sampTab
XIE.expDat        <- XIE.seuratfile$expDat
XIE.sampTab       <- droplevels(XIE.sampTab)
XIE.sampTab$cell  <- rownames(XIE.sampTab)

# Find genes in common to the data sets and limit analysis to these

genes.ntph <- rownames(ntph.singlet.cl)

commonGenes <- intersect(rownames(XIE.expDat), genes.ntph)
length(commonGenes)
# [1] 10250

XIE.expDat <- XIE.expDat[commonGenes,]

# Split for training and assessment, and transform training data
set.seed(123456789)
stList   <- splitCommon(sampTab=XIE.sampTab, ncells=100, dLevel="cluster")
stTrain  <- stList[[1]]
expTrain <- XIE.expDat[,rownames(stTrain)]

# Train the classifier using the Xie et al. data
class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, 
                        nTopGenes = 100, nTopGenePairs = 50, nRand = 50, nTrees = 1000, 
                        dLevel = "cluster", colName_samp = "cell")
# There are 397 top gene pairs

########################## TESTING ##############################
# Assessing the classifier with heldout data Apply to held out data
stTestList <- splitCommon(sampTab=stList[[2]], ncells=100, dLevel="cluster") # normalize validation data so that the assessment is as fair as possible
stTest     <- stTestList[[1]]
expTest    <- XIE.expDat[commonGenes,rownames(stTest)]

# Predict on held out data
classRes_val_all <- scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 100)

# Assess classifier performance on held out data
tm_heldoutassessment <- assess_comm(ct_scores  = classRes_val_all, 
                                    stTrain    = stTrain, 
                                    stQuery    = stTest, 
                                    dLevelSID  = "cell", 
                                    classTrain = "cluster", 
                                    classQuery = "cluster", 
                                    nRand = 100)

pdf(paste0(Sys.Date(),"_PR_curves_SingleCellNet_Xie_model_performance_ntph.pdf") )
plot_PRs(tm_heldoutassessment)
dev.off()

pdf(paste0(Sys.Date(),"_AUPRC_curves_SingleCellNet_Xie_model_performance_ntph.pdf") )
plot_metrics(tm_heldoutassessment)
dev.off()

# Classification result heatmap
# Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand          = 100
sla            = as.vector(stTest$cluster)
names(sla)     = as.vector(stTest$cell)
slaRand        = rep("rand", nrand)
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla            = append(sla, slaRand) # include in the random cells profile created

# Attribution plot
pdf(paste0(Sys.Date(),"_Classification_result_barplot_SingleCellNet_Tabula_Muris_Senis_IMMUNE_model_performance_ntph.pdf"), height = 5, width = 10 )
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="cluster", sid="cell")
dev.off()
# Pretty good

########################## QUERY ##############################
# Extract info from Seurat object for query data
# exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
ntph.seuratfile <- extractSeurat(ntph.singlet.cl, exp_slot_name = "counts")
ntph.sampTab    <- ntph.seuratfile$sampTab
ntph.expDat     <- ntph.seuratfile$expDat

classRes_ntph <- scn_predict(class_info[['cnProc']], ntph.expDat, nrand=50)

# Classification annotation assignment
# This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold of your choosing.
# The annotation result can be found in a column named category in the query sample table.

# Need to modify dLevel to SingleR.labels
stNtph <- get_cate(classRes = classRes_ntph, sampTab = ntph.sampTab, dLevel = "SingleR.labels", sid = "Sample_ID", nrand = 50)

table(stNtph$category,stNtph$Sample_ID)
#     3m_F_1 3m_F_2 3m_M_1 3m_M_2
# G2       8     14     20     25
# G3     364    325    429    377
# G4     697    525    838    783
# G5a    261    255    376    343
# G5b     22     27     27     16
# G5c     44     59    112     78


# Transfer annotation to SingleR data
sum(rownames(ntph.singlet.cl@meta.data) == rownames(stNtph)) # 6025
# sum(rownames(ntph.singlet.cl@meta.data) == sample(rownames(stNtph))) # 0

# dim(ntph.singlet.cl@meta.data) # 6025   15

ntph.singlet.cl@meta.data$SingleCellNet_Xie <- stNtph$category

pdf(paste0(Sys.Date(),"_BM_Neutrophils_10X_SingleCellNet_Xie_UMAP_annotated.pdf"), height = 5, width = 6 )
DimPlot(ntph.singlet.cl, reduction = "umap", group.by = "SingleCellNet_Xie",
        cols = c("blue3",
                   "dodgerblue",
                   "chartreuse3",
                   "darkgoldenrod1",
                   "firebrick1",
                   "darkviolet"))
dev.off()

save(ntph.singlet.cl, file = paste(Sys.Date(),"10x_BM_Ntph_USC_Xie_SingleCellNet_predictions_Seurat_object.RData",sep = "_"))

###########################
# Assess annotation

my.freq.table <- prop.table(x = table(ntph.singlet.cl@meta.data$SingleCellNet_Xie, ntph.singlet.cl@meta.data$Sample_ID), margin = 2)
my.freq.table.av <- apply(my.freq.table,1,mean)
my.freq.table.av.sort <- sort(my.freq.table.av, decreasing = T,index.return = T)

my.cols <- c("blue3",
             "dodgerblue",
             "chartreuse3",
             "darkgoldenrod1",
             "firebrick1",
             "darkviolet")

pdf(paste(Sys.Date(),"10XGenomics_BM_Neutrophils_USC_Xie_SingleCellNet_predictions_barplot_FREQUENCY.pdf",sep = "_"), height = 10, width = 15)
barplot(my.freq.table[my.freq.table.av.sort$ix,], las = 2,
        legend.text = rownames(my.freq.table)[my.freq.table.av.sort$ix],
        col = my.cols[my.freq.table.av.sort$ix],
        ylab = "Cell frequency in 10xGenomics single cell library")
box()
dev.off()

pdf(paste(Sys.Date(),"10XGenomics_BM_Neutrophils_USC_Xie_SingleCellNet_predictions_barplot_FREQUENCY_by_sex.pdf",sep = "_"), height = 10, width = 15)
barplot(my.freq.table[my.freq.table.av.sort$ix,], las = 2,
        legend.text = rownames(my.freq.table)[my.freq.table.av.sort$ix],
        col = my.cols[my.freq.table.av.sort$ix],
        ylab = "Cell frequency in 10xGenomics single cell library")
box()
dev.off()


# make beeswarm representation
my.freqs <- as.data.frame(as.matrix(t(as.matrix(my.freq.table))))
# my.freqs

apply(my.freq.table,1,mean)
#         G2         G3         G4        G5a        G5b        G5c
# 0.01096519 0.25023811 0.47068613 0.20467634 0.01575342 0.04768081

my.cols <- rep("deeppink",nrow(my.freqs))
my.cols[grep("_M_",as.character(my.freqs$Var1))] <- "dodgerblue"

pdf(paste(Sys.Date(),"10XGenomics_BM_Neutrophils_USC_Xie_SingleCellNet_predictions_BEESWARM_FREQUENCY.pdf",sep = "_"), height = 5, width = 7)
beeswarm::beeswarm(Freq ~ Var2, data = my.freqs, pc = 16, pwcol = my.cols, 
                   xlab = "Xie annotation", ylab = "Fraction of cells in 10xGenomics single cell library",
                   ylim = c(0,0.6))
abline(v = 1.5, lty = "dashed", col = "grey")
abline(v = 2.5, lty = "dashed", col = "grey")
abline(v = 3.5, lty = "dashed", col = "grey")
text(1,0.6,"preNeu (~1%)"  , cex = 0.6)
text(2,0.6,"immNeu (~25%)" , cex = 0.6)
text(3,0.6,"mNeu (~47%)"   , cex = 0.6)
text(5,0.6,"circPMN (~27%)", cex = 0.6)
dev.off()


########################################################################################################################################################
# 5. Assess marker genes for neutrophil subpopulations (G2-G5c) - Xie et al. dataset

# ntph.singlet.cl

ntph.singlet.cl@meta.data$Sex <- "NA"
ntph.singlet.cl@meta.data$Sex[grep("3m_F",ntph.singlet.cl@meta.data$Sample_ID)] <- "F"
ntph.singlet.cl@meta.data$Sex[grep("3m_M",ntph.singlet.cl@meta.data$Sample_ID)] <- "M"

ntph.singlet.cl <- SetIdent(ntph.singlet.cl, value = "SingleCellNet_Xie")

# Get Xie genes for the identified clusters
# The following requested variables were not found: Gm5483
my.xie.genes <- c("Gngt2","Gm2a","Fgl2",
                  "Ifit3","Rsad2","Isg15",
                  "Stfa2l1",  "Ccl6",       # "Gm5483",
                  "Cxcl2","Mmp8","Retnlg",
                  "Camp","Ngp","Ltf",
                  "Fcnb","Tuba1b","Chil3")

# Reorder clusters fot dotplot
# https://github.com/satijalab/seurat/issues/3500
levels(ntph.singlet.cl) <- c("G2", "G3" , "G4" , "G5a" ,  "G5b",  "G5c")

pdf(paste0(Sys.Date(),"10XGenomics_BM_Neutrophils_USC_Xie_SingleCellNet_predictions_Dotplot_markers.pdf"), width = 9, height = 5)
DotPlot(ntph.singlet.cl, features = my.xie.genes, cols = c("deeppink", "dodgerblue"), dot.scale = 8, split.by = "Sex") + RotatedAxis()
dev.off()

pdf(paste0(Sys.Date(),"10XGenomics_BM_Neutrophils_USC_Xie_SingleCellNet_predictions_Dotplot_markers_total.pdf"), width = 9, height = 3)
DotPlot(ntph.singlet.cl, features = my.xie.genes, dot.scale = 8) + RotatedAxis()
dev.off()

png(paste0(Sys.Date(),"_10XGenomics_BM_Neutrophils_USC_Xie_SingleCellNet_predictions_KNOWN_MARKER_heatmap.png"),width = 20, height = 15, units = "cm", res = 300)
DoHeatmap(ntph.singlet.cl, 
          features = my.xie.genes)  + theme(axis.text.y = element_text(size = 8))
dev.off()

pdf(paste0(Sys.Date(),"_10XGenomics_BM_Neutrophils_USC_Xie_SingleCellNet_predictions_KNOWN_MARKER_heatmap.pdf"),width = 20, height = 15)
DoHeatmap(ntph.singlet.cl, 
          features = my.xie.genes)  + theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("darkblue", "white", "red"))
dev.off()

########################################################################################################################################################
# 6. Pseudobulk analysis with muscat - sex-specific gene expression

# Convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html

# ntph.singlet.cl

ntph.sce <- as.SingleCellExperiment(ntph.singlet.cl)

# Test plots
p1 <- plotExpression(ntph.sce, features = "Ly6g", x = "SingleCellNet_Xie") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- scater::plotPCA(ntph.sce, colour_by = "SingleCellNet_Xie")
p1 + p2

pdf(paste0(Sys.Date(),"_scater_PCA_clean_ntph.pdf"))
p2
dev.off()

# Data preparation
ntph.sce.cl <- prepSCE(ntph.sce, 
                       kid = "SingleCellNet_Xie",  # subpopulation assignments
                       gid = "Sex"              ,  # group IDs (ctrl/stim)
                       sid  = "Sample_ID"       ,  # sample IDs (ctrl/stim.1234)
                       drop       = TRUE)          # drop all other colData columns


# For consistency and easy accession throughout this vignette, we will store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(ntph.sce.cl$cluster_id))
ns  <- length(sids <- levels(ntph.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(ntph.sce.cl$cluster_id, ntph.sce.cl$sample_id))
#         G2  G3  G4 G5a G5b G5c
# 3m_F_1   8 364 697 261  22  44
# 3m_F_2  14 325 525 255  27  59
# 3m_M_1  20 429 838 376  27 112
# 3m_M_2  25 377 783 343  16  78

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(ntph.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# One sheet per subpopulation
assayNames(pb)
# [1]  "G2"  "G3"  "G4"  "G5a" "G5b" "G5c" 

# Pseudobulk-level MDS plot
pb_mds <- pbMDS(pb)

pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_MDS.pdf"))
pb_mds
dev.off()

# run DS analysis
res <- pbDS(pb, method = "edgeR", min_cells = 8, verbose = TRUE)

# access results table
tbl <- res$table[[1]]

# one data.frame per cluster
names(tbl)
# "G2"  "G3"  "G4"  "G5a" "G5b" "G5c"

# print results to files
write.table(tbl[[1]][,1:7], file = paste0(Sys.Date(),"_Muscat_EdgeR_NTPH_cluster_G2_genes_statistics.txt") , sep = "\t" , row.names = F, quote=F)
write.table(tbl[[2]][,1:7], file = paste0(Sys.Date(),"_Muscat_EdgeR_NTPH_cluster_G3_genes_statistics.txt") , sep = "\t" , row.names = F, quote=F)
write.table(tbl[[3]][,1:7], file = paste0(Sys.Date(),"_Muscat_EdgeR_NTPH_cluster_G4_genes_statistics.txt") , sep = "\t" , row.names = F, quote=F)
write.table(tbl[[4]][,1:7], file = paste0(Sys.Date(),"_Muscat_EdgeR_NTPH_cluster_G5a_genes_statistics.txt"), sep = "\t" , row.names = F, quote=F)
write.table(tbl[[5]][,1:7], file = paste0(Sys.Date(),"_Muscat_EdgeR_NTPH_cluster_G5b_genes_statistics.txt"), sep = "\t" , row.names = F, quote=F)
write.table(tbl[[6]][,1:7], file = paste0(Sys.Date(),"_Muscat_EdgeR_NTPH_cluster_G5c_genes_statistics.txt"), sep = "\t" , row.names = F, quote=F)

# top-6 DS genes per cluster
pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_top5_per_cluster_heatmap.pdf"))
pbHeatmap(ntph.sce.cl, res, lfc = 0.25, top_n = 6, 
          col = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
          assay = "logcounts")
dev.off()

# top-20 DS genes for single cluster
pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_top20_for_each_cluster_heatmaps.pdf"))
pbHeatmap(ntph.sce.cl, res, lfc = 0.25, k = "G2", col = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
pbHeatmap(ntph.sce.cl, res, lfc = 0.25, k = "G3", col = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
pbHeatmap(ntph.sce.cl, res, lfc = 0.25, k = "G4", col = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
pbHeatmap(ntph.sce.cl, res, lfc = 0.25, k = "G5a", col = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
pbHeatmap(ntph.sce.cl, res, lfc = 0.25, k = "G5b", col = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
pbHeatmap(ntph.sce.cl, res, lfc = 0.25, k = "G5c", col = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

# Expression of sex-specific genes across all clusters
pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_Xist_Ddx3y_for_each_cluster_heatmaps.pdf"))
pbHeatmap(ntph.sce.cl, res, g = c("Xist","Ddx3y"), lfc = 0, fdr = 1, sort_by = "logFC", col = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

save(pb, res, file = paste0(Sys.Date(),"_SCE_with_Muscat_PsuedoBulk_ntph_object.RData"))

# Sex differences (MDS) for each neutrophil subpopulation (G2-G4, G5a-c)

ntph.G2.sce <- as.SingleCellExperiment(subset(ntph.singlet.cl, subset = SingleCellNet_Xie == "G2"))
ntph.G3.sce <- as.SingleCellExperiment(subset(ntph.singlet.cl, subset = SingleCellNet_Xie == "G3"))
ntph.G4.sce <- as.SingleCellExperiment(subset(ntph.singlet.cl, subset = SingleCellNet_Xie == "G4"))
ntph.G5a.sce <- as.SingleCellExperiment(subset(ntph.singlet.cl, subset = SingleCellNet_Xie == "G5a"))
ntph.G5b.sce <- as.SingleCellExperiment(subset(ntph.singlet.cl, subset = SingleCellNet_Xie == "G5b"))
ntph.G5c.sce <- as.SingleCellExperiment(subset(ntph.singlet.cl, subset = SingleCellNet_Xie == "G5c"))

ntph.G2.sce.cl <- prepSCE(ntph.G2.sce       , 
                          gid = "Sex"       , 
                          sid  = "Sample_ID", 
                          drop       = TRUE)    
ntph.G3.sce.cl <- prepSCE(ntph.G3.sce, 
                          kid = "Sex"       ,
                          gid = "Sex"       ,
                          sid  = "Sample_ID",
                          drop       = TRUE)
ntph.G4.sce.cl <- prepSCE(ntph.G4.sce       , 
                          kid = "Sex"       ,
                          gid = "Sex"       , 
                          sid  = "Sample_ID",
                          drop       = TRUE)
ntph.G5a.sce.cl <- prepSCE(ntph.G5a.sce, 
                           kid = "Sex"       ,
                           gid = "Sex"       ,
                           sid  = "Sample_ID",
                           drop       = TRUE)
ntph.G5b.sce.cl <- prepSCE(ntph.G5b.sce, 
                           kid = "Sex"       ,
                           gid = "Sex"       ,
                           sid  = "Sample_ID",
                           drop       = TRUE)
ntph.G5c.sce.cl <- prepSCE(ntph.G5c.sce      , 
                           kid = "Sex"       ,
                           gid = "Sex"       ,
                           sid  = "Sample_ID",
                           drop       = TRUE)

pb.G2 <- aggregateData(ntph.G2.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
pb.G3 <- aggregateData(ntph.G3.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
pb.G4 <- aggregateData(ntph.G4.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
pb.G5a <- aggregateData(ntph.G5a.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
pb.G5b <- aggregateData(ntph.G5b.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
pb.G5c <- aggregateData(ntph.G5c.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

pb.G2_mds <- pbMDS(pb.G2)
pb.G3_mds <- pbMDS(pb.G3)
pb.G4_mds <- pbMDS(pb.G4)
pb.G5a_mds <- pbMDS(pb.G5a)
pb.G5b_mds <- pbMDS(pb.G5b)
pb.G5c_mds <- pbMDS(pb.G5c)

pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_MDS_G2_by_sex.pdf"))
pb.G2_mds
dev.off()

pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_MDS_G3_by_sex.pdf"))
pb.G3_mds
dev.off()

pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_MDS_G4_by_sex.pdf"))
pb.G4_mds
dev.off()

pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_MDS_G5a_by_sex.pdf"))
pb.G5a_mds
dev.off()

pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_MDS_G5b_by_sex.pdf"))
pb.G5b_mds
dev.off()

pdf(paste0(Sys.Date(),"_10x_NTPH_Muscat_PB_MDS_G5c_by_sex.pdf"))
pb.G5c_mds
dev.off()

########################################################################################################################################################
# 7. Construct single-cell trajectory - monocle3

# ntph.singlet.cl

# Extract G2 cells (root cells)
ntph.G2.cells <- Cells(subset(ntph.singlet.cl, SingleCellNet_Xie == "G2"))
write.table(ntph.G2.cells, file = paste0(Sys.Date(),"_NTPH_G2_cell_list.txt") , sep = "\n" , row.names = F, quote=F, col.names = F)
G2 <- readLines("./2022-04-12_NTPH_G2_cell_list.txt")

# Learn graph, order cells and plot cells

ntph.cds <- as.cell_data_set(ntph.singlet.cl)
ntph.cds <- cluster_cells(cds = ntph.cds, reduction_method = "UMAP")
ntph.cds <- learn_graph(ntph.cds, use_partition = TRUE)

ntph.cds.G2 <- order_cells(ntph.cds, reduction_method = "UMAP", root_cells = G2)

pdf(paste0(Sys.Date(), "ntph_UMAP_monocle3_pseudotime_G2.pdf", sep = "_"), height = 5, width = 10)
plot_cells(cds = ntph.cds.G2,
  	 	   color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE)
dev.off()

########################################################################################################################################################
# (Not part of manuscript) 8. Run augur to figure the most sex-dimorphic subset

ntph.singlet.cl@meta.data$Sex <- "NA"
ntph.singlet.cl@meta.data$Sex[grep("3m_F",ntph.singlet.cl@meta.data$Sample_ID)] <- "F"
ntph.singlet.cl@meta.data$Sex[grep("3m_M",ntph.singlet.cl@meta.data$Sample_ID)] <- "M"

augur.ntph <-  calculate_auc(as.matrix(ntph.singlet.cl@assays$SCT@data),
                             ntph.singlet.cl@meta.data, 
                             cell_type_col = "SingleCellNet_Xie", 
                             label_col = "Sex",
                             n_threads = 1)
augur.ntph
# $AUC
# # A tibble: 6 × 2
# cell_type   auc
# <chr>     <dbl>
# 1 G5c       0.540
# 2 G5a       0.538
# 3 G5b       0.518
# 4 G3        0.516
# 5 G4        0.510
# 6 G2        0.447

save(augur.ntph, file = paste0(Sys.Date(),"_Augur_ntph_object.RData"))

pdf(paste0(Sys.Date(),"_Augur_ntph_UMAP.pdf"), height = 5, width = 6)
plot_umap(augur.ntph,ntph.singlet.cl, cell_type_col = "SingleCellNet_Xie")
dev.off()

pdf(paste0(Sys.Date(),"_Augur_ntph_UMAP_Red_Blue.pdf"), height = 5, width = 6)
plot_umap(augur.ntph,ntph.singlet.cl, cell_type_col = "SingleCellNet_Xie", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_ntph_Lollipop.pdf"))
plot_lollipop(augur.ntph)
dev.off()
#### not much here (not very stable)

########################################################################################################################################################
# 9. Differential expression analysis with Seurat

# Set identity
Idents(ntph.singlet.cl) <- "Sex"

# Get data
avg.macro <- log1p(AverageExpression(ntph.singlet.cl, verbose = FALSE)$SCT)
avg.macro$gene <- rownames(avg.macro)

pdf(paste(Sys.Date(),"BM_NTPH_Female_vs_male.pdf", sep = "_"), height = 5, width = 10)
plot(avg.macro$F,avg.macro$M, pch = 16, col = "grey", log = 'xy')
smoothScatter(avg.macro$F,avg.macro$M)
dev.off()

# MAST
sex.response <- FindMarkers(ntph.singlet.cl, ident.1 = "F", ident.2 = "M", test.use = "MAST", logfc.threshold = 0.2)
sum(sex.response$p_val_adj < 0.05, na.rm = T)

sex.response <- sex.response[!is.na(sex.response$p_val_adj),]
sex.response[sex.response$p_val_adj < 0.05,]
#                  p_val avg_log2FC pct.1 pct.2     p_val_adj
# Xist      0.000000e+00  0.6197908 0.354 0.000  0.000000e+00
# Eif2s3y  2.762694e-201 -0.3555775 0.000 0.217 3.093941e-197
# Ddx3y    7.740711e-162 -0.2965546 0.000 0.181 8.668823e-158
# G0s2     8.383518e-112  0.5786645 0.626 0.390 9.388702e-108
# Uty      7.943555e-108 -0.2376009 0.000 0.120 8.895987e-104
# Hdc       4.546197e-52  0.4294456 0.878 0.768  5.091286e-48
# Rap1gap2  1.380331e-50  0.3579984 0.850 0.728  1.545832e-46
# Lrg1      2.557613e-47  0.4134569 0.593 0.428  2.864270e-43
# Pim1      5.186573e-47  0.3700524 0.529 0.368  5.808443e-43
# Dgat1     1.875326e-46  0.3655076 0.639 0.509  2.100177e-42

write.table(sex.response[sex.response$p_val_adj < 0.05,], file = paste0(Sys.Date(),"_female_vs_male_ntph_MAST.txt"), sep = "\t", quote = F)

########################################################################################################################################################
# 10. get stats for data descriptor

table(ntph.singlet.cl@meta.data$Sample_ID)
# 3m_F_1 3m_F_2 3m_M_1 3m_M_2 
# 1396   1205   1802   1622 

# Median UMI counts of QC neutrophils
aggregate(ntph.singlet.cl@meta.data$nCount_RNA, by = list(ntph.singlet.cl@meta.data$Sample_ID), FUN = "median")
#   Group.1      x
# 1  3m_F_1 5551.5
# 2  3m_F_2 5980.0
# 3  3m_M_1 4715.5
# 4  3m_M_2 4765.0

# Median detected gene counts of QC neutrophils
aggregate(ntph.singlet.cl@meta.data$nFeature_RNA, by = list(ntph.singlet.cl@meta.data$Sample_ID), FUN = "median")
# Group.1      x
# 1  3m_F_1 1662.0
# 2  3m_F_2 1737.0
# 3  3m_M_1 1540.5
# 4  3m_M_2 1540.0

# Median HTO counts of QC neutrophils
aggregate(ntph.singlet.cl@meta.data$nCount_HTO, by = list(ntph.singlet.cl@meta.data$Sample_ID), FUN = "median")
# Group.1   x
# 1  3m_F_1 331
# 2  3m_F_2 208
# 3  3m_M_1 273
# 4  3m_M_2 260

# Mean mt read percentage
aggregate(ntph.singlet.cl@meta.data$percent.mito, by = list(ntph.singlet.cl@meta.data$Sample_ID), FUN = "mean")
# Group.1         x
# 1  3m_F_1 0.3765691
# 2  3m_F_2 0.3621518
# 3  3m_M_1 0.3854459
# 4  3m_M_2 0.3684778

summary(ntph.singlet.cl@meta.data$percent.mito,)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1552  0.2865  0.3742  0.4768 16.7742 

########################################################################################################################################################
sink(paste0(Sys.Date(),"_NTPH_scRNAseq_analysis_sessionInfo.txt"))
sessionInfo()
sink()
