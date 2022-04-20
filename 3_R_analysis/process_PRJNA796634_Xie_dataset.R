setwd('/Volumes/BB_Home/Neutrophils/Public_data/Neutrophil_maturation_data/scRNA-seq_dataset')
options(stringsAsFactors = F)

library('Seurat')

# Process PRJNA796634 (Xie et al. dataset) for neutrophil annotation

########################################################################################################################################################
my.meta.data <- read.csv("GSE137539_Mouse_WT_blank_meta.txt" , "\t", header = T)
my.barcodes  <- read.csv("GSE137539_Cell_index_Barcode.txt"  , "\t", header = T)
my.wt_bm1    <- read.csv("GSM4081545_wt_ctl_bm1.txt"         , " ", header = T)
my.wt_bm2    <- read.csv("GSM4081546_wt_ctl_bm2.txt"         , " ", header = T)
my.wt_pb2    <- read.csv("GSM4081547_wt_ctl_pb2.txt"         , " ", header = T)
my.wt_sp2    <- read.csv("GSM4081548_wt_ctl_sp2.txt"         , " ", header = T)

colnames(my.wt_bm1)    <- paste0("wt.ctl.bm1_", colnames(my.wt_bm1) )
colnames(my.wt_bm2)    <- paste0("wt.ctl.bm2_", colnames(my.wt_bm2) )
colnames(my.wt_pb2)    <- paste0("wt.ctl.pb2_", colnames(my.wt_pb2) )
colnames(my.wt_sp2)    <- paste0("wt.ctl.sp2_", colnames(my.wt_sp2) )

my.wt_bm1.Seurat <- CreateSeuratObject(my.wt_bm1, meta.data = my.meta.data[my.meta.data$orig.ident %in% "wt.ctl.bm1",] )
my.wt_bm2.Seurat <- CreateSeuratObject(my.wt_bm2, meta.data = my.meta.data[my.meta.data$orig.ident %in% "wt.ctl.bm2",] )
my.wt_pb2.Seurat <- CreateSeuratObject(my.wt_pb2, meta.data = my.meta.data[my.meta.data$orig.ident %in% "wt.ctl.pb2",] )
my.wt_sp2.Seurat <- CreateSeuratObject(my.wt_sp2, meta.data = my.meta.data[my.meta.data$orig.ident %in% "wt.ctl.sp2",] )

Xie.combined <- merge(x = my.wt_bm1.Seurat, y = list(my.wt_bm2.Seurat,my.wt_pb2.Seurat,my.wt_sp2.Seurat), project = "Xie_scNeutrophils")

my.non.na <- Xie.combined@meta.data[rownames(Xie.combined@meta.data[!is.na(Xie.combined@meta.data$cell_type),]),]
my.ntph   <- my.non.na[!is.na(my.non.na$cluster),]

Xie.combined.filt <- subset(Xie.combined, cells = rownames(my.ntph))
# An object of class Seurat 
# 27998 features across 12285 samples within 1 assay 
# Active assay: RNA (27998 features, 0 variable features)


save(Xie.combined.filt, file = paste0(Sys.Date(),"_NTPH_Xie_paper_Seurat_obj.RData"))	# 2020-08-31_NTPH_Xie_paper_Seurat_obj.RData

########################################################################################################################################################

sink(paste0(Sys.Date(),"_PRJNA796634_Xie_dataset_processing_sessionInfo.txt"))
sessionInfo()
sink()
