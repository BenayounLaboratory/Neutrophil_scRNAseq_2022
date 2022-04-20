##########################################
###############   README   ###############
##########################################


*** R analysis ***

- Retrieve/reconstitute the annotated single cell RNA-seq object from Xie et al (GSE137539)
	 Using version: R 3.6.3 and Seurat 3.2.2

	% Files retrieved from GEO:
				GSE137539_Mouse_WT_blank_meta.txt
				GSE137539_Cell_index_Barcode.txt
				GSM4081545_wt_ctl_bm1.txt        
				GSM4081546_wt_ctl_bm2.txt        
				GSM4081547_wt_ctl_pb2.txt        
				GSM4081548_wt_ctl_sp2.txt        
	% process_PRJNA796634_Xie_dataset.R: R script to generate Seurat object reference from Xie et al data
	% 2020-08-31_NTPH_Xie_paper_Seurat_obj.RData: Seurat object used for annotation
  

- Process, annotate and Quality check the scRNA-seq dataset
  
	 Using version: R 4.1.2
	                Seurat 4.1.0
	                SingleR 1.8.1
	                singleCellNet 0.1.0 
	                monocle3 1.0.0
	                muscat 1.5.2
	                edgeR 3.36.0
  
	% process_10x_BM_Ntph_v7.R: script to analyze the dataset in R
