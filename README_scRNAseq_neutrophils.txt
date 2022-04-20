###############################
##########  README   ########## 

**Single cell RNA-seq of female and male primary mouse bone marrow neutrophils**

Minhoo Kim, Ryan J. Lu and Bérénice A. Benayoun


This code was used to process, demultiplex, annotate and quality control 10xGenomics 
single cell RNA-seq data from young female and male murine bone marrow neutrophils.

**********************************
*****   Standalone software  *****
R 4.1.2
cellranger 6.0.2
CITE-seq-Count 1.4.5
**********************************

**********************************
*****   R Package versions   *****
Seurat 4.1.0
SingleR 1.8.1
singleCellNet 0.1.0 
monocle3 1.0.0
muscat 1.5.2
edgeR 3.36.0
**********************************


     - 1_Cell_Ranger   : contains the code used to run cellranger from fastq files (RNA quantification)
     - 2_CITEseq_tools : contains the code used to run CITEseq tools from fastq files (HTO quantification)
     - 3_R_analysis    : contains the code to process and annotate the data 


The corresponding raw sequencing data has been deposited under PRJNA796634.
