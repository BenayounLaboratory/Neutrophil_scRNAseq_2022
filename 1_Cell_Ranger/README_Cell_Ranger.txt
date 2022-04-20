##########################################
###############   README   ###############
##########################################


*** 10xGenomics library quantification ***

- Quantify 10x chromium 3' RNA-seq libraries using Cell Ranger (v6.0.2). 
  Use "intron" mode as per 10xGenomics recommendations since neutrophils have very high
  intron-retention levels.
  (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/tutorials/neutrophils?src=social&lss=wechat&cnm=soc-we-ra_g-program&cid=7011P000000oOU6)
  
	% run_ntph_cell_ranger_v2.sh: bash script for the Cell Ranger Count step
	% Output: folder containing the compressed output folders from quantification