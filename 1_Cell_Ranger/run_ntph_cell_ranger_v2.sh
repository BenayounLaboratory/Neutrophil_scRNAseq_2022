#### cellranger 6.0.2

cellranger count --id 3m_F_NTPH_v2 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /home/benayoun/Projects/10xGenomics/2021-05-18_scNTPH/210517_NB551671_0100_AHCTMHBGXJ_BASECALLS/ \
                 --sample 3m_F_NTPH \
                 --expect-cells 5000 \
                 --localmem 64 \
                 --localcores 12 \
                 --include-introns
                 
                 
cellranger count --id 3m_M_NTPH_v2 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /home/benayoun/Projects/10xGenomics/2021-05-18_scNTPH/210517_NB551671_0100_AHCTMHBGXJ_BASECALLS/ \
                 --sample 3m_M_NTPH \
                 --expect-cells 5000 \
                 --localmem 64 \
                 --localcores 12 \
                 --include-introns
