# Hemp-Genotype-variation-in-the-microbiome
Paper in review: Ahmad W., Coffman L., Flavier A., Crawford K., Weerasooriya A., Khan A. (2023) Microbiome Diversity of Industrial CBD and Fiber Hemp Genotypes. Frontiers in Plant Science.

# Data 
The data folder contains the following files: \
• 16S/filtered-phyllosphere-16S.qza : A QIIME2 artifact with the 16S ASV abundances for phyllosphere. \
• 16S/filtered-table-root16S.qza : A QIIME2 artifact with the 16S ASV abundances for root. \
• 16S/filtered-table-soil16S.qza : A QIIME2 artifact with the 16S ASV abundances for soil. \
• 16S/taxonomy16S.qza: A QIIME2 artifact the taxonomic assignment for 16S ASV. \
• 16S/rooted-tree16S.qza: A QIIME2 artifact containing the phylogenetic tree for 16S ASVs. \
• 16S/shannon_vector_F001.qza: A QIIME2 artifact containing the shannon diversity for the soil 16S ASVs. \
• 16S/shannon_Rhizosphere16S.qza: A QIIME2 artifact containing the shannon diversity for the rhizosphere 16S ASVs. \
• 16S/shannon_phyllosphere16S.qza: A QIIME2 artifact containing the shannon diversity for the phyllosphere 16S ASVs. \
• 16S/Hempmetadata-upd-16S.tsv: A table with the 16S sample information. \

• ITS/table-phyllosphere-ITS.qza : A QIIME2 artifact with the ITS ASV abundances for phyllosphere. \
• ITS/table-rhizosphere-ITS.qza : A QIIME2 artifact with the ITS ASV abundances for rhizosphere. \
• ITS/table-trim_ITS.qza : A QIIME2 artifact with the trimmed ITS ASV abundances. \
• ITS/taxonomy_ITS.qza: A QIIME2 artifact the taxonomic assignment for ITS ASV. \
• ITS/rooted-tree-ITS.qza: A QIIME2 artifact containing the phylogenetic tree for ITS ASVs. \
• ITS/shannon_phyllosphere_ITS.qza: A QIIME2 artifact containing the shannon diversity for the phyllosphere ITS ASVs. \
• ITS/shannon_rhizosphere_ITS.qza: A QIIME2 artifact containing the shannon diversity for the rhizosphere ITS ASVs. \
• ITS/Hempmetadata_ITS.tsv: A table with the ITS sample information. \

# Code 
The src folder contains the following R scripts for figure creation: \
• phyloseq.R: Moved QIIME2 .qza files into a phyloseq object \
• core.microbiome.R: Using phylseq object found core microbiome \
• venn-diagram.R: Using asv table created venn diagrams \
• Qiime2Rstudio.R: Create Shannon diversity plots using qiime2R package and tidyverse \
• ANCOM.R: Using phylseq object and running ANCOM-BC2 differential abundance analysis \

# Requirements 
To run the code, you will need R-4.2.1 and the following R packages: \
• tidyverse \
• qiime2R \
• phyloseq \
• microbiome \
• dplyr \
• vegan \
• HMP \
• ggplot2 \
• ggvenn \
• microViz \
• eulerr \
• dendextend \
• microbiomeSeq \
• microbial \
• plotly \
• picante \
• ANCOMBC \
• kableExtra \
• breakaway \
• apeglm \
• knitr \
• ggpubbr \
• MicrobeR \
• microbiomeutilities \
• RColorBrewer \
