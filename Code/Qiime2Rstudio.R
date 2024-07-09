library(qiime2R)
library(tidyverse)

#################Shannon_Bar Plot#################
##################################################

metadata <- read_q2metadata("Hempmetadata_ITS_CBDvsHemp.tsv")
shannon <- read_qza("shannon-CDvsHemp-phyllo-ITS.qza")

shannon <- shannon$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(shannon)


metadata %>%
  dplyr::filter(!is.na(shannon_entropy)) %>%
  ggplot(aes(x=`Genotype`, y=shannon_entropy, fill=`Genotype`)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(0,8)) + # adjust y-axis
  facet_grid(~`Organ`) + # create a panel for each body site
  ggtitle("rhizoITS") +
  xlab("Genotype") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(values=c("cornflowerblue","indianred", "green", "yellow", "purple", "orange", "dark red", "light blue")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed
ggsave("Shannon_rhizosphere_ITS.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches
ggsave("Shannon_rhizosphere_ITS.png", height=4, width=4, device="png", dpi=600)


#################Shannon_Bar Plot#################
#################Use comparison#################
##################################################

metadata <- read_q2metadata("Hempmetadata-upd-16S.tsv")
shannon <- read_qza("shannon-CDvsHemp-soil-16s.qza")

shannon <- shannon$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(shannon)


metadata %>%
  dplyr::filter(!is.na(shannon_entropy)) %>%
  ggplot(aes(x=`Use`, y=shannon_entropy, fill=`Use`)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(0,15)) + # adjust y-axis
  #facet_grid(~`Replicates`) + # create a panel for each body site
  ggtitle("soil-cbd vs hemp16S") +
  xlab("Use") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(values=c("cornflowerblue","indianred", "green", "yellow", "purple", "orange", "dark red", "light blue")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed
  ggsave("Shannon_CBDvsHemp-soil-16s.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches
  ggsave("Shannon_CBDvsHemp-soil-16s.png", height=4, width=4, device="png", dpi=600)

################Shannon_PCoA Plot################
###############Bray Curtis################
#################################################

braycurtis <-read_qza("bray-curtis-CBDvsHemp-phyllo-ITS.qza")

braycurtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Genotype`, shape=`Organ`, size=shannon_entropy)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  ggtitle("BCRphylloITS") +
  theme_q2r() +
  scale_shape_manual(values=c(16,1,2,3), name="Genotype") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Genotype")
ggsave("bray-curtis-PCoA-phyllosphere-ITS.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave("bray-curtis-PCoA-phyllosphere-ITS.png", height=4, width=5, device="png", dpi=300)

################# Plot Genus Level#################
##########################################################

metadata_q2r<-read_q2metadata("Hempmetadata_ITS_CBDvsHemp.tsv")
SVs_q2r<-read_qza("table-rhizosphere_ITS.qza")$data
taxonomy_q2r<-read_qza("taxonomy_ITS.qza")$data %>% parse_taxonomy()

taxasums_Genus<-summarize_taxa(SVs_q2r, taxonomy_q2r)$Genus

taxa_heatmap(taxasums_Genus, metadata_q2r, "Genotype")

ggsave("heatmap-genus-rhizosphere-ITS.pdf", height=8, width=16, device="pdf") # save a PDF 4 inches by 8 inches
ggsave("heatmap-genus-rhizosphere-ITS.png", height=8, width=14, device="png", dpi=600)

taxa_barplot(taxasums_Genus, metadata_q2r, "Replicates")

ggsave("barplot-genus-rhizosphere-ITS.pdf", height=4, width=10, device="pdf") # save a PDF 4 inches by 8 inches
ggsave("barplot-genus-rhizosphere-ITS.png", height=4, width=10, device="png", dpi=600)

################# Plot Phylum Level#################
###########################################################

taxasums_Phylum <-summarize_taxa(SVs_q2r, taxonomy_q2r)$Phylum

taxa_heatmap(taxasums_Phylum, metadata_q2r, "Genotype")

ggsave("heatmap-phylum-rhizosphere-ITS.pdf", height=6, width=10, device="pdf") # save a PDF 4 inches by 8 inches
ggsave("heatmap-phylum-rhizosphere-ITS.png", height=6, width=10, device="png", dpi=600)

taxa_barplot(taxasums_Phylum, metadata_q2r, "Organ")

ggsave("barplot-phylum-rhizosphere-ITS.pdf", height=4, width=5, device="pdf") # save a PDF 4 inches by 8 inches
ggsave("barplot-phylum-rhizosphere-ITS.png", height=4, width=5, device="png", dpi=600)

################# Family Level#################
#####################################################

taxasums_family<-summarize_taxa(SVs_q2r, taxonomy_q2r)$Family

taxa_heatmap(taxasums_family, metadata_q2r, "Genotype")

ggsave("heatmap-family-rhizosphere-ITS.pdf", height=8, width=14, device="pdf") # save a PDF 4 inches by 8 inches
ggsave("heatmap-family-rhizosphere-ITS.png", height=8, width=14, device="png", dpi=600)

taxa_barplot(taxasums_family, metadata_q2r, "Genotype")

ggsave("barplot-family-rhizosphere-ITS.pdf", height=4, width=10, device="pdf") # save a PDF 4 inches by 8 inches
ggsave("barplot-family-rhizosphere-ITS.png", height=4, width=10, device="png", dpi=600)

