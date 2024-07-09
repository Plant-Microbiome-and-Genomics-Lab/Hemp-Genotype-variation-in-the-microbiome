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
