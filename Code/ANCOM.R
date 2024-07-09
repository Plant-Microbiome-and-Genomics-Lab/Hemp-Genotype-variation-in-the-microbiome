#Differential abundance analysis

library(tidyverse)
library(phyloseq)
library(DESeq2)
library(microbiome)
library(vegan)
library(picante)
library(ALDEx2)
library(metagenomeSeq)
library(HMP)
library(dendextend)
library(selbal)
library(rms)
library(breakaway)
library(ANCOMBC)
library(kableExtra)

#sorts samples based on total read counts
sort(phyloseq::sample_sums(ps))

#Visualize relative abundance
#get count of phyla
table(phyloseq::tax_table(ps)[, "Phylum"])

#Convert to relative abundance
ps_rel_abund <- phyloseq::transform_sample_counts(ps, function(x){x/sum(x)})
phyloseq::otu_table(ps)[1:5, 1:5]
phyloseq::otu_table(ps_rel_abund)[1:5, 1:5]

#Agglomerate to phylum level and rename
ps_phylum <- phyloseq::tax_glom(ps, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

#Metl and plot
phyloseq::psmelt(ps_phylum) %>% 
  ggplot(data = ., aes(x = Replicates, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = OTU), height = 0, width = 0.2) + 
  labs(x = "", y = "Abundances/n") +
  facet_wrap(~ OTU, scales = "free")
ggsave("Abd-meltplot-phylum_Complete_ITS.pdf", height=8, width=20, device="pdf") # save a PDF 3 inches by 4 inches
ggsave("Abd-meltplot-phylum_phyllosphere_ITS.png", height=8, width=10, device="png", dpi=600)


####################### ANCOM-BC2 ##########################

ancom_data <- ancombc2(ps, tax_level = "Genus", group = "Genotype", p_adj_method = "fdr", 
                       fix_formula = "Genotype", pairwise = TRUE, alpha = 0.05, global = TRUE,
                       dunnet = TRUE)

ancom_res <- ancom_data$res
write.csv(ancom_res,"ancom_16S.csv")
