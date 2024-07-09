library(vegan)
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(microViz)
library(eulerr)
library(microbiomeSeq)
library(microbial)


#Load data into phyloseq
ps <- qza_to_phyloseq("table-phyllosphere16S.qza", "rooted-tree16S.qza", 
                              "taxonomy16S.qza", "Hempmetadata-upd-16S.tsv")

#####################Agglomerating and subsetting taxa #################
########################################################################
ps_Phylum <- tax_glom(ps, "Phylum")
otu_Phylum <- otu_table(ps_Phylum)
otu_Phylum <- t(otu_Phylum)
df_Phylum <- as.data.frame(otu_Phylum)
head(tax_table(ps_Phylum))
identical(colnames(df_Phylum), taxa_names(ps_Phylum))
colnames(df_Phylum) <- as.data.frame(tax_table(ps_Phylum))$Phylum
df_Phylum <- t(df_Phylum)
write.csv(otu_table(ps_Phylum), "df_phylum_CBDvsFiber-Complete_ITS.csv")

ps_Genus <- tax_glom(ps, "Genus")
otu_Genus <- otu_table(ps_Genus)
otu_Genus <- t(otu_Genus)
df_Genus <- as.data.frame(otu_Genus)
head(tax_table(ps_Genus))
identical(colnames(df_Genus), taxa_names(ps_Genus))
colnames(df_Genus) <- as.data.frame(tax_table(ps_Genus))$Genus
df_Genus <- t(df_Genus)
print(df_Genus)
write.csv(df_Genus, "df_genus_CBDvsFiber-Complete_ITS.csv")

#Agglomerating and subsetting taxa- !Relative  abundance!
physeq_rel <- microbiome::transform(ps, "compositional")

ps_rel_Phylum <- tax_glom(physeq_rel, "Phylum")
otu_rel_Phylum <- otu_table(ps_rel_Phylum)
otu_rel_Phylum <- t(otu_rel_Phylum)
df_rel_Phylum <- as.data.frame(otu_rel_Phylum)
head(tax_table(ps_rel_Phylum))
identical(colnames(df_rel_Phylum), taxa_names(ps_Phylum))
colnames(df_rel_Phylum) <- as.data.frame(tax_table(ps_Phylum))$Phylum
df_rel_Phylum <- t(df_rel_Phylum)
write.csv(otu_table(ps_Phylum), "df_phylum_CBDvsFiber-Complete_ITS.csv")


ps_rel_Genus <- tax_glom(physeq_rel, "Genus")
otu_rel_Genus <- otu_table(ps_rel_Genus)
otu_rel_Genus <- t(otu_rel_Genus)
df_rel_Genus <- as.data.frame(otu_rel_Genus)
head(tax_table(ps_rel_Genus))
identical(colnames(df_rel_Genus), taxa_names(ps_Genus))
colnames(df_rel_Genus) <- as.data.frame(tax_table(ps_Genus))$Genus
df_rel_Genus <- t(df_rel_Genus)
#print(df_Genus)
write.csv(df_rel_Genus, "rel_abd_genus_CBDvsFiber-Complete_ITS.csv")


###############################################
###############Bray Curtis#####################
################################################
GP.NMDS = ordinate(ps, "PCoA", "bray")
ptitle = "PCoA of Bray-Curtis"
p = plot_ordination(ps, GP.NMDS, type = 'samples', color = "Replicates",
                    shape='Use', title = ptitle, point_size = 10)
print(p)



ggsave("bray-curtis-PCoA-CBDvsFiber-soil_ITS.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave("bray-curtis-PCoA-CBDvsFiber-soil_ITS.png", height=4, width=5, device="png", dpi=600)

#################################################
#get names of most abd phyla and use for subsetting
top.phyla = sort(tapply(taxa_sums(ps), tax_table(phyloseq)[,"Phylum"], sum), TRUE)
top.phyla = top.phyla[1:5] #Prune to just the most abd 5 phyla
physeq_top_5_Phyla = subset_taxa(ps, Phylum%in% names(top.phyla))
get_taxa_unique(physeq_top_5_Phyla, "Phylum")
physeq_top_5_Phyla = prune_taxa(names(sort(taxa_sums(physeq_top_5_Phyla), TRUE)[1:100]), physeq_top_5_Phyla)#Prune futher to top 200 most abd taxa
p2= plot_ordination(physeq_top_50, ordinate(physeq_top_50, "CCA"), type = "samples", color = "Genotype")
p2 + geom_point(size = 5) + geom_polygon(aes(fill = treatment))

#Plot heatmap
plot_heatmap(ps, "NMDS", "bray", "Genotype", "Family")
ggsave("p.pdf", height=5, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave("phyloseq-CBDvsFiber_ITS.png", height=5, width=5, device="png", dpi=600)

#Microbiome 
physeq_rel <- microbiome::transform(ps, "compositional")

plot_net(physeq_rel, maxdist = 0.8, color = "Genotype")

ig <- make_network(physeq_rel, dist.fun="bray", max.dist=0.8)
plot_network(ig, physeq_rel, color="Genotype", line_weight=0.4, label=NULL)

##########################################################
#######################Taxanomic Barplot##########################
########################phylum#########################

plotbar(ps, level = "Phylum", color = NULL, group = "Replicates", top = 50,
        fontsize.x = 15, fontsize.y = 12)


ggsave("barplot-Phylum-CBDvsFiber-Complete-ITS.pdf", height=4, width=10, device="pdf") # save a PDF 4 inches by 8 inches
ggsave("barplot-Phylum-CBDvsFiber-Complete-ITS.png", height=4, width=10, device="png", dpi=600)

########################Family#########################

plotbar(ps, level = "Family", color = NULL, group = "Replicates", top = 15, 
        fontsize.x = 15, fontsize.y = 12)

ggsave("barplot-family-CBDvsFiber-Complete-ITS.pdf", height=4, width=10, device="pdf") # save a PDF 4 inches by 8 inches
ggsave("barplot-family-CBDvsFiber-Complete-ITS.png", height=4, width=10, device="png", dpi=600)

########################Genus#########################


plotbar(ps, level = "Genus", color = NULL, group = "Replicates", top = 15,
        fontsize.x = 15, fontsize.y = 12)

ggsave("barplot-Genus-CBDvsFiber-Complete-ITS.pdf", height=4, width=10, device="pdf") # save a PDF 4 inches by 8 inches
ggsave("barplot-Genus-CBDvsFiber-Complete-ITS.png", height=4, width=10, device="png", dpi=600)



