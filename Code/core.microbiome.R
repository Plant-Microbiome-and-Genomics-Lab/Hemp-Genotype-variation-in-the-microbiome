library(microbiome)

ps.prune <- prune_taxa(taxa_sums(ps) > 0, ps)
print(ps.prune)

pseq.rel <- microbiome::transform(ps.prune, "compositional")

core.taxa.standard <- core_members(pseq.rel, detection = 0.0001, prevalence = 50/100)
core.taxa.standard

core.taxa <- taxa(pseq.core)
class(core.taxa)

pseq.core <- core(pseq.rel, detection = 0.0001, prevalence = .5)

# get the taxonomy data
tax.mat <- tax_table(pseq.core)
tax.df <- as.data.frame(tax.mat)

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
write.csv(core.taxa.class, "core-Complete_ITS.csv")

