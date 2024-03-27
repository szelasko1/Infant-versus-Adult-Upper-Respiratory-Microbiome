library("phyloseq")
library("ggplot2")
library("readr")
library("decontam")
library("gridExtra")
library("vegan")
library("ggvegan")
library("dplyr")
library("reshape2")
library("microbiome")
library("tidyverse")
library("ggpubr")
library("ggthemes")
library("microViz")

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#import data from csv to phyloseq object
otumat_nasaloralHMP_all <- read.csv("MG0002_NasalOral_HMP_counts.csv", row.names=1)
otumat_nasaloralHMP_all <- data.matrix(otumat_nasaloralHMP_all)

taxmat_nasaloralHMP_all <- read.csv("HMP_NasalOral_Taxonomy.csv", row.names=8, na = "NA")
taxmat_nasaloralHMP_all <- as.matrix(taxmat_nasaloralHMP_all)
colnames(taxmat_nasaloralHMP_all) <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

sampledata_nasaloralHMP_all <- read.csv("MG0002_SubmissionPlate.csv")
rownames(sampledata_nasaloralHMP_all) <- colnames(otumat_nasaloralHMP_all)
SAM_nasaloralHMP_all <- sample_data(sampledata_nasaloralHMP_all, errorIfNULL = TRUE)

OTU_nasaloralHMP_all = otu_table(otumat_nasaloralHMP_all, taxa_are_rows = TRUE)
TAX_nasaloralHMP_all = tax_table(taxmat_nasaloralHMP_all)
physeq_nasaloralHMP_all = phyloseq(OTU_nasaloralHMP_all, TAX_nasaloralHMP_all, SAM_nasaloralHMP_all, package="decontam")

#identify contaminants by prevalence
sample_data(physeq_nasaloralHMP_all)$is.neg <- sample_data(physeq_nasaloralHMP_all)$Sample_or_Control == "Negative Control"
contamdf.prev_nasaloralHMP_all <- isContaminant(physeq_nasaloralHMP_all, method="prevalence", neg="is.neg", threshold = 0.1)
table(contamdf.prev_nasaloralHMP_all$contaminant)

ps.noncontam.nasaloralHMP_all <- prune_taxa(!contamdf.prev_nasaloralHMP_all$contaminant, physeq_nasaloralHMP_all)

nc.nasaloralHMP_all.clean = subset_samples(ps.noncontam.nasaloralHMP_all, Sample_or_Control == "Sample")

#prevalence table
prevelancedf_nasaloralHMP_all = apply(X = otu_table(nc.nasaloralHMP_all.clean),
                                      MARGIN = 1,
                                      FUN = function(x){sum(x > 0)})

prevelancedf_nasaloralHMP_all = data.frame(Prevalence = prevelancedf_nasaloralHMP_all,
                                           TotalAbundance = taxa_sums(nc.nasaloralHMP_all.clean),
                                           tax_table(nc.nasaloralHMP_all.clean))
prevelancedf_nasaloralHMP_all[1:10,]

#whole phylum filtering
nc.nasaloralHMP_all.clean.1 <- subset_taxa(nc.nasaloralHMP_all.clean, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "<NA>", "NA"))

plyr::ddply(prevelancedf_nasaloralHMP_all, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

phyla2Filter = c("NA", "<NA>", "", "Artverviricota", "Chordata", "Cossaviricota", "Cressdnaviricota",
                 "Dividoviricota", "Duplornaviricota", "Hofneiviricota", "Kitrinoviricota", "Negarnaviricota",
                 "Phixviricota", "Pisuviricota", "Preplasmiviricota", "Saleviricota", "Taleaviricota")

nc.nasaloralHMP_all.clean.p = subset_taxa(nc.nasaloralHMP_all.clean.1, !Phylum %in% phyla2Filter)

#agglomerate at Species level, remove all without Genus-level assignment
length(get_taxa_unique(nc.nasaloralHMP_all.clean.p, taxonomic.rank = "Species"))
nc.nasaloralHMP_all.clean.2 = tax_glom(nc.nasaloralHMP_all.clean.p, "Species", NArm = TRUE)

#rarify
standf = function(x) round(1E5 * (x / sum(x)))
nc.nasaloralHMP_all.clean.3.depth = transform_sample_counts(nc.nasaloralHMP_all.clean.2, standf)

#subset samples
nc.nasaloralHMP_all.clean.3.depth_adult = subset_samples(nc.nasaloralHMP_all.clean.3.depth, Age == "Adult")
nc.nasaloralHMP_all.clean.3.depth_infant = subset_samples(nc.nasaloralHMP_all.clean.3.depth, Age == "Infant")
nc.nasaloralHMP_all.clean.3.depth_nasal = subset_samples(nc.nasaloralHMP_all.clean.3.depth, Site == "Nasal")
nc.nasaloralHMP_all.clean.3.depth_oral = subset_samples(nc.nasaloralHMP_all.clean.3.depth, Site == "Oral")

#bacteria only
bac.nc.nasaloralHMP_all.clean.3 = subset_taxa(nc.nasaloralHMP_all.clean.3.depth, Superkingdom == "Bacteria")
bac.nasaloralHMP_all.clean.3.depth_nasal = subset_samples(bac.nc.nasaloralHMP_all.clean.3, Site == "Nasal")
bac.nasaloralHMP_all.clean.3.depth_oral = subset_samples(bac.nc.nasaloralHMP_all.clean.3, Site == "Oral")
bac.nasaloralHMP_all.clean.3.depth_adult = subset_samples(bac.nc.nasaloralHMP_all.clean.3, Age == "Adult")
bac.nasaloralHMP_all.clean.3.depth_infant = subset_samples(bac.nc.nasaloralHMP_all.clean.3, Age == "Infant")

#ordination
dist = phyloseq::distance(nc.nasaloralHMP_all.clean.3.depth_oral, method="bray")
ordination = ordinate(nc.nasaloralHMP_all.clean.3.depth_oral, method="NMDS", distance= dist)

#beta-diversity statistics
dist = phyloseq::distance(nc.nasaloralHMP_all.clean.3.depth_nasal, method="bray")
metadata <- data.frame(sample_data(nc.nasaloralHMP_all.clean.3.depth_nasal))

#PERMANOVA
adonis2(dist ~ Age, data = metadata, permutations = 999)

metadata <- data.frame(sample_data(nc.nasaloralHMP_all.clean.3.depth))
cbn <- combn(x=unique(metadata$Source_Name), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(nc.nasaloralHMP_all.clean.3.depth, Source_Name %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method="jaccard", binary = TRUE) ~ Source_Name, 
                                data = metadata_sub, permutations = 999)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#betadisper
beta <- betadisper(dist, metadata$Age, type = "median")
permutest(beta)

#ANOSIM
anosim(dist, metadata$Site, permutations = 999)

metadata <- data.frame(sample_data(nc.nasaloralHMP_all.clean.3.depth))
cbn <- combn(x=unique(metadata$Source_Name), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(nc.nasaloralHMP_all.clean.3.depth, Source_Name %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- anosim(phyloseq::distance(ps.subs, method="bray"), 
                               metadata_sub$Source_Name, permutations = 999)
  p <- c(p, permanova_pairwise$signif[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#ordination plot
plot_ord <- plot_ordination(nc.nasaloralHMP_all.clean.3.depth, ordination, type = "samples",
                            color="Source_Name") + 
  stat_ellipse(type = "t", geom = "polygon", show.legend = FALSE, linetype = 0,
               aes(fill = Source_Name, alpha = 0.001)) +
  geom_point(size = 6, alpha = 0.75, aes(shape = Age))+
  scale_shape_manual(values = c(17, 16))+
  facet_grid(~factor(Age, levels=c("Infant", "Adult")))

#alpha diversity
alpha_div <- estimate_richness(physeq = na.omit(nc.nasaloralHMP_all.clean.3.depth))
metadata <- sample_data(object = nc.nasaloralHMP_all.clean.3.depth) %>% 
  data.frame(.) 
stopifnot( all( rownames(metadata) == rownames(alpha_div) ) )
alpha_div_metadata <- cbind(alpha_div, metadata)

alpha_div_metadata_plot <- select(alpha_div_metadata, Observed, Shannon, InvSimpson, Fisher, Age, Site)
nasaloralHMP_all_alpha <- melt(alpha_div_metadata_plot) 

ggplot(nasaloralHMP_all_alpha, aes(x = Site, y = value))+
  geom_jitter(alpha=0.3, size = 4, aes(color = Site), width = 0.25, height = 0.25)+
  geom_boxplot(alpha=0.5, aes(fill = Site), outlier.alpha = 0)+
  facet_grid(variable~Age, scales = "free_y"))






top20G.names.bac.n = sort(tapply(taxa_sums(bac.nasaloralHMP_all.clean.3.depth_nasal), 
                                 tax_table((bac.nasaloralHMP_all.clean.3.depth_nasal))[, "Genus"], sum), TRUE)[1:20]
top20G.bac.n = subset_taxa(bac.nasaloralHMP_all.clean.3.depth_nasal, Genus %in% names(top20G.names.bac.n))

transform <- microbiome::transform
pseq_g <- transform(top20G.bac.n, "compositional") %>%
  aggregate_taxa(level = "Genus") %>%
  subset_samples(Age == "Adult")

pst <- pseq_g %>% ps_arrange(desc(Cutibacterium), .target = "otu_table")

plot_composition(pst, otu.sort = "abundance")+
  labs(title = "Adult")+
  xlab("")+
  ylab("Abundance (%)")+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c("0", "25", "50", "75", "100"))+
  scale_fill_manual(values = nasal_adult)+
  theme_Publication()+
  theme(legend.text = element_text(size = 24, face = "italic"),
        legend.title = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key.size = unit(1, units = "cm"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 34, face = "bold"),
        plot.title = element_text(size = 36, face = "bold"),
        plot.subtitle = element_text(size = 38, hjust = 0, face = "italic"),
        strip.background = element_blank(),
        strip.text = element_text(size = 36))


#boxplots relative abundance - nasal genera
top10G.names.bac.n = sort(tapply(taxa_sums(bac.nasaloralHMP_all.clean.3.depth_nasal), 
                                 tax_table((bac.nasaloralHMP_all.clean.3.depth_nasal))[, "Genus"], sum), TRUE)[1:10]
top10G.bac.n = subset_taxa(bac.nasaloralHMP_all.clean.3.depth_nasal, Genus %in% names(top10G.names.bac.n))

relab_genera = transform_sample_counts(top10G.bac.n, function(x) 100*x/sum(x))
relab_genera = tax_glom(relab_genera, taxrank = "Genus")
datframe_top10G.bac.n_relab = psmelt(relab_genera)
datframe_top10G.bac.n_relab$Log10_Abundance <- log10(datframe_top10G.bac.n_relab$Abundance+1)

ggplot(datframe_top10G.bac.n_relab, aes(x = Age, y=Log10_Abundance)) +
  geom_jitter(alpha = 0.25, color = "grey", size = 2, width = 0.3)+
  geom_boxplot(stat = "boxplot", outlier.shape = NA, linewidth = 1.5,
               aes(fill = Genus, alpha = Age))+
  facet_grid(~Genus, scales = "free", space = "free")+
  stat_compare_means(method = "wilcox.test",
                     paired = FALSE)+

#boxplots relative abundance - nasal species
top20S.names.bac.n = sort(tapply(taxa_sums(top10G.bac.n), 
                                 tax_table((top10G.bac.n))[, "Species"], sum), TRUE)[1:20]
top20S.bac.n = subset_taxa(top10G.bac.n, Species %in% names(top20S.names.bac.n))

relab_species= transform_sample_counts(top20S.bac.n, function(x) 100*x/sum(x))
relab_species = tax_glom(relab_species, taxrank = "Species")
datframe_top20S.bac.n_relab = psmelt(relab_species)
datframe_top20S.bac.n_relab$Log10_Abundance <- log10(1+datframe_top20S.bac.n_relab$Abundance)

ggplot(datframe_top20S.bac.n_relab, aes(x = Age, y=Log10_Abundance)) +
  geom_jitter(alpha = 0.25, color = "grey", size = 2, width = 0.3)+
  geom_boxplot(stat = "boxplot", outlier.shape = NA, linewidth = 1.5,
               aes(fill = Genus, alpha = Age))+
  facet_grid(~Species, scales = "free", space = "free")+
  stat_compare_means(method = "wilcox.test",
                     paired = FALSE)
  
#boxplots relative abundance - oral genera
top10G.names.bac.o = sort(tapply(taxa_sums(bac.nasaloralHMP_all.clean.3.depth_oral), 
                                 tax_table((bac.nasaloralHMP_all.clean.3.depth_oral))[, "Genus"], sum), TRUE)[1:10]
top10G.bac.o = subset_taxa(bac.nasaloralHMP_all.clean.3.depth_oral, Genus %in% names(top10G.names.bac.o))

relab_genera = transform_sample_counts(top10G.bac.o, function(x) 100*x/sum(x))
relab_genera = tax_glom(relab_genera, taxrank = "Genus")
datframe_top10G.bac.o_relab = psmelt(relab_genera)
datframe_top10G.bac.o_relab$Log10_Abundance <- log10(1+datframe_top10G.bac.o_relab$Abundance)

ggplot(datframe_top10G.bac.o_relab, aes(x = Age, y=Log10_Abundance)) +
  geom_jitter(alpha = 0.25, color = "grey", size = 2, width = 0.3)+
  geom_boxplot(stat = "boxplot", outlier.shape = NA, linewidth = 1.5,
               aes(fill = Genus, alpha = Age))+
  scale_alpha_discrete(range = c(0.4, 0.8))+
  facet_grid(~Genus, scales = "free", space = "free")+
  stat_compare_means(method = "wilcox.test",
                     paired = FALSE)

#boxplots relative abundance - oral species
top20S.names.bac.o = sort(tapply(taxa_sums(top10G.bac.o), 
                                 tax_table((top10G.bac.o))[, "Species"], sum), TRUE)[1:20]
top20S.bac.o = subset_taxa(top10G.bac.o, Species %in% names(top20S.names.bac.o))

relab_species= transform_sample_counts(top20S.bac.o, function(x) 100*x/sum(x))
relab_species = tax_glom(relab_species, taxrank = "Species")
datframe_top20S.bac.o_relab = psmelt(relab_species)
datframe_top20S.bac.o_relab$Log10_Abundance <- log10(1+datframe_top20S.bac.o_relab$Abundance)

ggplot(datframe_top20S.bac.o_relab, aes(x = Age, y=Log10_Abundance)) +
  geom_jitter(alpha = 0.25, color = "grey", size = 2, width = 0.3)+
  geom_boxplot(stat = "boxplot", outlier.shape = NA, linewidth = 1.5,
               aes(fill = Genus, alpha = Age))+
  facet_grid(~Species, scales = "free", space = "free")+
  stat_compare_means(method = "wilcox.test",
                     paired = FALSE)
 

