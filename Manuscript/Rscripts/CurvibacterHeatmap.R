# Porthmeus
# 16.03.21

require(data.table)
require(DESeq2)
if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R)
require(ggplot2)
require(vegan)
require(phyloseq)
require(umap)
require(pals)
require(lme4)


AbFilt_FT = "../data/dada2_AbFilt_FT.csv"
AbFilt_Tax = "../data/dada2_AbFilt_Tax.csv"
AbFilt_meta = "../data/dada2_AbFilt_meta.csv"
AbFilt_fasta_tree = "../data/dada2_AbFilt_fasta_tree.qza"
MetaLakes = "../MetaDataLakes.csv"

otu <- read.csv(AbFilt_FT, row.names = 1)
# remove the smallest 5% (<3368 reads) from the samples - ~14
sumRC <- apply(otu,2,sum)
otu <- otu[,sumRC > quantile(sumRC,0.05)]
tax <- read.csv(AbFilt_Tax, row.names=1)
meta <- read.csv(AbFilt_meta, row.names=1)
metaLakes <- read.csv(MetaLakes,row.names=1)

otu <- read.csv(AbFilt_FT, row.names = 1)
tax <- read.csv(AbFilt_Tax, row.names=1)
meta <- read.csv(AbFilt_meta, row.names=1)


# simplify the reproduction vector
map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]

# add the correct naming for the nutritional state of the water bodies
map_nut <- c("Hypereutr.", "Eutr.","Mesoeut.")
meta[["Nutrient_load"]] <- map_nut[metaLakes[meta$PopID,"NutrientLoad"]]
meta[["Water_body2"]] <- metaLakes[meta$PopID, "Waterbody2"]
meta[["Water_body"]] <- metaLakes[meta$PopID, "waterbody"]
meta[["Water_body"]] <- gsub("River","running",gsub("Lake","standing",meta[["Water_body"]]))
meta[["LocationIDshort"]] <- sapply(meta$LocationID, function(x) strsplit(x, split ="/")[[1]][2])

# reorder the sampling sites, so that they occur in a similar order as no the map
siteOrder <- c("M34","M72","M71","M70","M31","M67","M89","M90","M86","M85","M84","M79","M78","M109","M28","M26","M44","M108","M107","M83","R12")
meta[["SiteID"]] <- factor(meta[["PopID"]], levels = siteOrder)

# add a single count to the OTU matrix, stabilize variance, and shift the values to positive values
otu2 <- otu+1
otu_vst <- varianceStabilizingTransformation(as.matrix(otu2), fitType = "local")
if(min(otu_vst) <0){
    otu_vst <- otu_vst + min(otu_vst)*-1
}
otu_depthNorm <- (otu+1) * estimateSizeFactorsForMatrix(otu+1)
otu_depthNorm <- otu_depthNorm - min(otu_depthNorm)

curvi <- c("59e1be396c0f4128232640255c7d994f","42bdceb495e3fd7c545d1a13a50ec1b7")
meta[["curviCounts"]] <- apply(otu_vst[curvi,rownames(meta)],2,sum)
meta[["curviCounts_raw"]] <- apply(otu_depthNorm[curvi,rownames(meta)],2,sum)
meta <- data.table(meta)


kruskal.test(data=meta,
   curviCounts ~ Species)
kruskal.test(data=meta,
   curviCounts ~ PopID)
kruskal.test(data=meta,
          curviCounts ~ Reproduction)
kruskal.test(data=meta,
   curviCounts ~ Nutrient_load)
kruskal.test(data=meta,
   curviCounts ~ Water_body)

with(meta,
   pairwise.wilcox.test(curviCounts, Species,p.adjust="BH"))
with(meta,
   pairwise.wilcox.test(curviCounts, PopID,p.adjust="BH"))
with(meta,
   pairwise.wilcox.test(curviCounts, Reproduction,p.adjust="BH"))


p_pop <- ggplot(meta, aes(x= PopID, y = curviCounts, fill = PopID)) +
    geom_boxplot() +
    theme_bw() +
    theme(legend.position = "None") +
    ylab("Curvibacter norm. RC")+
#    ylab("Variance stabalized and normalized RC vor Curvibacter") +
    xlab("Site ID")
    #theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))
p_waterBody <- ggplot(meta, aes(x= Water_body, y = curviCounts, fill = Water_body)) +
    geom_boxplot() +
    theme_bw() +
    theme(legend.position = "None") +
    ylab("Curvibacter norm. RC")+
    #ylab("Variance stabalized and normalized RC vor Curvibacter") +
    xlab("Water body")
    #theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))
p_nutrientLoad <- ggplot(meta, aes(x= Nutrient_load, y = curviCounts, fill = Nutrient_load)) +
    geom_boxplot() +
    theme_bw() +
    theme(legend.position = "Nutrient load") +
    ylab("Curvibacter norm. RC")+
    #ylab("Variance stabalized and normalized RC vor Curvibacter") +
    xlab("Nutrient load")
    #theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))
p_species <- ggplot(meta, aes(x= Species, y = curviCounts, fill = Species)) +
    geom_boxplot() +
    theme_bw()+
    theme(legend.position = "None") +
    ylab("Curvibacter norm. RC")
    #ylab("Variance stabalized and normalized RC vor Curvibacter")+
p_reproduction <- ggplot(meta, aes(x= Reproduction, y = curviCounts, fill = Reproduction)) +
    geom_boxplot() +
    theme_bw()+
    theme(legend.position = "None") +
    ylab("Curvibacter norm. RC")
    #ylab("Variance stabalized and normalized RC vor Curvibacter")+

p_bottom <- cowplot::plot_grid(p_waterBody, p_nutrientLoad, p_species, p_reproduction, ncol = 2, rel_widths = c(3,3,3,3), labels = LETTERS[2:5])
pp <- cowplot::plot_grid(p_pop, p_bottom, ncol = 1, labels = c("A",NA), rel_heights= c(1,2))

ggsave(pp, file = "figures/Fig6_Curvibacter_counts.pdf",width =8 , height = 7)
