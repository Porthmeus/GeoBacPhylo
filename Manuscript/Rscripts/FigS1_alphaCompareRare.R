# Porthmeus
# 29.07.21

require(data.table)
require(DESeq2)
if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R)
require(ggplot2)
require(vegan)
require(phyloseq)
require(lme4)
require(lmerTest)
require(car)
require(lsmeans)
require(multcomp)
require(ggpubr)
require(cowplot)

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

# simplify the reproduction vector
map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]

# add the correct naming for the nutritional state of the water bodies
map_nut <- c("Hypereutrophic", "Eutrophic","Mesoeutrophic")
meta[["Nutrient_load"]] <- map_nut[metaLakes[meta$PopID,"NutrientLoad"]]
meta[["Water_body2"]] <- metaLakes[meta$PopID, "Waterbody2"]
meta[["Water_body"]] <- metaLakes[meta$PopID, "waterbody"]
meta[["Water_body"]] <- gsub("River","running",gsub("Lake","standing",meta[["Water_body"]]))
meta[["LocationIDshort"]] <- sapply(meta$LocationID, function(x) strsplit(x, split ="/")[[1]][2])

# reorder the sampling sites, so that they occur in a similar order as no the map
siteOrder <- c("M34","M72","M71","M70","M31","M67","M89","M90","M86","M85","M84","M79","M78","M109","M28","M26","M44","M108","M107","M83","R12")
meta[["SiteID"]] <- factor(meta[["PopID"]], levels = siteOrder)

# remove chloroplast ESVs from (probably) algea
sel <- which(tax[["Order"]] != "Chloroplast")
tax <- tax[sel,]
otu <- otu[rownames(tax),]

# create a rarecurve plot
samplePoints <- log2(range(colSums(otu)))
samplePoints <- 2^seq(from = log2(100), to = samplePoints[2], length= 100)
rarecurveMat <- rarefy(t(otu), sample = samplePoints)
rarecurveTab <- cbind(as.data.frame(rarecurveMat), ID = colnames(otu))
rarecurveTab <- reshape2::melt(rarecurveTab,id.var = "ID")
rarecurveTab[["N"]] <- as.numeric(gsub("^N","",rarecurveTab[["variable"]]))
# remove those values which exceed the actual sample size of the sample
rarecurveTab <- data.table(rarecurveTab)
sums <- colSums(otu)
for(id in names(sums)){
    rarecurveTab[ID==id & sums[id] < N, value := NA]
}
rarecurveTab <- merge(rarecurveTab, data.table(meta, keep.rownames = TRUE), by.x="ID", by.y = "rn")

p_rarecurv <- ggplot(rarecurveTab, aes(x= N/1000,y=value, group = ID, color = Species)) +
    geom_line() +
    facet_wrap(~SiteID,ncol = 7) +
    ylab("Rarefraction")+
    xlab("Number of rarified samples x1000")+
    theme_bw()
p_rarecurv

# create a boxplot of the rarified values
rare <- meta
rare[["rare"]] <- as.data.frame(rarefy(t(otu), min(colSums(otu))))[rownames(meta),1]

p_rareBox <- ggplot(rare, aes(x = SiteID, y = rare, fill = SiteID)) +
    geom_boxplot() + 
    ylab("Rarefraction") +
    theme_bw()+
    theme(legend.position = "none") 
p_rareBox

# check the correlation between Shannon, Chao and Rarefraction
phy <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
                tax_table(as.matrix(tax)),
                sample_data(meta))

alpha <- estimate_richness(phy,measures = c("Shannon","Chao1"))
alpha <- merge(alpha, rare, by=0)

p_corShan <- ggplot(alpha, aes(x = rare, y = Shannon)) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor(method = "pearson",label.y = 6)+
    xlab("Rarefraction")+
    theme_bw()
p_corShan

p_corChao <- ggplot(alpha, aes(x = rare, y = Chao1)) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor(method = "pearson")+
    xlab("Rarefraction")+
    theme_bw()
p_corChao


p_bottom <- plot_grid(p_corShan, p_corChao, labels = c("C","D"))
pp <- plot_grid(p_rarecurv, p_rareBox, p_bottom, ncol =1, rel_heights = c(2,1,1), labels=c("A","B",NA))
pp
ggsave(pp, file = "figures/FigS1_alphaCompareRare.pdf", width = 8, height = 10)
