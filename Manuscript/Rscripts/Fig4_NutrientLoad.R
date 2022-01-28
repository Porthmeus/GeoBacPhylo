# Porthmeus
# 12.05.21

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
map_nut <- c("Hypereutr.", "Eutr.","Mesoeutr.")
meta[["Nutrient_load"]] <- factor(map_nut[metaLakes[meta$PopID,"NutrientLoad"]], levels = map_nut)
meta[["Water_body2"]] <- metaLakes[meta$PopID, "Waterbody2"]
meta[["Water_body"]] <- metaLakes[meta$PopID, "waterbody"]
meta[["Water_body"]] <- gsub("River","running",gsub("Lake","standing",meta[["Water_body"]]))

# reorder the sampling sites, so that they occur in a similar order as no the map
siteOrder <- c("M34","M72","M71","M70","M31","M67","M89","M90","M86","M85","M84","M79","M78","M109","M28","M26","M44","M108","M107","M83","R12")
meta[["SiteID"]] <- factor(meta[["PopID"]], levels = siteOrder)

# remove chloroplast ESVs from (probably) algea
sel <- which(tax[["Order"]] != "Chloroplast")
tax <- tax[sel,]
otu <- otu[rownames(tax),]

# add a single count to the OTU matrix, stabilize variance, and shift the values to positive values
otu2 <- otu+1
otu_vst <- varianceStabilizingTransformation(as.matrix(otu2), fitType = "local")
if(min(otu_vst) <0){
    otu_vst <- otu_vst + min(otu_vst)*-1
}

phy_vst <- phyloseq(otu_table(otu_vst,
                              taxa_are_rows = TRUE),
                              tax_table(as.matrix(tax)),
                              sample_data(meta))


# calculate bray curtis
brayCurt <- distance(phy_vst, method = "bray")

# calculate UMAP
set.seed(4524)
umapBray <- umap(as.matrix(brayCurt), input = "dist")
set.seed(-1)

# plot umap
layoutBray <- as.data.frame(umapBray$layout)
colnames(layoutBray) <- c("UMAP1","UMAP2")
p_data <- merge(meta,layoutBray, by=0)
rownames(p_data) <- p_data[["Row.names"]]

p_umap <- ggplot(p_data, aes(x=UMAP1, y=UMAP2, color = Nutrient_load)) +
    geom_point(size = 2) +
    theme_bw() +
    theme(legend.position = "none")
p_umap

# plot PCA
pca <- prcomp(brayCurt,scale=TRUE,center=TRUE)
pca_var <- (pca$sdev)^2/sum(pca$sdev^2)
p_data <- merge(meta, pca$x[,1:2], by=0)
rownames(p_data) <- p_data[["Row.names"]]

p_pca <- ggplot(p_data, aes(x=PC1, y=PC2, color = Nutrient_load)) +
    geom_point(size = 2) +
    ylab(paste0("PC2 (",round(pca_var[2]*100,2),"%)")) +
    xlab(paste0("PC1 (",round(pca_var[1]*100,2),"%)")) +
    theme_bw() +
    theme(legend.position = "none")
p_pca

# plot distance boxplot
testMeta <- data.table(meta, keep.rownames=TRUE)
setkey(testMeta, "rn")
testMeta <- testMeta[rownames(as.matrix(brayCurt)),]

perm <- how(nperm = 999, blocks = testMeta$Species, plots = Plots(strata= testMeta$Reproduction))
adn <- adonis2(brayCurt~Nutrient_load+Water_body, data=testMeta, by="margin", permutations = perm)
print(adn)

brayCurt_long <-as.vector(brayCurt)
n <- as.data.frame(t(combn(rownames(as.matrix(brayCurt)),2)))
brayCurt_long <- cbind(n,brayCurt = brayCurt_long)
colnames(brayCurt_long) <- c("Sample1","Sample2","Bray_Curtis")
brayCurt_long <- data.table(brayCurt_long)

# add the information if within or between population IDs
brayCurt_long[,Nutrient_load := c("within","between")[as.integer(meta[Sample1, "Nutrient_load"] != meta[Sample2, "Nutrient_load"])+1]]
brayCurt_long[,Nutrient_load := factor(Nutrient_load, levels = c("within","between"))]
Nutrient_load_long <- paste(meta[brayCurt_long[,Sample1],"Nutrient_load"],meta[brayCurt_long[,Sample2], "Nutrient_load"], sep = "-")
n <- combn(unique(meta[["Nutrient_load"]]),2)
crct <- apply(n,2,paste,collapse="-")
wrng <- apply(n[c(2,1),],2,paste,collapse="-")
for(i in 1:ncol(n)){
    Nutrient_load_long <- gsub(wrng[i],crct[i],Nutrient_load_long)
}
brayCurt_long[,Nutrient_load_long:=Nutrient_load_long]


p_boxDist <- ggplot(data = brayCurt_long, aes(x=Nutrient_load_long, y=Bray_Curtis)) +
    geom_boxplot() +
    ylab("Bray-Curtis distance") +
    theme_bw() +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
p_boxDist

# plot alpha diversity
# create a phyloseq object
phy <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
                tax_table(as.matrix(tax)),
                sample_data(meta))

# get the richness estimates and merge them into a k
measures <- c("Chao1","Shannon")
alphaDiv <- estimate_richness(phy,measures=measures)
alphaDiv <- alphaDiv[,measures]
alphaDiv[["Sample"]] <- rownames(alphaDiv)
alphaDiv <- reshape2::melt(alphaDiv, id.vars = "Sample")
alphaDiv <- merge(alphaDiv, meta, by.x="Sample", by.y=0, all.x=TRUE)

p_alpha <- ggplot(alphaDiv, aes(x=Nutrient_load, y=value, fill = Nutrient_load)) +
    geom_boxplot() +
    facet_wrap(~variable, scale ="free_y", ncol = 1) +
    theme_bw() +
    ylab("Value")+
    xlab("Nutrient load")+
    labs(fill = "Nutrient load")+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
p_alpha

#p_toprow <- cowplot::plot_grid(p_umap,p_pca,p_boxDist,ncol=3,labels=LETTERS[1:3], rel_widths=c(2,2,1))
p_toprow <- cowplot::plot_grid(p_umap,p_pca,ncol=2,labels=LETTERS[1:2], rel_widths=c(2,2))
p_leftrow <- cowplot::plot_grid(p_umap,p_pca,ncol=1,labels=LETTERS[1:2], rel_widths=c(2,2))
p_botrow <- cowplot::plot_grid(p_boxDist,p_alpha, ncol = 2, labels = LETTERS[3:4], rel_widths = c(1,1.3))
p_rightrow <- cowplot::plot_grid(p_boxDist,p_alpha, ncol = 1, labels = LETTERS[3:4], rel_heights = c(1,1))
#p_all <- cowplot::plot_grid(p_toprow, p_alpha,ncol=1,labels=c("","D"), rel_heights = c(1,1))
p_all <- cowplot::plot_grid(p_leftrow,p_rightrow,ncol=2, rel_heights = c(1,1))
p_all <- cowplot::plot_grid(p_toprow,p_botrow,ncol=1, rel_heights = c(1,1.3))
p_all

ggsave(p_all,file = "figures/Fig4_NutrientLoad.pdf", width = 8, height= 8)
