# Porthmeus
# 20.04.21

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
require(plotly)
require(data.table)


AbFilt_FT = "../data/dada2_AbFilt_FT.csv"
AbFilt_Tax = "../data/dada2_AbFilt_Tax.csv"
AbFilt_meta = "../data/dada2_AbFilt_meta.csv"
AbFilt_fasta_tree = "../data/dada2_AbFilt_fasta_tree.qza"
MetaLakes = "../MetaDataLakes.csv"


otu <- read.csv(AbFilt_FT, row.names = 1)
tax <- read.csv(AbFilt_Tax, row.names=1)
meta <- read.csv(AbFilt_meta, row.names=1)

# simplify the reproduction vector
map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]

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
meta <- data.table(meta, keep.rownames=TRUE)

# calculate bray-curtis and jaccard distance
brayCurt <- phyloseq::distance(phy_vst, method = "bray")
jac <- phyloseq::distance(phy_vst, method = "jaccard")



# calculate the PCA for all populations where the separation is visible
corPop <- c("M70","M83","M107","M79","M67")
tab <- meta[PopID %in% corPop,.(PC1 = 0, PC2 = 0, Species, Reproduction, rn, PopID, Distance = "Jaccard")]
tab2 <- meta[PopID %in% corPop,.(PC1 = 0, PC2 = 0, Species, Reproduction, rn, PopID, Distance = "Bray-Curtis")]
setkey(tab, "rn")
setkey(tab2, "rn")
for(pop in corPop){
    sel <- tab[PopID == pop, rn]
    pca <- prcomp(as.matrix(jac)[sel,sel])
    tab[sel, PC1 := pca[["x"]][,"PC1"]]
    tab[sel, PC2 := pca[["x"]][,"PC2"]]
    pca <- prcomp(as.matrix(brayCurt)[sel,sel])
    tab2[sel, PC1 := pca[["x"]][,"PC1"]]
    tab2[sel, PC2 := pca[["x"]][,"PC2"]]
}
tab <- rbind(tab,tab2)

p <- ggplot(tab, aes(x=PC1,y=PC2, color = Species, shape = Reproduction)) +
    geom_point(size =3) +
    facet_grid(PopID~Distance)+
    theme_bw()
ggsave(p,
       file = "figures/PCA_speciesCorrectClustering.pdf",
       width = 8,
       height = 10)
ggsave(p,
       file = "Panels/PCA_speciesCorrectClustering.png",
       dpi = 600,
       width = 8,
       height = 10)


icorPop <- c("M90","M28")
tab <- meta[PopID %in% icorPop,.(PC1 = 0, PC2 = 0, Species, Reproduction, rn, PopID, Distance = "Jaccard")]
tab2 <- meta[PopID %in% icorPop,.(PC1 = 0, PC2 = 0, Species, Reproduction, rn, PopID, Distance = "Bray-Curtis")]
setkey(tab, "rn")
setkey(tab2, "rn")
for(pop in icorPop){
    sel <- tab[PopID == pop, rn]
    pca <- prcomp(as.matrix(jac)[sel,sel])
    tab[sel, PC1 := pca[["x"]][,"PC1"]]
    tab[sel, PC2 := pca[["x"]][,"PC2"]]
    pca <- prcomp(as.matrix(brayCurt)[sel,sel])
    tab2[sel, PC1 := pca[["x"]][,"PC1"]]
    tab2[sel, PC2 := pca[["x"]][,"PC2"]]
}
tab <- rbind(tab,tab2)

p <- ggplot(tab, aes(x=PC1,y=PC2, color = Species, shape = Reproduction)) +
    geom_point(size =3) +
    facet_grid(PopID~Distance)+
    theme_bw()
ggsave(p,
       file = "figures/PCA_speciesIncorrectClustering.pdf",
       width = 8,
       height = 4)
ggsave(p,
       file = "panels/PCA_speciesIncorrectClustering.png",
       dpi = 600,
       width = 8,
       height = 4)

