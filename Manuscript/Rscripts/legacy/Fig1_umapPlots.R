# Porthmeus
# 29.01.21

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


AbFilt_FT = "../data/dada2_AbFilt_FT.csv"
AbFilt_Tax = "../data/dada2_AbFilt_Tax.csv"
AbFilt_meta = "../data/dada2_AbFilt_meta.csv"
AbFilt_fasta_tree = "../data/dada2_AbFilt_fasta_tree.qza"
MetaLakes = "../MetaDataLakes.csv"


otu <- read.csv(AbFilt_FT, row.names = 1)
tax <- read.csv(AbFilt_Tax, row.names=1)
meta <- read.csv(AbFilt_meta, row.names=1)
metaLakes <- read.csv(MetaLakes,row.names=1)

# simplify the reproduction vector
map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]
map_nut <- c("Hypertrophic", "Eutrophic","Mesotrophic")
meta[["Nutrient_load"]] <- map_nut[metaLakes[meta$PopID,"NutrientLoad"]]


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
                   
# calculate bray-curtis and jaccard distance
brayCurt <- distance(phy_vst, method = "bray")
jac <- distance(phy_vst, method = "jaccard", binary=TRUE)

# calculate UMAP
set.seed(4524)
umapBray <- umap(as.matrix(brayCurt), input = "dist")
umapJac <- umap(as.matrix(jac), input = "dist")
set.seed(-1)


# plot ordination
layoutBray <- as.data.frame(umapBray$layout)
colnames(layoutBray) <- c("UMAP1","UMAP2")

layoutJac <- as.data.frame(umapJac$layout)
colnames(layoutJac) <- c("UMAP1","UMAP2")

p_data <- merge(meta,layoutBray, by=0)
rownames(p_data) <- p_data[["Row.names"]]
p_data <- merge(p_data[,-1], layoutJac, by = 0, suffixes = c("_bray","_jac"))

# plot the overview
p_bray <- ggplot(p_data, aes(x=UMAP1_bray, y=UMAP2_bray)) +
    geom_point(aes(color = PopID)) +
    theme_bw() +
    ggtitle("Bray-Curtis") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 

ggsave(p_bray , 
       file = "figures/UMAP_braycurtis.pdf",
       width = 5,
       height= 4
)

p_jac <- ggplot(p_data, aes(x=UMAP1_jac, y=UMAP2_jac)) +
    geom_point(aes(color = PopID)) +
    theme_bw() +
    ggtitle("Jaccard") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 

ggsave(p_jac , 
       file = "figures/UMAP_jaccard.pdf",
       width = 5,
       height= 4
)

p_bray <- ggplot(p_data, aes(x=UMAP1_bray, y=UMAP2_bray)) +
    geom_point(aes(color = Nutrient_load)) +
    theme_bw() +
    ggtitle("Bray-Curtis") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 

ggsave(p_bray , 
       file = "figures/UMAP_braycurtis_NutrientLoad.pdf",
       width = 5,
       height= 4
)

p_jac <- ggplot(p_data, aes(x=UMAP1_jac, y=UMAP2_jac)) +
    geom_point(aes(color = Nutrient_load)) +
    theme_bw() +
    ggtitle("Jaccard") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 

ggsave(p_jac , 
       file = "figures/UMAP_jaccard_NutrientLoad.pdf",
       width = 5,
       height= 4
)

# quickly do a PCA on it
pca_jac <- prcomp(as.matrix(jac))
pca_bray <- prcomp(as.matrix(brayCurt))
pca_jac_var <- (pca_jac[["sdev"]])^2/sum((pca_jac[["sdev"]])^2)
pca_bray_var <- (pca_bray[["sdev"]])^2/sum((pca_bray[["sdev"]])^2)
pca_jacPC <- pca_jac[["x"]][,c("PC1","PC2")]
pca_jacPC <- data.table(pca_jacPC, keep.rownames=TRUE)
pca_brayPC <- pca_bray[["x"]][,c("PC1","PC2")]
pca_brayPC <- data.table(pca_brayPC, keep.rownames=TRUE)
p_data2 <- merge(data.table(p_data), pca_jacPC, by.x="Row.names", by.y="rn")
p_data3 <- merge(data.table(p_data), pca_brayPC, by.x="Row.names", by.y="rn")
p_data2[,Distance := "Jaccard"]
p_data3[,Distance := "Bray-Curtis"]
p_data2 <- rbind(p_data2,p_data3)

p_pca <- ggplot(p_data2[Distance == "Jaccard",], aes(x=PC1, y=PC2, color = Nutrient_load)) + 
    geom_point(size=3) +
    ylab(paste0("PC2 (",round(pca_jac_var[2]*100,digits=2),"%)"))+
    xlab(paste0("PC1 (",round(pca_jac_var[1]*100,digits=2),"%)"))+
    ggtitle("Jaccard") + 
    theme_bw()
ggsave(file = "figures/PCA_nutrientLoadAll_jaccard.pdf",
       width = 4.8,
       height = 4,
       p_pca)

p_pca <- ggplot(p_data2[Distance == "Bray-Curtis",], aes(x=PC1, y=PC2, color = Nutrient_load)) + 
    geom_point(size=3) +
    ylab(paste0("PC2 (",round(pca_bray_var[2]*100,digits=2),"%)"))+
    xlab(paste0("PC1 (",round(pca_bray_var[1]*100,digits=2),"%)"))+
    ggtitle("Bray-Curtis") + 
    theme_bw()
ggsave(file = "figures/PCA_nutrientLoadAll_bray.pdf",
       width = 4.8,
       height = 4,
       p_pca)


p_bray <- ggplot(p_data, aes(x=UMAP1_bray, y=UMAP2_bray)) +
    geom_point(aes(color = Species)) +
    theme_bw() +
    ggtitle("Bray-Curtis") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 

ggsave(p_bray , 
       file = "figures/UMAP_braycurtis_SpeciesAll.pdf",
       width = 4.4,
       height= 4
)

p_jac <- ggplot(p_data, aes(x=UMAP1_jac, y=UMAP2_jac)) +
    geom_point(aes(color = Species)) +
    theme_bw() +
    ggtitle("Jaccard") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 

ggsave(p_jac , 
       file = "figures/UMAP_jaccard_SpeciesAll.pdf",
       width = 4.4,
       height= 4
)


# find the populations with all species and plot the umap of these
p_dt <- data.table(p_data)
p_sel <- p_dt[, length(unique(Species)), by =PopID]
pops <- p_sel[V1 >1, PopID]

p_jac_single <- ggplot(p_dt[PopID %in% pops,], aes(x=UMAP1_jac, y=UMAP2_jac)) +
    geom_point(aes(color = Species)) +
    theme_bw() +
    facet_wrap(~PopID, scale = "free")+
    ggtitle("Jaccard") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 

ggsave(p_jac_single,
       file = "figures/UMAP_jaccard_Species.pdf",
       width = 7,
       height = 6)

p_bray_single <- ggplot(p_dt[PopID %in% pops,], aes(x=UMAP1_bray, y=UMAP2_bray)) +
    geom_point(aes(color = Species)) +
    theme_bw() +
    facet_wrap(~PopID, scale = "free")+
    ggtitle("Bray-Curtis") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 


ggsave(p_jac_single,
       file = "figures/UMAP_bray_Species.pdf",
       width = 7,
       height = 6)


# do the same for the sexual reproduction
p_sel <- p_dt[Species=="OLI", length(unique(Reproduction)), by =PopID]
pops <- p_sel[V1 >1, PopID]


p_jac_single <- ggplot(p_dt[PopID %in% pops,], aes(x=UMAP1_jac, y=UMAP2_jac)) +
    geom_point(aes(color = Reproduction)) +
    theme_bw() +
    facet_wrap(~PopID, scale = "free")+
    ggtitle("Jaccard") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 

ggsave(p_jac_single,
       file = "figures/UMAP_jaccard_Reproduction.pdf",
       width = 7,
       height = 6)

p_bray_single <- ggplot(p_dt[PopID %in% pops,], aes(x=UMAP1_bray, y=UMAP2_bray)) +
    geom_point(aes(color = Reproduction)) +
    theme_bw() +
    facet_wrap(~PopID, scale = "free")+
    ggtitle("Bray-Curtis") +
    ylab ("UMAP2") +
    xlab ("UMAP1") 


ggsave(p_jac_single,
       file = "figures/UMAP_bray_Reproduction.pdf",
       width = 7,
       height = 6)


# do the statistics

cat("\n\n\n Population effect on beta diversity\n")
meta <- meta[rownames(as.matrix(jac)),]
adonis(formula = jac ~ meta[,"PopID"], perm = 500)
adonis(formula = brayCurt ~ meta[,"PopID"], perm = 500)

cat("\n\n\n Species effect on beta diversity\n")
adonis(formula = jac ~ meta[,"Species"], strata = meta[,"PopID"], perm = 500)
adonis(formula = brayCurt ~ meta[,"Species"],strata = meta[,"PopID"], perm = 500)

cat("\n\n\n Nutrient load effect on beta diversity\n")
adonis(formula = jac ~ meta[,"Nutrient_load"], perm = 500)
adonis(formula = brayCurt ~ meta[,"Nutrient_load"], perm = 500)

# remove M89

cat("\n\n\n Species effect without Population M89\n")
sel <- rownames(meta[meta$PopID != "M89",])

adonis(formula = as.dist(as.matrix(jac)[sel,sel]) ~ meta[sel,"Species"], strata = meta[sel,"PopID"], perm = 500)
adonis(formula = as.dist(as.matrix(brayCurt)[sel,sel]) ~ meta[sel,"Species"], strata = meta[sel,"PopID"], perm = 500)

# remove H. circ

cat("\n\n\n Species effect without H. circumcincta\n")
sel <- rownames(meta[meta$Species != "CIR",])

adonis(formula = as.dist(as.matrix(jac)[sel,sel]) ~ meta[sel,"Species"], strata = meta[sel,"PopID"], perm = 500)
adonis(formula = as.dist(as.matrix(brayCurt)[sel,sel]) ~ meta[sel,"Species"], strata = meta[sel,"PopID"], perm = 500)



# test Species effect in specific Populations
for(pop in unique(meta$PopID)){
    sel <- meta$PopID == pop
    tblSpec <- table(meta[sel,"Species"])
    if(min(tblSpec) > 2 & length(tblSpec) >=2){
        cat(paste0("\n\n\n Species effect on beta diversity within ",pop,"\n"))
        print(adonis(formula = as.dist(as.matrix(jac)[sel,sel]) ~ meta[sel,"Species"], strata = meta[sel,"PopID"], perm = 500))
        print(adonis(formula = as.dist(as.matrix(brayCurt)[sel,sel]) ~ meta[sel,"Species"], strata = meta[sel,"PopID"], perm = 500))
    }
}


