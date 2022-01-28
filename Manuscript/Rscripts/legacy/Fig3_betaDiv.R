# Porthmeus
# 12.02.21

require(ggplot2)
require(phyloseq)
require(qiime2R)
require(DESeq2)
require(data.table)
require(cowplot)
require(Hmisc)


otu <- read.csv("../data/dada2_AbFilt_FT.csv", row.names = 1)
tax <- read.csv("../data/dada2_AbFilt_Tax.csv", row.names =1)
meta <- read.csv("../data/dada2_AbFilt_meta.csv", row.names=1)
tree <- read_qza("../data/dada2_AbFilt_fasta_tree.qza")$data
meta_lakes <- read.csv("../MetaDataLakes.csv")
rownames(meta_lakes) <- meta_lakes[["SiteID"]]

map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]
map_nut <- c("Hypertrophic","Eutrophic","Mesotrophic")
meta[["Nutrient_load"]] <- factor(map_nut[meta_lakes[meta[["PopID"]],"NutrientLoad"]], levels = map_nut)
meta_short <- meta[,c("Species","Reproduction","PopID")]
meta_short[["sample"]] <- rownames(meta_short) 
otu2 <- otu+1

otu_vst <- varianceStabilizingTransformation(as.matrix(otu2), fitType = "local")
if(min(otu_vst) <0){
    otu_vst <- otu_vst + min(otu_vst)*-1
}
phy_vst <- phyloseq(otu_table(otu_vst,
                              taxa_are_rows = TRUE),
                    tax_table(as.matrix(tax)),
                    sample_data(meta),
                    phy_tree(tree))




# calculate the combinations to remove duplicates from distance matrix
cmbs <- combn(meta_short[,"sample"], 2)
cmbs <- apply(cmbs,2,paste, collapse ="")

distMs <- c( "bray","jaccard")
mtrxs <- list()
for(meas in distMs){
         dst <- as.matrix(phyloseq::distance(phy_vst, method = meas))
         dst <- reshape2::melt(dst, varnames= c("sample.x","sample.y"))
         dst <- data.table(dst)
         dst <- dst[paste0(sample.x,sample.y) %in% cmbs,]
         dst <- merge(dst, data.table(meta_short), by.x = "sample.x", by.y="sample")
         dst <- merge(dst, data.table(meta_short), by.x = "sample.y", by.y="sample", suffixes= c(".x",".y"))
         dst <- cbind(dst,
              dst[,.(Population_inner = as.character(PopID.x == PopID.y), 
                     Species_inner = as.character(Species.x == Species.y),
                     Reproduction_inner = as.character(Reproduction.x == Reproduction.y))])
         dst[Population_inner =="FALSE",Species_inner := "between sites"]
         dst[Population_inner =="FALSE",Reproduction_inner :="between sites"]
         dst[Population_inner == "TRUE", Population_inner := "within"]
         dst[Population_inner == "FALSE", Population_inner := "between"]
         dst[Species_inner == "TRUE", Species_inner := "within"]
         dst[Species_inner == "FALSE", Species_inner := "between"]
         dst[Reproduction_inner == "TRUE", Reproduction_inner := "within"]
         dst[Reproduction_inner == "FALSE", Reproduction_inner := "between"]
         dst[,Species_inner := factor(Species_inner, levels = c("within","between","between sites"))]
         dst[,Reproduction_inner := factor(Reproduction_inner, levels = c("within","between","between sites"))]
         dst[,Population_inner := factor(Population_inner, levels = c("within","between"))]
#
        dst <- cbind(dst,
                     Species_dist = "between sites")
        dst[Population_inner == "within" & Species.x == "OLI" & Species.y == "OLI", Species_dist := "OLI-OLI"]
        dst[Population_inner == "within" & Species.x == "CIR" & Species.y == "CIR", Species_dist := "CIR-CIR"]
        dst[Population_inner == "within" & Species.x == "VUL" & Species.y == "VUL", Species_dist := "VUL-VUL"]
        dst[Population_inner == "within" & Species.x == "OLI" & Species.y == "CIR", Species_dist := "OLI-CIR"]
        dst[Population_inner == "within" & Species.y == "OLI" & Species.x == "CIR", Species_dist := "OLI-CIR"]
        dst[Population_inner == "within" & Species.x == "OLI" & Species.y == "VUL", Species_dist := "OLI-VUL"]
        dst[Population_inner == "within" & Species.y == "OLI" & Species.x == "VUL", Species_dist := "OLI-VUL"]
        dst[Population_inner == "within" & Species.x == "CIR" & Species.y == "VUL", Species_dist := "VUL-CIR"]
        dst[Population_inner == "within" & Species.y == "CIR" & Species.x == "VUL", Species_dist := "VUL-CIR"]
        dst[,Species_dist := factor(Species_dist, levels = c("OLI-OLI","VUL-VUL","CIR-CIR","OLI-VUL","OLI-CIR","VUL-CIR","between sites"))]
#
        dst <- cbind(dst,
                     Reproduction_dist = "between sites")
        dst[Population_inner == "within" & Reproduction.x == "ASEX" & Reproduction.y == "ASEX", Reproduction_dist := "ASEX-ASEX"]
        dst[Population_inner == "within" & Reproduction.x == "NR" & Reproduction.y == "NR", Reproduction_dist := "NR-NR"]
        dst[Population_inner == "within" & Reproduction.x == "SEX" & Reproduction.y == "SEX", Reproduction_dist := "SEX-SEX"]
        dst[Population_inner == "within" & Reproduction.x == "ASEX" & Reproduction.y == "NR", Reproduction_dist := "ASEX-NR"]
        dst[Population_inner == "within" & Reproduction.y == "ASEX" & Reproduction.x == "NR", Reproduction_dist := "ASEX-NR"]
        dst[Population_inner == "within" & Reproduction.x == "ASEX" & Reproduction.y == "SEX", Reproduction_dist := "ASEX-SEX"]
        dst[Population_inner == "within" & Reproduction.y == "ASEX" & Reproduction.x == "SEX", Reproduction_dist := "ASEX-SEX"]
        dst[Population_inner == "within" & Reproduction.x == "NR" & Reproduction.y == "SEX", Reproduction_dist := "SEX-NR"]
        dst[Population_inner == "within" & Reproduction.y == "NR" & Reproduction.x == "SEX", Reproduction_dist := "SEX-NR"]
        dst[,Reproduction_dist := factor(Reproduction_dist, levels = c("ASEX-ASEX","SEX-SEX","NR-NR","ASEX-SEX","ASEX-NR","SEX-NR","between sites"))]
#
         mtrxs[[meas]] <- dst
}

plots <- list()
for(mtrx in distMs){
    if(mtrx == "bray"){
        ttl <- "Bray-Curtis"
    } else {
        ttl <- capitalize(mtrx)
    }
    inner_plots <- list()
    for(meas in grep("_inner", colnames(mtrxs[[mtrx]]), value =TRUE)){
        pp <- ggplot(mtrxs[[mtrx]], aes_string(x=meas, y = "value")) +
            geom_boxplot() +
            theme_bw() +
            ggtitle(ttl) +
            ylab("Distance") +
            xlab(gsub("_inner","", meas)) +
            theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1))
        inner_plots[[meas]] <- pp
    }
    plots[[mtrx]] <- plot_grid(plotlist=inner_plots, ncol = 3)
    ggsave(plots[[mtrx]],
           file = paste0("figures/Fig3_distance_",mtrx,".pdf"),
           height = 4, width = 6)
}

plots2 <- list()
for(mtrx in distMs){
    if(mtrx == "bray"){
        ttl <- "Bray-Curtis"
    } else {
        ttl <- capitalize(mtrx)
    }
    inner_plots2 <- list()
    for(meas in grep("_dist", colnames(mtrxs[[mtrx]]), value =TRUE)){
        pp <- ggplot(mtrxs[[mtrx]], aes_string(x=meas, y = "value")) +
            geom_boxplot() +
            theme_bw() +
            ggtitle(ttl) +
            ylab("Distance") +
            xlab(gsub("_inner","", meas)) +
            theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1))
        inner_plots2[[meas]] <- pp
    }
    plots2[[mtrx]] <- plot_grid(plotlist=inner_plots2, ncol = 1)
    ggsave(plots2[[mtrx]],
           file = paste0("figures/Fig3_distanceSingular_",mtrx,".pdf"),
           height = 8, width = 6)
}


p_popDist_jac <- ggplot(mtrxs[["jaccard"]], aes(x=Population_inner,y=value)) +
    geom_boxplot() +
    theme_bw() +
    ylab("Distance") +
    xlab("Sampling site") +
    ggtitle("Jaccard")
p_popDist_bray <- ggplot(mtrxs[["bray"]], aes(x=Population_inner,y=value)) +
    geom_boxplot() +
    theme_bw() +
    ylab("Distance") +
    xlab("Sampling site") +
    ggtitle("Bray-Curtis")
pp <- plot_grid(p_popDist_jac, p_popDist_bray, ncol = 2)
ggsave(file ="figures/PopID_distances.pdf",
       width = 8,
       height = 3,
       pp)


p_speciesDist_jac <- ggplot(mtrxs[["jaccard"]], aes(x=Species_dist,y=value)) +
    geom_boxplot() +
    theme_bw() +
    ylab("Distance") +
    xlab("Species distance") +
    ggtitle("Jaccard")+
    theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1))
p_speciesDist_bray <- ggplot(mtrxs[["bray"]], aes(x=Species_dist,y=value)) +
    geom_boxplot() +
    theme_bw() +
    ylab("Distance") +
    xlab("Species distance") +
    ggtitle("Bray-Curtis") +
    theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1))
pp <- plot_grid(p_speciesDist_jac, p_speciesDist_bray, ncol = 2)
pp

ggsave(file ="figures/Species_distancesSingular.pdf",
       width = 8,
       height = 3,
       pp)
