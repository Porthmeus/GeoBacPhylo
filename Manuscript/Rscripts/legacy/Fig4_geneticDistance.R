# Porthmeus
# 15.02.21

# test for association between genetic distance and beta diversity
require(ggplot2)
require(phyloseq)
require(qiime2R)
require(DESeq2)
require(data.table)
require(cowplot)
require(Hmisc)
require(vegan)


otu <- read.csv("../data/dada2_AbFilt_FT.csv", row.names = 1)
tax <- read.csv("../data/dada2_AbFilt_Tax.csv", row.names =1)
tree <- read_qza("../data/dada2_AbFilt_fasta_tree.qza")$data

meta <- read.csv("../data/dada2_AbFilt_meta.csv", row.names=1)
map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]
meta[["Sample"]] <- rownames(meta)


otu2 <- otu+1

otu_vst <- varianceStabilizingTransformation(as.matrix(otu2), fitType = "local")
if(min(otu_vst) <0){
        otu_vst <- otu_vst + min(otu_vst)*-1
}
phy_vst <- phyloseq(otu_table(otu_vst, taxa_are_rows = TRUE),
                    tax_table(as.matrix(tax)),
                    sample_data(meta),
                    phy_tree(tree))


rdat <- load("Rscripts/genetic_distance2019.Rdata")
genDist <- get(rdat)
meta2 <- meta
rownames(meta2) <- meta2[["Isolation_ID"]]
rownames(genDist) <- meta2[rownames(genDist), "Sample"]
colnames(genDist) <- meta2[colnames(genDist), "Sample"]


distMs <- c( "bray","jaccard","wunifrac")
mtrxs <- list()
plots <- list()
stats <- list()
for(meas in distMs){
    mtrxs[[meas]] <- as.matrix(phyloseq::distance(phy_vst, method = meas))
    stts <- mantel(mtrxs[[meas]][rownames(genDist),colnames(genDist)],
                   genDist,
                   strata = meta[rownames(genDist),"PopID"],
                   method = "pearson")
    stats[[meas]] <- stts
    frm <- data.frame(
                      betaDiv = as.vector(as.dist(mtrxs[[meas]][rownames(genDist),colnames(genDist)])),
                      genDist = as.vector(as.dist(genDist)))
    if(meas == "bray"){
        measT <- "Bray-Curtis"
    } else {
        measT <- capitalize(meas)
    }
    plots[[meas]] <- ggplot(frm, aes(x = betaDiv, y = genDist))+
        geom_point() +
        geom_smooth(method="lm")+
        theme_bw() +
        ggtitle(measT) + 
        xlab("Beta-Diversity") +
        ylab("Genetic distance") + 
        geom_text(aes(x=min(betaDiv)+((max(betaDiv) -min(betaDiv))/100*5),
                      y=max(genDist), 
                      label = paste0("r = ",
                                     round(stats[[meas]][["statistic"]], digits = 2),
                                     "\np = ",
                                     round(stats[[meas]][["signif"]], digits = 3)
                                     )
                      ))
}

pp <- plot_grid(plotlist=plots, ncol = 3)
ggsave(pp, file="figures/GeneticDistance_BetaDiv.pdf", width = 10, height=5)
