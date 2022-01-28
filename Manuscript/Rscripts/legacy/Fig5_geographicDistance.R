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
require(geosphere)


# load and transform data
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

lakes <- read.csv("../MetaDataLakes.csv")
rownames(lakes) <- lakes[["SiteID"]]
meta2 <- merge(meta, lakes, by.x = "PopID", by.y = "SiteID", all.x= TRUE)
map_nut <- c("Hypertrophic", "Eutrophic","Mesotrophic")
meta2[["NutrientLoad"]] <-factor(map_nut[lakes[meta2[["PopID"]],"NutrientLoad"]],
                                levels = map_nut)


# calculate distance matrices, do stats and plot
distMs <- c("jaccard","bray")
mtrxs <- list()
plots <- list()
plots_singular <- list()
cmbs <- combn(meta2[,"Sample"], 2)
cmbs <- apply(cmbs,2,paste, collapse ="")
for(meas in distMs){
    mtrxs[[meas]] <- as.matrix(phyloseq::distance(phy_vst, method = meas))
    stts_adonis <- adonis(formula = mtrxs[[meas]] ~ factor(meta2$NutrientLoad))
    if(meas == "bray"){
        ttl <- "Bray-Curtis"
    } else {
        ttl <- capitalize(meas)
    }
    dst <- mtrxs[[meas]]
    dst <- reshape2::melt(dst, varnames= c("sample.x","sample.y"))
    dst <- data.table(dst)
    dst <- dst[paste0(sample.x,sample.y) %in% cmbs,]
    dst <- merge(dst, data.table(meta2), by.x = "sample.x", by.y="Sample")
    dst <- merge(dst, data.table(meta2), by.x = "sample.y", by.y="Sample", suffixes= c(".x",".y"))
    dst <- cbind(dst,
         dst[,.(Nutrition_inner = as.character(NutrientLoad.x == NutrientLoad.y))])
    dst[Nutrition_inner == "TRUE", Nutrition_inner := "within"]
    dst[Nutrition_inner == "FALSE", Nutrition_inner := "between"]
    dst[, Nutrition_dist := "Site_dist"]
    dst[Nutrition_inner == "within" & NutrientLoad.x == "Hypertrophic" & NutrientLoad.y == "Hypertrophic", Nutrition_dist := "Hyper-Hyper"]
    dst[Nutrition_inner == "within" & NutrientLoad.x == "Mesotrophic" & NutrientLoad.y == "Mesotrophic", Nutrition_dist := "Meso-Meso"]
    dst[Nutrition_inner == "within" & NutrientLoad.x == "Eutrophic" & NutrientLoad.y == "Eutrophic", Nutrition_dist := "Eu-Eu"]
    dst[Nutrition_inner == "between" & NutrientLoad.x == "Mesotrophic" & NutrientLoad.y == "Eutrophic", Nutrition_dist := "Eu-Meso"]
    dst[Nutrition_inner == "between" & NutrientLoad.y == "Mesotrophic" & NutrientLoad.x == "Eutrophic", Nutrition_dist := "Eu-Meso"]
    dst[Nutrition_inner == "between" & NutrientLoad.x == "Hypertrophic" & NutrientLoad.y == "Eutrophic", Nutrition_dist := "Hyper-Eu"]
    dst[Nutrition_inner == "between" & NutrientLoad.y == "Hypertrophic" & NutrientLoad.x == "Eutrophic", Nutrition_dist := "Hyper-Eu"]
    dst[Nutrition_inner == "between" & NutrientLoad.x == "Hypertrophic" & NutrientLoad.y == "Mesotrophic", Nutrition_dist := "Hyper-Meso"]
    dst[Nutrition_inner == "between" & NutrientLoad.y == "Hypertrophic" & NutrientLoad.x == "Mesotrophic", Nutrition_dist := "Hyper-Meso"]
    dst[,Nutrition_dist := factor(Nutrition_dist, levels = c("Hyper-Hyper","Eu-Eu","Meso-Meso","Hyper-Eu","Eu-Meso","Hyper-Meso"))]
#
    pp <- ggplot(dst, aes(x = Nutrition_inner, y = value)) +
        geom_boxplot() +
        theme_bw() +
        ggtitle(ttl) +
        ylab("Distance") +
        xlab("Nutrition level of water body") 
    plots[[meas]] <- pp
    #
    pp <- ggplot(dst, aes(x = Nutrition_dist, y = value)) +
        geom_boxplot() +
        theme_bw() +
        ggtitle(ttl) +
        ylab("Distance") +
        xlab("Nutrition level of water body") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    plots_singular[[meas]] <- pp
}

p_dist <- plot_grid(plotlist=plots, ncol = 2)
ggsave(p_dist, file ="figures/Distance_NutritionLakes.pdf", height =4, width =6)
p_dist <- plot_grid(plotlist=plots_singular, ncol = 2)
ggsave(p_dist, file ="figures/Distance_NutritionLakesSingular.pdf", height =4, width =8)



# geographic and microbial distance
geoDist <- distm(lakes[,c("lon","lat")])
colnames(geoDist) <- lakes[["SiteID"]]
rownames(geoDist) <- lakes[["SiteID"]]
#geoDist[upper.tri(geoDist, diag=FALSE)] <- NA

# create a data frame which associates distance of microbiome to geographic distance
#geoDistLong <- data.table(reshape2::melt(geoDist))
#geoDistLong <- geoDistLong[!is.na(value),]

mtrxsLong <- list()
for(mat in distMs){
    bacDist <- as.matrix(mtrxs[[mat]])
    bacDist[upper.tri(bacDist, diag=TRUE)] <- NA
    bacDistLong <- data.table(reshape2::melt(bacDist))
    bacDistLong <- bacDistLong[!is.na(value),]
    # add the PopID information to associate it back to the geographic distance
    bacDistLong <- cbind(bacDistLong, 
                          Var1PopID = meta[bacDistLong[,Var1], "PopID"],
                          Var2PopID = meta[bacDistLong[,Var2], "PopID"],
                          GeoDist = -1)
    bacDistLong <- data.frame(bacDistLong)
    for(i in 1:nrow(bacDistLong)){
        bacDistLong[i,"GeoDist"] <- geoDist[bacDistLong[i,"Var1PopID"],
                                          bacDistLong[i,"Var2PopID"]]
    }
    bacDistLong <- data.table(bacDistLong)
    mtrxsLong[[mat]] <- bacDistLong
}


# create a distance matrix for the distance in geo-location of the same dimensions like the beta-diversity to perform a mantel test
# 1. cast into wide format
geoDist_Samples <- dcast(mtrxsLong[[1]], Var1~Var2, value.var = "GeoDist")
# 2. format as data.frame, correct rownames, which are missing in data.table
geoDist_Samples_names <- geoDist_Samples[, Var1]
geoDist_Samples <- as.data.frame(geoDist_Samples)
rownames(geoDist_Samples) <- geoDist_Samples_names
geoDist_Samples <- geoDist_Samples[,-1]
# 3. add an additional column and row, because I left out the diagonal by the first conversion from distance to long format
geoDist_Samples <- cbind(geoDist_Samples, added = NA)
colnames(geoDist_Samples)[ncol(geoDist_Samples)] <- geoDist_Samples_names[length(geoDist_Samples_names)]
geoDist_Samples_names <- colnames(geoDist_Samples)
geoDist_Samples <- rbind(added = NA, geoDist_Samples)
rownames(geoDist_Samples)[1] <- geoDist_Samples_names[1]
# 4. reformat into distance object
geoDist_Samples2 <- geoDist_Samples
geoDist_Samples2[geoDist_Samples2==0] <- NA
geoDist_Samples <- as.dist(geoDist_Samples)
geoDist_Samples2 <- as.dist(geoDist_Samples2)



geoDistPlots <- list()
geoDistModels <- list()
for(mat in names(mtrxsLong)){
    stts <- mantel(geoDist_Samples, mtrxs[[mat]])
    if(mat == "bray"){
        nm <- "Bray-Curtis"
    } else {
        nm <- capitalize(mat)
    }
    p <- ggplot(mtrxsLong[[mat]], aes(y=value, x=GeoDist)) +
        geom_point(alpha = 0.1) +
        geom_text(aes(y=min(value)+ ((max(value)-min(value))*0.05),
                      x = min(GeoDist) + ((max(geoDist)-min(geoDist))*0.95),
                      label = paste0("r = ",
                                     round(stts[["statistic"]], digits=2),
                                     "\np = ",
                                     round(stts[["signif"]], digits=3)))) +
        xlab("Geographic distance") +
        ylab("Beta diversity") +
        ggtitle(nm) +
        theme_bw()
    geoDistPlots[[mat]] <- p
}

pp <- cowplot::plot_grid(plotlist = geoDistPlots, ncol = 1)
ggsave(pp, file="figures/GeographicDistance_BetaDiv.pdf", width = 6, height=6)
