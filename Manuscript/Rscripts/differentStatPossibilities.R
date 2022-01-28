# Porthmeus
# 14.05.21

# here is some code which can be considered when dealing with the statistical problems of microbial association in the Hydra samples


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
meta[["Nutrient_load"]] <- map_nut[metaLakes[meta$PopID,"NutrientLoad"]]
meta[["Water_body2"]] <- metaLakes[meta$PopID, "Waterbody2"]
meta[["Water_body"]] <- metaLakes[meta$PopID, "waterbody"]
meta[["Water_body"]] <- gsub("River","flowing",gsub("Lake","standing",meta[["Water_body"]]))
meta[["LocationIDshort"]] <- sapply(meta$LocationID, function(x) strsplit(x, split ="/")[[1]][2])

# remove rare species (lower 25% of the data)
countsPerESV <- apply(otu,1,sum)
thrs <- quantile(countsPerESV,0.25)
summary(countsPerESV)/sum(otu)
plot(ecdf(log10(countsPerESV+1)))
abline(v = log10(thrs+1))
hist(log10(countsPerESV))
abline(v=log10(thrs))

ESVs <- names(countsPerESV)[countsPerESV > thrs]
otu <- otu[ESVs,]
tax <- tax[ESVs,]

# reorder the sampling sites, so that they occur in a similar order as no the map
siteOrder <- c("M34","M72","M71","M70","M31","M67","M89","M90","M86","M85","M84","M79","M78","M109","M28","M26","M44","M108","M107","M83","R12")
meta[["SiteID"]] <- factor(meta[["PopID"]], levels = siteOrder)
meta[["Nutrient_load"]] <- factor(meta[["Nutrient_load"]], levels = rev(map_nut))

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
smpls <- rownames(as.matrix(brayCurt))
meta <- meta[smpls,]

# these are the final stats!

# PERMANOVA is sensitive to the order of terms and basically only the last term is what is really the effect of the factor after all other variability has been explained by the other factors - thus I permute the order to get the "real" p-values and sum of squares
perm <- how(nperm  = 9999, blocks=meta$SiteID)
adn <- adonis2(data = meta, formula = brayCurt ~ Water_body+Nutrient_load+Species+Reproduction, permutations = perm, method = "Terms")
res <- as.data.frame(adn)[4,]
adn <- adonis2(data = meta, formula = brayCurt ~ Water_body+Reproduction+Nutrient_load+Species, permutations = perm, method = "Terms")
res <- rbind(res, as.data.frame(adn)[4,])
adn <- adonis2(data = meta, formula = brayCurt ~ Water_body+Species+Reproduction+Nutrient_load, permutations = perm, method = "Terms")
res <- rbind(res, as.data.frame(adn)[4,])
adn <- adonis2(data = meta, formula = brayCurt ~ Nutrient_load+Species+Reproduction+Water_body, permutations = perm, method = "Terms")
res <- rbind(res, as.data.frame(adn)[4:6,])
res["Residual","SumOfSqs"] <- res["Total","SumOfSqs"] - sum(res[1:4,"SumOfSqs"])
res["Residual","R2"] <- res["Total","R2"] - sum(res[1:4,"R2"])
write.csv(res, file ="tables/BetaDiv_AllInOne_adonisTab.csv")

# for technical reasons we decided to test all factors individually
adn <- adonis2(data = meta, formula = brayCurt ~ SiteID, permutations = 9999, method = "Terms")
res <- cbind(TestedModel = "SiteID",Stat = rownames(adn), adn)
adn <- adonis2(data = meta, formula = brayCurt ~ Water_body, permutations = 9999, method = "Terms")
res <- rbind(res,cbind(TestedModel = "Water_body",Stat = rownames(adn), adn))
adn <- adonis2(data = meta, formula = brayCurt ~ Nutrient_load, permutations = 9999, method = "Terms")
res <- rbind(res,cbind(TestedModel = "Nutrient_load",Stat = rownames(adn), adn))
adn <- adonis2(data = meta, formula = brayCurt ~ Species, permutations = 9999, method = "Terms")
res <- rbind(res,cbind(TestedModel = "Species",Stat = rownames(adn), adn))
adn <- adonis2(data = meta, formula = brayCurt ~ Reproduction, permutations = 9999, method = "Terms")
res <- rbind(res,cbind(TestedModel = "Reproduction",Stat = rownames(adn), adn))
write.csv(res, file ="tables/BetaDiv_Individual_adonisTab.csv")

# old stuff which we tested at some point
#perm <- how(nperm  = 999, blocks=meta$Species, plots = Plots(strata=meta$Reproduction))
#adn <- adonis2(data = meta, formula = brayCurt ~ Water_body+Nutrient_load, permutations = perm, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDiv_NutrientLoad_adonisTab.csv")
#
#adn <- adonis2(data = meta, formula = brayCurt ~ Nutrient_load+Water_body, permutations = perm, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDiv_WaterBody_adonisTab.csv")
#
#
#adn <- adonis2(data = meta, formula = brayCurt ~ SiteID, permutations = perm, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDiv_SiteID_adonisTab.csv")
#
#perm <- how(nperm  = 999, blocks=meta$SiteID)
#adn <- adonis2(data = meta, formula = brayCurt ~ Species+Reproduction, permutations = perm, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDiv_Reproduction_adonisTab.csv")
#adn <- adonis2(data = meta, formula = brayCurt ~ Reproduction+Species, permutations = perm, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDiv_Species_adonisTab.csv")
#
## the last one gives amount variability explained to the different features
#adn <- adonis2(data = meta, formula = brayCurt ~ Species+Reproduction+Nutrient_load+Water_body+SiteID+LocationID, permutations =2, method = "Margin")
#write.csv(as.data.frame(adn), file = "tables/BetaDiv_VariationSplitAll_adonisTab.csv")
#
#
## pairwise tests
#meta2 <- data.table(meta, keep.rownames=TRUE)
#combs <- combn(unique(meta2$Nutrient_load),2)
#adnRes <- list()
#for(i in 1:ncol(combs)){
#    contrast <- paste(combs[,i], collapse = "-")
#    samples <- meta2[Nutrient_load %in% combs[,i],rn] 
#    perm <- how(nperm  = 999, blocks=meta2[Nutrient_load %in% combs[,i], Species], plots = Plots(strata=meta2[Nutrient_load %in% combs[,i],Reproduction]))
#    adn<- adonis2(data = meta2[Nutrient_load %in% combs[,i],], formula = as.dist(as.matrix(brayCurt)[samples,samples]) ~ Water_body+Nutrient_load, permutations = perm, method = "Margin")
#    adn <- as.data.frame(adn)
#    adn <- cbind(factor=rownames(adn), adn)
#    adn <- cbind(contrast = contrast, adn[1:(nrow(adn)-1),])
#    adnRes[[contrast]] <- adn
#}
#adnRes <- do.call(rbind, adnRes)
#adnRes[["padj"]] <- p.adjust(adnRes[,ncol(adnRes)], method = "BH")
#write.csv(adnRes, file = "tables/BetaDiv_NutrientLoad_PairwiseTab.csv", row.names=FALSE)
#
## pairwise tests
#meta2 <- data.table(meta, keep.rownames=TRUE)
#combs <- combn(unique(meta2$Species),2)
#adnRes <- list()
#for(i in 1:ncol(combs)){
#    contrast <- paste(combs[,i], collapse = "-")
#    samples <- meta2[Species %in% combs[,i],rn] 
#    perm <- how(nperm  = 999, blocks=meta2[Species %in% combs[,i],SiteID])
#    adn <- adonis2(data = meta2[Species %in% combs[,i],], formula = as.dist(as.matrix(brayCurt)[samples,samples]) ~ Reproduction+Species, permutations = perm, method = "Terms")
#    adn <- as.data.frame(adn)
#    adn <- cbind(factor=rownames(adn), adn)
#    adn <- cbind(contrast = contrast, adn[1:(nrow(adn)-1),])
#    adnRes[[contrast]] <- adn
#}
#adnRes <- do.call(rbind, adnRes)
#adnRes[["padj"]] <- p.adjust(adnRes[,ncol(adnRes)], method = "BH")
#write.csv(adnRes, file = "tables/BetaDiv_Species_PairwiseTab.csv", row.names=FALSE)
#
#meta2 <- data.table(meta, keep.rownames=TRUE)
#combs <- combn(unique(meta2$Reproduction),2)
#adnRes <- list()
#for(i in 1:ncol(combs)){
#    contrast <- paste(combs[,i], collapse = "-")
#    samples <- meta2[Reproduction %in% combs[,i],rn] 
#    perm <- how(nperm  = 999, blocks=meta2[Reproduction %in% combs[,i],SiteID])
#    adn <- adonis2(data = meta2[Reproduction %in% combs[,i],], formula = as.dist(as.matrix(brayCurt)[samples,samples]) ~ Species+Reproduction, permutations = perm, method = "Terms")
#    adn <- as.data.frame(adn)
#    adn <- cbind(factor=rownames(adn), adn)
#    adn <- cbind(contrast = contrast, adn[1:(nrow(adn)-1),])
#    adnRes[[contrast]] <- adn
#}
#adnRes <- do.call(rbind, adnRes)
#adnRes[["padj"]] <- p.adjust(adnRes[,ncol(adnRes)], method = "BH")
#adnRes
#write.csv(adnRes, file = "tables/BetaDiv_Reproduction_PairwiseTab.csv", row.names=FALSE)
#
## do the simple statistics
#adn <- adonis2(data = meta, formula = brayCurt ~ SiteID, permutations = 999, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDivSimple_SiteID_adonisTab.csv")
#adn <- adonis2(data = meta, formula = brayCurt ~ Water_body, permutations = 999, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDivSimple_Water_body_adonisTab.csv")
#adn <- adonis2(data = meta, formula = brayCurt ~ Nutrient_load, permutations = 999, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDivSimple_NutrientLoad_adonisTab.csv")
#adn <- adonis2(data = meta, formula = brayCurt ~ Species, permutations = 999, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDivSimple_Species_adonisTab.csv")
#adn <- adonis2(data = meta, formula = brayCurt ~ Reproduction, permutations = 999, method = "Terms")
#write.csv(as.data.frame(adn), file = "tables/BetaDivSimple_Reproduction_adonisTab.csv")
#
#
#
#perm <- how(nperm  = 999, blocks=meta$SiteID)
#adn <- adonis2(data = meta, formula = brayCurt ~ Species+Reproduction+Nutrient_load+Water_body, permutations = perm, method = "Terms")
#adonis2(data = meta, formula = brayCurt ~ Species+Reproduction+Water_body+Nutrient_load, permutations = perm, method = "Terms")
#adonis2(data = meta, formula = brayCurt ~ Species+Water_body+Nutrient_load+Reproduction, permutations = perm, method = "Terms")
#adonis2(data = meta, formula = brayCurt ~ Water_body+Nutrient_load+Reproduction+Species, permutations = perm, method = "Terms")
#
