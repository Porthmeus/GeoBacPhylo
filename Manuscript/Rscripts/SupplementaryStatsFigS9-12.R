# Porthmeus
# 16.12.21

# do the statistics as before only on the oli samples
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
meta[["Water_body"]] <- gsub("River","running",gsub("Lake","standing",meta[["Water_body"]]))
meta[["LocationIDshort"]] <- sapply(meta$LocationID, function(x) strsplit(x, split ="/")[[1]][2])

# reorder the sampling sites, so that they occur in a similar order as no the map
siteOrder <- c("M34","M72","M71","M70","M31","M67","M89","M90","M86","M85","M84","M79","M78","M109","M28","M26","M44","M108","M107","M83","R12")
meta[["SiteID"]] <- factor(meta[["PopID"]], levels = siteOrder)

# remove all non-Oli samples
meta <- meta[meta[["Species"]] == "OLI",]

# remove rare species (lower 25% of the data)
countsPerESV <- apply(otu,1,sum)
thrs <- quantile(countsPerESV,0.25)
summary(countsPerESV)/sum(otu)
plot(ecdf(log10(countsPerESV+1)))
abline(v = log10(thrs))
hist(log10(countsPerESV))
abline(v=log10(thrs))

ESVs <- names(countsPerESV)[countsPerESV > thrs]
otu <- otu[ESVs,]
tax <- tax[ESVs,]

# remove chloroplast ESVs from (probably) algea
sel <- which(tax[["Order"]] != "Chloroplast")
tax <- tax[sel,]
otu <- otu[rownames(tax),]

# keep just the Oli samples
otu <- otu[, colnames(otu) %in% rownames(meta)]

# create phyloseq object
phy <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
                tax_table(as.matrix(tax)),
                sample_data(meta))


alpha <- estimate_richness(phy,measures = "Shannon")
alpha <- merge(alpha, meta, by=0)
alpha[["Nutrient_load"]] <- relevel(as.factor(alpha[["Nutrient_load"]]), "Mesoeutr.")
 
# simple 
mods <- list()
anovas <- list()
tukeys <- list()
for(fac in c("SiteID","Water_body","Nutrient_load","Reproduction")){
    fom <- as.formula(paste0("Shannon ~ ",fac))
    mod <- lm(data = alpha, fom)
    mods[[fac]] <- mod
    anva <- as.data.frame(anova(mod))
    anovas[[fac]] <- cbind(TestedModel = fac,
                           Factor = rownames(anva),
                           anva)
    tuk <- as.data.frame(TukeyHSD(aov(mod))[[1]])
    tukeys[[fac]] <- cbind(TestedModel = fac,
                           Comparison = rownames(tuk),
                           tuk)
}

saveRDS(mods, file = "tables/AlphaDivSupplementFig9-12_Simple_models.RDS")

anovas <- do.call(rbind, anovas)
anovas[["padj"]] <- p.adjust(anovas[["Pr(>F)"]], method = "BH")
write.csv(file="tables/AlphaDivSupplementFig9-12_Simple_AnovaTab_typ1.csv", anovas, row.names = FALSE)

tukeys <- do.call(rbind, tukeys)
write.csv(file="tables/AlphaDivSupplementFig9-12_Simple_TukeyTab.csv", tukeys, row.names = FALSE)





# All in one
mod <- lmer(data = alpha,
            Shannon ~ Water_body + Nutrient_load + Reproduction + (1|PopID)
            )
summary(mod)
Anova(mod,type =2)
pairwiseTab_Reproduction <- lsmeans(mod, pairwise~Reproduction)
pairwiseTab_water <- lsmeans(mod, pairwise~Water_body)
pairwiseTab_nutrient <- lsmeans(mod, pairwise~Nutrient_load)
#plot(mod)
#qqPlot(resid(mod))
#hist(as.data.frame(cbind(resid(mod),alpha[,2])))
saveRDS(mod, file = "tables/AlphaDivSupplementFig9-12_AllInOne_lmMod.RDS")
write.csv(file="tables/AlphaDivSupplementFig9-12_AllInOne_AnovaTab_typ2.csv", as.data.frame(Anova(mod,type=2)))

pairwiseTab_AllInOne <- rbind(
                              as.data.frame(pairwiseTab_water$contrasts),
                              as.data.frame(pairwiseTab_nutrient$contrasts),
                              as.data.frame(pairwiseTab_Reproduction$contrasts)
)

pairwiseTab_AllInOne[["padj"]] <- p.adjust(pairwiseTab_AllInOne[,"p.value"], method = "BH")
write.csv(file="tables/AlphaDivSupplementFig9-12_AllInOne_TukeyPairwiseTab.csv", pairwiseTab_AllInOne ,row.names=FALSE)



# now Beta diversity

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
adn <- adonis2(data = meta, formula = brayCurt ~ Water_body+Nutrient_load+Reproduction, permutations = perm, method = "Terms")
res <- as.data.frame(adn)[3,]
adn <- adonis2(data = meta, formula = brayCurt ~ Water_body+Reproduction+Nutrient_load, permutations = perm, method = "Terms")
res <- rbind(res, as.data.frame(adn)[3,])
adn <- adonis2(data = meta, formula = brayCurt ~ Nutrient_load+Reproduction+Water_body, permutations = perm, method = "Terms")
res <- rbind(res, as.data.frame(adn)[3:5,])
res["Residual","SumOfSqs"] <- res["Total","SumOfSqs"] - sum(res[1:3,"SumOfSqs"])
res["Residual","R2"] <- res["Total","R2"] - sum(res[1:3,"R2"])
write.csv(res, file ="tables/BetaDivSupplementFig9-12_AllInOne_adonisTab.csv")

# for technical reasons we decided to test all factors individually
adn <- adonis2(data = meta, formula = brayCurt ~ SiteID, permutations = 9999, method = "Terms")
res <- cbind(TestedModel = "SiteID",Stat = rownames(adn), adn)
adn <- adonis2(data = meta, formula = brayCurt ~ Water_body, permutations = 9999, method = "Terms")
res <- rbind(res,cbind(TestedModel = "Water_body",Stat = rownames(adn), adn))
adn <- adonis2(data = meta, formula = brayCurt ~ Nutrient_load, permutations = 9999, method = "Terms")
res <- rbind(res,cbind(TestedModel = "Nutrient_load",Stat = rownames(adn), adn))
adn <- adonis2(data = meta, formula = brayCurt ~ Reproduction, permutations = 9999, method = "Terms")
res <- rbind(res,cbind(TestedModel = "Reproduction",Stat = rownames(adn), adn))
write.csv(res, file ="tables/BetaDivSupplementFig9-12_Individual_adonisTab.csv")
