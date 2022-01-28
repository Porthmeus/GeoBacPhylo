# Porthmeus
# 20.05.21

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
for(fac in c("SiteID","Water_body","Nutrient_load","Species","Reproduction")){
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

saveRDS(mods, file = "tables/AlphaDiv_Simple_models.RDS")

anovas <- do.call(rbind, anovas)
anovas[["padj"]] <- p.adjust(anovas[["Pr(>F)"]], method = "BH")
write.csv(file="tables/AlphaDiv_Simple_AnovaTab_typ1.csv", anovas, row.names = FALSE)

tukeys <- do.call(rbind, tukeys)
write.csv(file="tables/AlphaDiv_Simple_TukeyTab.csv", tukeys, row.names = FALSE)





# All in one
mod <- lmer(data = alpha,
            Shannon ~ Water_body + Nutrient_load + Species + Reproduction + (1|PopID)
            )
summary(mod)
Anova(mod,type =2)
pairwiseTab_Species <- lsmeans(mod, pairwise~Species)
pairwiseTab_Reproduction <- lsmeans(mod, pairwise~Reproduction)
pairwiseTab_water <- lsmeans(mod, pairwise~Water_body)
pairwiseTab_nutrient <- lsmeans(mod, pairwise~Nutrient_load)
#plot(mod)
#qqPlot(resid(mod))
#hist(as.data.frame(cbind(resid(mod),alpha[,2])))
saveRDS(mod, file = "tables/AlphaDiv_AllInOne_lmMod.RDS")
write.csv(file="tables/AlphaDiv_AllInOne_AnovaTab_typ2.csv", as.data.frame(Anova(mod,type=2)))

pairwiseTab_AllInOne <- rbind(
                              as.data.frame(pairwiseTab_water$contrasts),
                              as.data.frame(pairwiseTab_nutrient$contrasts),
                              as.data.frame(pairwiseTab_Species$contrasts),
                              as.data.frame(pairwiseTab_Reproduction$contrasts)
)

pairwiseTab_AllInOne[["padj"]] <- p.adjust(pairwiseTab_AllInOne[,"p.value"], method = "BH")
write.csv(file="tables/AlphaDiv_AllInOne_TukeyPairwiseTab.csv", pairwiseTab_AllInOne ,row.names=FALSE)

# older stuff I have done, but which was too complicated

#mod <- lm(data=alpha, Shannon ~  Nutrient_load + Water_body )
#aov_mod <- Anova(mod,type=3)
#aov_mod 
#aov_mod[["Sum Sq"]] /sum(aov_mod[["Sum Sq"]])
#summary(mod)
#
#mod <- lmer(data=alpha, Shannon ~ Water_body+Nutrient_load+ (1|SiteID))
#summary(mod)
#Anova(mod,type=2)
#
#mod <- lm(data=alpha, Shannon ~ Water_body+Nutrient_load)
#summary(mod)
#Anova(mod,type=2)
#
#mod <- lmer(data=alpha, Shannon ~ Species+Reproduction+(1|SiteID)+(1|SiteID:LocationID))
#summary(mod)
#Anova(mod, type=3)
#plot(mod)
#qqPlot(resid(mod))
#hist(as.data.frame(cbind(resid(mod),alpha[,2])))

#mod <- lmer(data=alpha, Shannon ~ SiteID + (1|Species)+(1|Species:Reproduction))
#summary(mod)
#Anova(mod, type=2)
#pairwiseTab <- lsmeans(mod, pairwise~SiteID)
#print(pairwiseTab)
#plot(mod)
#qqPlot(resid(mod))
#hist(as.data.frame(cbind(resid(mod),alpha[,2])))
#saveRDS(mod, file = "tables/AlphaDiv_SiteID_lmerMod.RDS")
#write.csv(file="tables/AlphaDiv_SiteID_AnovaTab.csv", as.data.frame(Anova(mod,type=2)))
#write.csv(file="tables/AlphaDiv_SiteID_TukeyPairwiseTab.csv", as.data.frame(pairwiseTab$contrasts),row.names=FALSE)
#
#mod <- lmer(data=alpha, Shannon ~ Nutrient_load + Water_body + (1|Species)+(1|Species:Reproduction))
#summary(mod)
#Anova(mod, type=2)
#pairwiseTab_water <- lsmeans(mod, pairwise~Water_body)
#pairwiseTab_nutrient <- lsmeans(mod, pairwise~Nutrient_load)
#plot(mod)
#qqPlot(resid(mod))
#hist(as.data.frame(cbind(resid(mod),alpha[,2])))
#saveRDS(mod, file = "tables/AlphaDiv_WaterBodyNutrientLoad_lmerMod.RDS")
#write.csv(file="tables/AlphaDiv_WaterBodyNutrientLoad_AnovaTab.csv", as.data.frame(Anova(mod,type=2)))
#write.csv(file="tables/AlphaDiv_NutrientLoad_TukeyPairwiseTab.csv", as.data.frame(pairwiseTab_nutrient$contrasts),row.names=FALSE)
#write.csv(file="tables/AlphaDiv_WaterBody_TukeyPairwiseTab.csv", as.data.frame(pairwiseTab_water$contrasts))
#
#mod <- lmer(data=alpha, Shannon ~ Nutrient_load + Water_body + (1|SiteID))
#summary(mod)
#Anova(mod, type=2)
#pairwiseTab_water <- lsmeans(mod, pairwise~Water_body)
#pairwiseTab_nutrient <- lsmeans(mod, pairwise~Nutrient_load)
#plot(mod)
#qqPlot(resid(mod))
#hist(as.data.frame(cbind(resid(mod),alpha[,2])))
#saveRDS(mod, file = "tables/AlphaDiv_WaterBodyNutrientLoad_RandomSiteID_lmerMod.RDS")
#write.csv(file="tables/AlphaDiv_WaterBodyNutrientLoad_RandomSiteID_AnovaTab.csv", as.data.frame(Anova(mod,type=2)))
#write.csv(file="tables/AlphaDiv_NutrientLoad_RandomSiteID_TukeyPairwiseTab.csv", as.data.frame(pairwiseTab_nutrient$contrasts),row.names=FALSE)
#write.csv(file="tables/AlphaDiv_WaterBody_RandomSiteID_TukeyPairwiseTab.csv", as.data.frame(pairwiseTab_water$contrasts))
#
#
#mod <- lmer(data = alpha, Shannon ~ Species+Reproduction+(1|SiteID))
#summary(mod)
#Anova(mod,type =2)
#pairwiseTab_Species <- lsmeans(mod, pairwise~Species)
#pairwiseTab_Reproduction <- lsmeans(mod, pairwise~Reproduction)
#print(pairwiseTab)
#plot(mod)
#qqPlot(resid(mod))
#hist(as.data.frame(cbind(resid(mod),alpha[,2])))
#saveRDS(mod, file = "tables/AlphaDiv_SpeciesReproduction_lmMod.RDS")
#write.csv(file="tables/AlphaDiv_SpeciesReproduction_AnovaTab_typ2.csv", as.data.frame(Anova(mod,type=2)))
#write.csv(file="tables/AlphaDiv_Species_TukeyPairwiseTab.csv", as.data.frame(pairwiseTab_Species$contrasts),row.names=FALSE)
#write.csv(file="tables/AlphaDiv_Reproduction_TukeyPairwiseTab.csv", as.data.frame(pairwiseTab_Reproduction$contrasts),row.names=FALSE)
#

#mod <- lm(data = alpha, Shannon ~ SiteID)
#anovaTab <- anova(mod)
#tukeyTab <- TukeyHSD(aov(mod))[[1]]
#write.csv(anovaTab,file = "tables/AlphaDivSimple_SiteID_AnovaTab.csv")
#write.csv(tukeyTab,file = "tables/AlphaDivSimple_SiteID_TukeyTab.csv")
#mod <- lm(data = alpha, Shannon ~ Species)
#anovaTab <- anova(mod)
#tukeyTab <- TukeyHSD(aov(mod))[[1]]
#write.csv(anovaTab,file = "tables/AlphaDivSimple_Species_AnovaTab.csv")
#write.csv(tukeyTab,file = "tables/AlphaDivSimple_Species_TukeyTab.csv")
#mod <- lm(data = alpha, Shannon ~ Reproduction)
#anovaTab <- anova(mod)
#tukeyTab <- TukeyHSD(aov(mod))[[1]]
#write.csv(anovaTab,file = "tables/AlphaDivSimple_Reproduction_AnovaTab.csv")
#write.csv(tukeyTab,file = "tables/AlphaDivSimple_Reproduction_TukeyTab.csv")
#mod <- lm(data = alpha, Shannon ~ Water_body)
#anovaTab <- anova(mod)
#tukeyTab <- TukeyHSD(aov(mod))[[1]]
#write.csv(anovaTab,file = "tables/AlphaDivSimple_WaterBody_AnovaTab.csv")
#write.csv(tukeyTab,file = "tables/AlphaDivSimple_WaterBody_TukeyTab.csv")
#mod <- lm(data = alpha, Shannon ~ Nutrient_load)
#anovaTab <- anova(mod)
#tukeyTab <- TukeyHSD(aov(mod))[[1]]
#write.csv(anovaTab,file = "tables/AlphaDivSimple_NutrientLoad_AnovaTab.csv")
#write.csv(tukeyTab,file = "tables/AlphaDivSimple_NutrientLoad_TukeyTab.csv")
