# Porthmeus
# 07.07.21

# The idea here is to show the significant difference in the alpha diversity by environmental factors corrected by the other fixed effects of the model


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
map_nut <- c("Hypereutr.", "Eutr.","Mesoeutr.")
meta[["Nutrient_load"]] <- factor(map_nut[metaLakes[meta$PopID,"NutrientLoad"]], levels= map_nut)
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

# create phyloseq object
phy <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
                tax_table(as.matrix(tax)),
                sample_data(meta))

alpha <- estimate_richness(phy,measures = "Shannon")
alpha <- merge(alpha, meta, by=0)

alpha <- data.table(alpha)
alpha[, Nutrient_load := factor(Nutrient_load, levels = map_nut)]

# correct the data for eutrophy and water body
mod <- lmer(data=alpha, Shannon ~ Nutrient_load + Water_body + (1|Species)+(1|Species:Reproduction))
summary(mod)
Anova(mod, type=2)
lsmeans(mod, pairwise~Water_body)
lsmeans(mod, pairwise~Nutrient_load)
plot(mod)
qqPlot(resid(mod))
hist(as.data.frame(cbind(resid(mod),alpha[,2])))

fixMod <- fixef(mod)
alpha[,ShannonWaterBody := Shannon]
alpha[,ShannonNutrient := Shannon]

for(n in names(fixMod)){
    if(n %in% grep("Nutrient_load", names(fixMod), value=TRUE)){
        n_n <- gsub("Nutrient_load","",n)
        alpha[Nutrient_load == n_n,ShannonWaterBody := Shannon * fixMod[n]]
    } else if (n %in% grep("Water_body", names(fixMod), value=TRUE)){
        n_n <- gsub("Water_body","",n)
        alpha[Water_body == n_n,ShannonNutrient:= Shannon * fixMod[n]]
    }
}



# correct data for Species and reproduction
mod <- lmer(data = alpha, Shannon ~ Species+Reproduction+(1|SiteID))
summary(mod)
Anova(mod, type=2)
lsmeans(mod, pairwise~Species)
lsmeans(mod, pairwise~Reproduction)
plot(mod)
qqPlot(resid(mod))
hist(as.data.frame(cbind(resid(mod),alpha[,2])))

fixMod <- fixef(mod)
alpha[,ShannonSpecies := Shannon]
alpha[,ShannonReproduction := Shannon]

for(n in names(fixMod)){
    if(n %in% grep("Species", names(fixMod), value=TRUE)){
        n_n <- gsub("Species","",n)
        alpha[Species == n_n,ShannonReproduction := Shannon * fixMod[n]]
    } else if (n %in% grep("Reproduction", names(fixMod), value=TRUE)){
        n_n <- gsub("Reproduction","",n)
        alpha[Reproduction == n_n,ShannonSpecies:= Shannon * fixMod[n]]
    }
}

# plot the corrected effects
boxWaterBody <- ggplot(alpha, aes(x = Water_body, y = ShannonWaterBody, fill = Water_body)) +
    geom_boxplot() +
    theme_bw() +
    labs(fill = "Water body", x = "Water body", y="Shannon index corrected\nby effect of nutrient load")+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
          legend.position = "None")
boxNutrient <- ggplot(alpha, aes(x = Nutrient_load, y = ShannonNutrient, fill = Nutrient_load)) +
    geom_boxplot() +
    theme_bw() +
    labs(fill = "Nutrient load", x = "Nutrient load", y="Shannon index corrected\nby effect of water body")+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
          legend.position = "None")

boxSpecies <- ggplot(alpha, aes(x = Species, y = ShannonSpecies, fill = Species)) +
    geom_boxplot() +
    theme_bw() +
    labs(fill = "Species", x = "Species", y="Shannon index corrected\nby effect of reproduction mode")+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
          legend.position = "None")
boxReproduction <- ggplot(alpha, aes(x = Reproduction, y = ShannonReproduction, fill = Reproduction)) +
    geom_boxplot() +
    theme_bw() +
    labs(fill = "Reproduction mode", x = "Reproduction mode", y="Shannon index corrected\nby effect of species")+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
          legend.position = "None")

pp <- plot_grid(boxWaterBody, boxNutrient,boxSpecies, boxReproduction, ncol = 2, labels = "AUTO")
ggsave(pp, file ="figures/FigS2_alphaCorByCofactor.pdf", width = 4, height = 6.4)
