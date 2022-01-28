# Porthmeus
# 29.08.21

# create a correlation matrix for all the relevant factors which are tested statistically

require(data.table)
require(mltools)
require(ggplot2)
require(cowplot)

AbFilt_meta = "../data/dada2_AbFilt_meta.csv"
MetaLakes = "../MetaDataLakes.csv"
meta <- read.csv(AbFilt_meta, row.names=1)
metaLakes <- read.csv(MetaLakes,row.names=1)

# simplify the reproduction vector
map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]

# add the correct naming for the nutritional state of the water bodies
map_nut <- c("Hypereutr.", "Eutr.","Mesoeutr.")
map_nut <- factor(map_nut, levels= map_nut)
meta[["Nutrient_load"]] <- map_nut[metaLakes[meta$PopID,"NutrientLoad"]]
meta[["Water_body2"]] <- metaLakes[meta$PopID, "Waterbody2"]
meta[["Water_body"]] <- metaLakes[meta$PopID, "waterbody"]
meta[["Water_body"]] <- gsub("River","flowing",gsub("Lake","standing",meta[["Water_body"]]))
meta[["LocationIDshort"]] <- sapply(meta$LocationID, function(x) strsplit(x, split ="/")[[1]][2])

# reorder the sampling sites, so that they occur in a similar order as no the map
siteOrder <- c("M34","M72","M71","M70","M31","M67","M89","M90","M86","M85","M84","M79","M78","M109","M28","M26","M44","M108","M107","M83","R12")
meta[["SiteID"]] <- factor(meta[["PopID"]], levels = siteOrder)
meta[["Nutrient_load"]] <- factor(meta[["Nutrient_load"]], levels = rev(map_nut))
meta[["SampleID"]] <- rownames(meta)

meta_s <- meta[,c("SiteID","Water_body","Nutrient_load","Species","Reproduction")]
meta_lst <-list()
for(cl in colnames(meta_s)){
    meta_lst[[cl]] <- as.factor(meta_s[,cl])
}
meta_oneHot <- one_hot(as.data.table(meta_lst))
meta_cor <- cor(as.matrix(meta_oneHot))

meta_s <- meta[,c("Water_body","Nutrient_load","Species","Reproduction")]
meta_lst <-list()
for(cl in colnames(meta_s)){
    meta_lst[[cl]] <- as.factor(meta_s[,cl])
}
meta_oneHot <- one_hot(as.data.table(meta_lst))
meta_cor <- cor(as.matrix(meta_oneHot), method = "kendall")
meta_corM <- reshape2::melt(meta_cor)
p_cor1 <- ggplot(meta_corM, aes(x=Var1,y=Var2, fill =value ))+
                 geom_tile() +
                 geom_text(aes(label = round(value,2)))+
                 labs(y="",x="", title = "Correlation matrix, Kendall", fill = "Tau") +
                 theme_bw() +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
                 scale_fill_distiller(palette = "RdBu")
p_cor1

meta_s <- meta[,c("Water_body","Nutrient_load","Species","Reproduction")]
meta_lst <-list()
for(cl in colnames(meta_s)){
    meta_lst[[cl]] <- as.integer(as.factor(meta_s[,cl]))
}

meta_cor2 <- cor(as.data.frame(meta_lst), method = "kendall")
meta_cor2M <- reshape2::melt(meta_cor2)
p_cor2 <- ggplot(meta_cor2M, aes(x=Var1,y=Var2, fill =value ))+
                 geom_tile() +
                 geom_text(aes(label = round(value,2), size = 0.5))+
                 labs(y="",x="", title = "Correlation matrix, Kendall", fill = "Tau") +
                 theme_bw() +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
                 scale_fill_distiller(palette = "RdBu")
p_both <- plot_grid(p_cor2, p_cor1,ncol = 1, labels = "AUTO")
p_both

ggsave(p_both, file = "figures/FigSX_correlationMatrix.pdf", width = 8, height = 8)


meta <- data.table(meta)
watNut <- meta[,.(Nutrient_load=Nutrient_load, Water_body=Water_body), by = PopID]
watNut <- as.matrix(watNut[!duplicated(PopID),2:3])

meta_lst <-list()
for(cl in colnames(watNut)){
    meta_lst[[cl]] <- as.integer(as.factor(watNut[,cl]))
}
meta_cor2 <- cor(as.data.frame(meta_lst), method = "kendall")
meta_cor2M <- reshape2::melt(meta_cor2)
p_cor2 <- ggplot(meta_cor2M, aes(x=Var1,y=Var2, fill =value ))+
                 geom_tile() +
                 geom_text(aes(label = round(value,2), size = 0.5))+
                 labs(y="",x="", title = "Correlation matrix, Kendall", fill = "Tau") +
                 theme_bw() +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
                 scale_fill_distiller(palette = "RdBu")
p_cor2

meta_lst <-list()
for(cl in colnames(watNut)){
    meta_lst[[cl]] <- as.factor(watNut[,cl])
}
meta_oneHot <- one_hot(as.data.table(meta_lst))
meta_cor <- cor(as.matrix(meta_oneHot))
meta_corM <- reshape2::melt(meta_cor)
p_cor2 <- ggplot(meta_corM, aes(x=Var1,y=Var2, fill =value ))+
                 geom_tile() +
                 geom_text(aes(label = round(value,2), size = 0.5))+
                 labs(y="",x="", title = "Correlation matrix, Kendall", fill = "Tau") +
                 theme_bw() +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
                 scale_fill_distiller(palette = "RdBu")
p_cor2
