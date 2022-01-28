# Porthmeus
# 25.05.20

require(data.table)
require(ggplot2)
require(biomformat)

biomFiles <- "../data/PE_cut_dada2_FT.qza"
statFiles <- "../reports/PE_cut_dada2_Stats.qza"

for(ft in biomFiles){
    # go through the biomeFiles and extract the number of ESVs
    ft_uz <- grep("feature-table.biom",unzip(ft,list=T)[,1], value = T)
    unzip(ft, ft_uz)
    ft_uz <- file.path(dirname(ft), ft_uz)
    ft_biom <- read_biom(ft_uz)
    # TODO add some function for combining the data frame of extract the data which might be important
}
# 

