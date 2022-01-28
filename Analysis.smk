# Porthmeus
# 22.05.20

# Analysis pipeline for Qiime and the environmental 16s sequences obtained by Jacint Tökölyi 

# define samples globally for snakemake
import pandas as pd
import os

df = pd.read_table("FileManifest.tab")
spl = list(df.iloc[:,1])
spl.extend(list(df.iloc[:,2]))
samples = [os.path.split(x)[1].strip(".fastq.gz") for x in spl]

# define the thresholds to evaluate how filtering affects the data
trunc_thr = list(range(150,315,15))
qualities = list(range(0,22,2))

# define denoising methods
dn = ["deblur","deblurSilva","dada2", "dada2NoFilt"]
rule all:
    input: 
        qiimeEnd = "data/PE_denoise_dada2_physeq.RDS",
        reports = ["Rmarkdown/Natrix_EvaluationAbundanceFilter.html",
            "Rmarkdown/Natrix_StatisticalAnalysis.html",
            "Rmarkdown/Natrix_IndicatorSpeciesDESeq2.html"],
        FT = "data/dada2_AbFilt_FT.csv",
        meta = "data/dada2_AbFilt_meta.csv",
        tax = "data/dada2_AbFilt_Tax.csv",
        fasta = "data/dada2_AbFilt_fasta.fasta"
'''
        ipt = expand(["data/PE_demux_trim.qza",
            "reports/multiqc_raw.html",
            "reports/PE_demux_summary.qzv",
            "reports/PE_demux_cut_summary.qzv",
            "reports/PE_demux_trim_summary.qzv",
            "reports/PE_denoise_deblur_Stats.qzv",
            "reports/PE_denoise_{dn}_seqTab.qzv",
            "reports/PE_demux_cut_joined.qzv",
            "reports/PE_denoise_{dn}_sumFT.qzv",
            "reports/{infile}_trim_summary.qzv",
            "reports/PE_denoise_{dn}_Stats.qzv",
            "results/distance/{dn}_bray_curtis_emperor.qzv",
            "reports/PE_denoise_{dn}_Tax.qzv",
            "reports/PE_denoise_{dn}_sumRS.qzv",
            "results/taxa_bar_{dn}.qzv",
            "reports/{dn}_ThrlOnAlphaDiv.csv",
            "reports/thresholdSurvivors_{dn}.csv",
            "plots/eval/thresholdSurvivors.svg",
            "plots/eval/ThrlOnAlphaDiv.svg",
            "plots/eval/ThrlOnBetaDiv.jpg",
            "data/dada2_AbFilt_fasta_tree.qza",
            "data/dada2_AbFilt_fasta.fasta",
            "Rmarkdown/Natrix_EvaluationAbundanceFilter.html",
            "Rmarkdown/Natrix_StatisticalAnalysis.html",
            "Rmarkdown/Natrix_IndicatorSpeciesDESeq2.html"
            ], 
                dn = dn,
                infile="PE_demux_joined")
'''


include:"rules/dada2PairedQiimeSimple.smk"
include:"rules/dada2NoFiltPairedQiimeSimple.smk"
include:"rules/deblurQiime.smk"
include:"rules/deblurSilvaQiime.smk"
include:"rules/extractRepSeqs.smk"
include:"rules/fastqc_raw.smk"
include:"rules/importFastq.smk"
include:"rules/JoinPairsQiime.smk"
include:"rules/multiqc_raw.smk"
include:"rules/summarizeFT_Qiime.smk"
include:"rules/SummarizeTrimmingQiime.smk"
include:"rules/trimAdapterQiime.smk"
include:"rules/trimQualityQiime.smk"
include:"rules/visualizeDenoiseStatsQiime.smk"
include:"rules/mapPhylogenyQiime.smk"
include:"rules/diversityCoreMetricsQiime.smk"
include:"rules/taxonomyBLASTQiime.smk"
include:"rules/summarizeRS_Qiime.smk"
include:"rules/taxBarplotQiime.smk"
include:"rules/thresholdSurvivorsDeblur.smk"
include:"rules/thresholdSurvivorsDada2.smk"
include:"rules/thresholdSurvivorsDada2NoFilt.smk"
include:"rules/thresholdSurvivorsDeblurSilva.smk"
include:"rules/combineSurvivors.smk"
include:"rules/createBiomAndRDS.smk"
include:"rules/testDada2Threshold_1.smk"
include:"rules/testDada2Threshold_2.smk"
include:"rules/testDada2NoFiltThreshold_1.smk"
include:"rules/testDada2NoFiltThreshold_2.smk"
include:"rules/testDeblurThreshold_1.smk"
include:"rules/testDeblurThreshold_2.smk"
include:"rules/testDeblurSilvaThreshold_1.smk"
include:"rules/thresholdBetaDiv.smk"
include:"rules/testDeblurSilvaThreshold_2.smk"
include:"rules/thresholdAlphaDiv.smk"
include:"rules/Natrix_EvaluationAbundanceFilter.smk"
include:"rules/AbundanceFilter.smk"
include:"rules/mapPhylogenyQiime_AbFilt.smk"
include:"rules/Natrix_StatisticalAnalysis.smk"
include:"rules/Natrix_IndicatorSpeciesDESeq2.smk"
