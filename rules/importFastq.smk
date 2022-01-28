# Porthmeus
# 22.05.20

# import the data from the illumina sequencing into a useful format for qiime2

rule importFastq:
    input: "FileManifest.tab"
    output: temp("data/PE_demux.qza")
    log: "logs/importFastq.log"
    conda: "../envs/qiime2.yaml"
    shell:
        "qiime tools import \
                --type 'SampleData[PairedEndSequencesWithQuality]' \
                --input-path {input} \
                --output-path {output} \
                --input-format PairedEndFastqManifestPhred33V2 &> {log}"
