# Porthmeus
# 05.06.20

# export the feature table and representative sequence table

rule exportBasicDataQiime:
    input:
        FT = "data/{data}_FT.qza",
        RS = "data/{data}_RS.qza"
    output:
        FT = "data/{data}_FT.biom",
        RS = "data/{data}_RS.fas"
    log: "logs/exportBasicDataQiime_{data}.smk"
    params:
        tempDir = "tempDirExportBasicQiime"
    conda:
        "../envs/qiime2.yaml"
    shell:
        """
        qiime tools export \
                --input-path {input.FT} \
                --output-path {params.tempDir} &> {log}
        mv {params.tempDir}/feature-table.biom {output.FT} &> {log}
        biom
        qiime tools export \
                --input-path {input.RS} \
                --output-path {output.RS} &> {log}
        mv {params.tempDir}/dna-sequences.fasta {output.RS} &> {log}
        rm -r {params.tempDir}
        """

