# Porthmeus
# 07.06.20

# summarize feature tables and the representative sequences of the data

rule summarizeRS_Qiime:
    input:
        RS = "data/PE_denoise_{dn}_RS.qza",
        meta = "MetaData.csv"
    output:
        RS = "reports/PE_denoise_{dn}_sumRS.qzv"
    log: "logs/summarizeRS_Qiime_{dn}.log"
    conda: "../envs/qiime2.yaml"
    shell:
        """
        qiime feature-table tabulate-seqs \
                --i-data {input.RS} \
                --o-visualization {output.RS} &> {log}
        """
