# Porthmeus
# 07.06.20

# add taxonomic information to the otus

rule taxonomyQiime:
    input:
        refSeq ="misc/silva132_refSeqs.qza",
        taxTable = "misc/silva132_taxonomy.qza",
        RS = "data/PE_denoise_{dn}_RS.qza"
    output:
        tax = "data/PE_denoise_{dn}_Tax.qza",
        taxVis = "reports/PE_denoise_{dn}_Tax.qzv"
    log: "logs/taxonomyQiime_{dn}.log"
    conda: "../envs/qiime2.yaml"
    shell:
        """
        qiime feature-classifier classify-consensus-blast \
                --i-query {input.RS} \
                --i-reference-reads {input.refSeq} \
                --i-reference-taxonomy {input.taxTable} \
                --o-classification {output.tax} &> {log}
        qiime metadata tabulate \
                --m-input-file {output.tax} \
                --o-visualization {output.taxVis} &> {log}
        """

