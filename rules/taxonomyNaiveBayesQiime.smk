# Porthmeus
# 07.06.20

# add taxonomic information to the otus

rule taxonomyQiime:
    input:
        classifier ="misc/silva-132-99-515-806-nb-classifier.qza", 
        #"misc/silva-132-99-nb-classifier.qza",
        #"misc/gg-13-8-99-515-806-nb-classifier.qza",
        RS = "data/PE_denoise_{dn}_RS.qza"
    output:
        tax = "data/PE_denoise_{dn}_Tax.qza",
        taxVis = "reports/PE_denoise_{dn}_Tax.qzv"
    log: "logs/taxonomyQiime_{dn}.log"
    conda: "../envs/qiime2.yaml"
    threads: 4
    shell:
        """
        qiime feature-classifier classify-sklearn \
                --i-classifier {input.classifier} \
                --i-reads {input.RS} \
                --p-n-jobs {threads} \
                --p-reads-per-batch 1000 \
                --o-classification {output.tax} &> {log}
        qiime metadata tabulate \
                --m-input-file {output.tax} \
                --o-visualization {output.taxVis} &> {log}
        """

