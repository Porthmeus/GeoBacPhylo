# Porthmeus
# 02.06.20

# calculate the different diversity scores

rule diversityCoreMetrics:
    input: 
        features = "data/PE_denoise_{dn}_FT.qza",
        phylo = "data/PE_denoise_{dn}_tree.qza",
        meta = "MetaData.csv"
    output:
        rareTable = "data/distance/{dn}_rarefied_table.qza",
        faithVec = "data/distance/{dn}_faith_pd_vector.qza",
        obsVec = "data/distance/{dn}_observed_otus_vector.qza",
        shnVec = "results/distance/{dn}_shannon_vector.qza",
        evnVec = "results/distance/{dn}_evenness_vector.qza",
        unifracTab_UW ="results/distance/{dn}_unweighted_unifrac_distance_matrix.qza",
        unifracTab_W = "results/distance/{dn}_weighted_unifrac_distance_matrix.qza",
        jaccardTab = "results/distance/{dn}_jaccard_distance_matrix.qza",
        brayCurtisTab = "results/distance/{dn}_bray_curtis_distance_matrix.qza",
        unifracPCA_UW = "results/distance/{dn}_unweighted_unifrac_pcoa_results.qza",
        unifracPCA_W = "results/distance/{dn}_weighted_unifrac_pcoa_results.qza",
        jaccardPCA = "results/distance/{dn}_jaccard_pcoa_results.qza",
        brayCurtisPCA = "results/distance/{dn}_bray_curtis_pcoa_results.qza",
        unifracEmp_UW = report("results/distance/{dn}_unweighted_unifrac_emperor.qzv",
                caption = "../captions/diverstiyCoreMetrics.rst",
                category = "Diversity Core metrics"),
        unifracEmp_W = report("results/distance/{dn}_weighted_unifrac_emperor.qzv",
                caption = "../captions/diverstiyCoreMetrics.rst",
                category = "Diversity Core metrics"),
        jaccardEmp = report("results/distance/{dn}_jaccard_emperor.qzv",
                caption = "../captions/diverstiyCoreMetrics.rst",
                category = "Diversity Core metrics"),
        brayCurtisEmp =report( "results/distance/{dn}_bray_curtis_emperor.qzv",
                caption = "../captions/diverstiyCoreMetrics.rst",
                category = "Diversity Core metrics")
    log: "logs/diversityCoreMetric_{dn}.log"
    params:
        depth = 1204, # 10^(mean(log10(FC)_dada2+1)-2*var(log10(FC)_dada2+1))+1
        outdir = "result/distance_{dn}/"
    conda:
        "../envs/qiime2.yaml"
    shell:
        "qiime diversity core-metrics-phylogenetic \
                --i-phylogeny {input.phylo} \
                --i-table {input.features} \
                --p-sampling-depth {params.depth} \
                --m-metadata-file {input.meta} \
                --o-rarefied-table {output.rareTable} \
                --o-faith-pd-vector {output.faithVec} \
                --o-observed-otus-vector {output.obsVec} \
                --o-shannon-vector {output.shnVec} \
                --o-evenness-vector {output.evnVec} \
                --o-unweighted-unifrac-distance-matrix {output.unifracTab_UW} \
                --o-weighted-unifrac-distance-matrix {output.unifracTab_W} \
                --o-jaccard-distance-matrix {output.jaccardTab} \
                --o-bray-curtis-distance-matrix {output.brayCurtisTab} \
                --o-unweighted-unifrac-pcoa-results {output.unifracPCA_UW} \
                --o-weighted-unifrac-pcoa-results {output.unifracPCA_W} \
                --o-jaccard-pcoa-results {output.jaccardPCA} \
                --o-bray-curtis-pcoa-results {output.brayCurtisPCA} \
                --o-unweighted-unifrac-emperor {output.unifracEmp_UW} \
                --o-weighted-unifrac-emperor {output.unifracEmp_W} \
                --o-jaccard-emperor {output.jaccardEmp} \
                --o-bray-curtis-emperor {output.brayCurtisEmp} &> {log}"

