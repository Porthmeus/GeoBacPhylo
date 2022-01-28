# Porthmeus
# 30.07.20

# recalculate the phylogenetic tree for the abundance filtered data

rule mapPhylogenyQiime_AbFilt:
    input:
        "data/dada2_AbFilt_fasta.fasta"
    output:
        qza = "data/dada2_AbFilt_fasta.qza",
        aln = "data/dada2_AbFilt_fasta_aln.qza",
        msk_aln = "data/dada2_AbFilt_fasta_mskAln.qza",
        ur_tree = "data/dada2_AbFilt_fasta_urTree.qza",
        tree = "data/dada2_AbFilt_fasta_tree.qza"
    log: "logs/dada2_AbFilt_fasta.log"
    conda: "../envs/qiime2.yaml"
    shell:
        """qiime tools import \
                --type FeatureData[Sequence] \
                --input-path {input} \
                --output-path {output.qza} &> {log}
        qiime phylogeny align-to-tree-mafft-fasttree \
                --i-sequences {output.qza} \
                --o-alignment {output.aln} \
                --o-masked-alignment {output.msk_aln} \
                --o-tree {output.ur_tree} \
                --o-rooted-tree {output.tree} &> {log}"""
