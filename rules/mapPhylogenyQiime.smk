# Porthmeus
# 26.05.20

# add the phylogeny to the calculated representative sequences

rule mapPhylogenyQiime:
    input: "data/PE_denoise_{dn}_RS.qza"
    output:
        aln = "data/PE_denoise_{dn}_aln.qza",
        msk_aln = "data/PE_denoise_{dn}_mskAln.qza",
        ur_tree = "data/PE_denoise_{dn}_urTree.qza",
        tree = "data/PE_denoise_{dn}_tree.qza"
    log: "logs/mapPhylogenyQiime_{dn}.log"
    conda: "../envs/qiime2.yaml"
    shell:
        "qiime phylogeny align-to-tree-mafft-fasttree \
                --i-sequences {input} \
                --o-alignment {output.aln} \
                --o-masked-alignment {output.msk_aln} \
                --o-tree {output.ur_tree} \
                --o-rooted-tree {output.tree} &> {log}"

