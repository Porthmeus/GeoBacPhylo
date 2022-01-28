# Porthmeus
# 26.05.20

# join the fwd and rev reads

rule JoinPairsQiime:
    input:"data/PE_demux_cut.qza"
    output: 
        data = temp("data/PE_demux_cut_joined.qza"),
        report = "reports/PE_demux_cut_joined.qzv"
    log: "logs/qiimeJoinPairs.log"
    conda: "../envs/qiime2.yaml"
    shell:
        """
        qiime vsearch join-pairs \
                --i-demultiplexed-seqs {input}\
                --o-joined-sequences {output.data} &> {log}
        qiime demux summarize \
                --i-data {output.data} \
                --o-visualization {output.report}
        """
