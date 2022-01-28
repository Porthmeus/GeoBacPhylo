The plot shows the effect size of the morisita difference (beta diversity) between the Hydra microbiom samples within and outside the same geosite after the filtering and denoising by the different algorithms, using different filters.

Filters for deblur: triming length after merging with values ranging from 150-300 bp in 15 bp steps.
Deblur further filters by checking against different 16S-RNA databases and whether it is likely that the sequence obtained in the experiment is originating from an 16S sequence. The default (deblur) is to check against the greengenes database, but since the last update on the data base is long ago, I used a second database, silva_13.1 (deblurSilva) which is newer and actively maintained.

For dada2 there is the possibility to check for base quality (phred) in the read before merging. The filter for this option ranged between 0-20 increasing by 2.
Furthermore, there is the possibility to check for divergence in the forward and reverse read while merging. Default value for this is maximum 2 bp difference (dada2) and I also compared by not filtering these reads (dadaNoFilt).
