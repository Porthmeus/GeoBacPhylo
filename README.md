# Analysis Pipeline for the Hydra samples in Hungary

This is the corresponding github project to a paper "Population Differences and Host Species Predict Variation in the Diversity of Host-Associated" (doi: 10.3389/fmicb.2022.799333). To run the pipeline:

- clone the github repository `git clone git@github.com:Porthmeus/GeoBacPhylo.git`
- install [snakmake](https://snakemake.readthedocs.io/en/stable/)
- download the raw files from the SRA ([PRJNA795254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA795254)) and place them in the `raw` directory of this project
- invoke snakemake with `snakemake -s Analysis.smk`

This will run the basic analysis in order to create the readcount matrix, 16S sequencing fasta and the taxonomic annotation file. The latter two are available as supplementary tables in the manuscript.


In order to recreate the plots and statistics of the manuscript one can run the scripts in the `Manuscript` directory after the relevant data has been generated.

## Data

You can find the most important data from the experiment, without rerunning the analysis, in `Manuscript/Data`. 

| File                         | Description                                                        |
| ---------------------------- | ------------------------------------------------------------------ |
| `dada2_AbFilt_fasta.fasta` | Contains the sequence information for the reconstructed ASVs       |
| `dada2_AbFilt_FT.csv`      | The feature count table, containing the abundances of the bacteria |
| `dada2_AbFilt_meta.csv`    | Metadata for the samples                                           |
| `dada2_AbFilt_Tax.csv`     | Taxonomic information for the ASVs                                 |

For further questions and help raise an issue here on github or contact Jan Taubenheim (corresponding author of the original manuscript).


