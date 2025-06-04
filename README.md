# 5eU-seq-pipelines
[![DOI](https://zenodo.org/badge/855379253.svg)](https://doi.org/10.5281/zenodo.14908465)

A collection of analyses for EU enriched sequence experiments, including
cleavage measurements and modifications.  Execution is specified as a snakemake
workflow and controlled with the `config/config.yaml`, including reference
locations and which analyses to run.  A combination of all references and
junctions are produced.

## Usage

Install [mamba](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install)
which is recommended for snakemake.

Create a snakemake environment, using version 7
```bash
mamba create -c conda-forge -c bioconda -n snake snakemake=7.32.4
```

Next update the `config.yaml` to match your experimental setup.  The workflow
can be run locally or through a cluster by choosing the appropriate profile,
```bash
snakemake --profile cluster
# OR
snakemake --profile local
```
Installation should take a few minutes.

## Overview
Three types of analyses can be accomplished with the workflow, end_cov, cleavage,
methyl.  They are included together due to the overlap between the analyses,
it is unlikely all will be run on the same input files.

In all cases, fastqs are trimmed with trimmomatic in the rule `trim`.  Trimmed
fastqs are then aligned to the specified refrerence using STAR and filtered for
quality.

For end_cov analysis, samtools is used to filter bams which contain read 1, 2
or both.  Read 2 is further filtered with bedtools for 3 prime only.

For the cleavage analysis, featureCounts is used to measure the number of sequences
overlapping and non-overlapping with the specified saf junctions.  Those values
are tabulated to determine the fraction of reads spanning a junction.  Additionally, the
sample name and time can be parsed from the file name, using the config setting
`cleavage_sample_regex`.

Analysis of methyl scores use the `scripts/rna_mod_score.py` function
which analyzes output of end coverage files to estimate a scores representing
the likelihood of a modification at that position.

## Inference of modification sites
### 2OMe
Visualize the 5prime.read1.only bedgraph (in end_cov output folder) in IGV and
look for dips in the read counts.
Be careful for +1 offset in IGV (look at the counts and nucleotide numbers in
bedgraph to confirm that the nucleotide you identified with the dip in IGV is
correct)

### PseudoU: visualize the bam files in IGV (must also download the bam.bai
file), you can find this in the align_bam output folder and look for sites with
high proportion of deletions. Alternatively, you can look at the genes.tsv
output file and set a cutoff for sites that have, for example, more than 20%
deletion fraction. Sites identified must be uridines!
Cutoff can be more lenient or strict depending on amount of noise tolerable For
some organisms, some sites at least are known (e.g. worm) so should positive
control check those.

Important note: when trying to find sites annotated in the
literature in IGV, the nucleotide numbering may not make sense! this is because
there are offsets as the genomes get updated. So (for example) site 2727 in the
literature may still exist in your genome, but offset by 10 nucleotides. To
determine what this offset is (which is not always the same across the same
genome e.g. in 28S): use numbering information between literature annotated
sites to find them in IGV genome (e.g. if I know this site 1 is 10 nt away
from site 2, and site 3 is 4 nt away from site 2, etc..)

## Additional Information
Project is under development, issues and pull requests are welcome to improve the
code.  Only tested on linux systems, container support requires singularity.
Runtimes vary based on sequencing depth resources available.
