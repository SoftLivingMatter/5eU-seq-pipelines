# 5eU-seq-pipelines

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
snakemake --profiler cluster
# OR
snakemake --profiler local
```
Installation should take a few minutes.

## Overview
Four types of analyses can be accomplished with the workflow, end_cov, cleavage,
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

## Additional Information
Project is under development, issues and pull requests are welcome to improve the
code.  Only tested on linux systems, container support requires singularity.
Runtimes vary based on sequencing depth resources available.
