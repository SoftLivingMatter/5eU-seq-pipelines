#!/bin/bash

BASE_DIR=/projects/BRANGWYNNE/projects/genomics/primary_data/2023/01_January/20230117-MCF10-EUseq/data/03_end_cov

echo 'testing ribometh'
python rna_mod_score.py \
	--output all_2methylPositions_new.csv \
	--column-name methylPosition \
	--exclude-score methyl_positions.txt \
	--query-sites methyl_positions.txt \
	$BASE_DIR/*5prime.read1only.bedgraph
# --query-sites <(echo 4454) \
python cmp_csv.py all_2methylPositions.csv all_2methylPositions_new.csv methylPosition
echo 'done ribometh'
