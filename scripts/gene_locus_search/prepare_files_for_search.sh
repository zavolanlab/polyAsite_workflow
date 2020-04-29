#!/bin/bash
# 24-JUN-2020, Alexander Kanitz, alexander.kanitz@alumni.ethz.ch

# This script will process primary clusters from the polyAsite annotation
# pipeline, extract the header, sort the clusters by genomic position and
# generate an index with `tabix`.

# Version requirements
# tabix 0.2.5 (r1005)
# R 3.6.0 (2019-04-26)
# R packages:
# - BiocGenerics_0.30.0
# - rtracklayer_1.44.0

# Note: `Rscript`, `bgzip` and `tabix` need to be available in $PATH

# Parameters
cluster_file="/scicore/home/zavolan/herrmchr/PolyASite/Data/atlas/HomoSapiens/GRCh38-96/20190808_3READS_5p_corr/clusters.support.tsv.gz"
anno_file="/scicore/home/zavolan/GROUP/PolyA/PolyASite_regularUpdates/annotation/Homo_sapiens.GRCh38.96.chr_patch_hapl_scaff.gtf"
design_file="/scicore/home/zavolan/herrmchr/PolyASite/PolyASite_Dev/samples_human_1perSeries.tsv"
output_prefix="clusters"
anno_filtering_script="/scicore/home/zavolan/herrmchr/PolyASite/PolyASite_regularUpdates/scripts/gene_locus_search/gtf_extract_genes.sh"
processing_script="/scicore/home/zavolan/herrmchr/PolyASite/PolyASite_regularUpdates/scripts/gene_locus_search/process_primary_clusters.R"

# Process annotations
echo "Processing annotations..."
anno_filtered="${output_prefix}.annotations_filtered.gtf"
"$anno_filtering_script" "$anno_file" "$anno_filtered"

# Process cluster file
echo "Processing clusters..."
processed_clusters="${output_prefix}.processed.tsv"
Rscript "$processing_script" \
	"$cluster_file" \
	"$anno_filtered" \
	"$design_file" \
	"$processed_clusters"

# Extract header
echo "Extracting header..."
header="${output_prefix}.header.tsv"
head -n 1 "$processed_clusters" > "$header"

# Sort clusters
echo "Sorting and compressing clusters..."
sorted_clusters="${output_prefix}.sorted.tsv.gz"
tail -n "+2" "$processed_clusters" | sort -k "1,1" -k "2,2n" | bgzip > "$sorted_clusters"

# Index clusters
echo "Indexing clusters..."
tabix -p "bed" -f "$sorted_clusters"
