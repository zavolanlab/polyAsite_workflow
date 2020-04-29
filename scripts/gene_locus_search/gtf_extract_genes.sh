#!/bin/bash
# 24-JUN-2020, Alexander Kanitz, alexander.kanitz@alumni.ethz.ch

# Usage: gtf_extraxt_genes.sh <INPUT_FILE> <OUTPUT_FILE>

# This script will extract `gene` entries from an uncompressed GTF file.

# Parameters
input_file="$1"
output_file="$2"

# Extract gene entries
awk '$3 == "gene"' "$input_file" > "$output_file"
