#!/usr/bin/env Rscript
# 24-JUN-2019, Alexander Kanitz, alexander.kanitz@alumni.ethz.ch

# Usage: Rscript process_primary_clusters.R <PRIMARY_CLUSTERS_FILE> <ANNOTATION_FILE> <DESIGN_FILE> <OUTPUT_FILE>

# Notes: No checking of input parameters is performed. Check the version info below and ensure that you have a compatible R version available.

# Required inputs (set below):
# 1. Primary cluster file (output of polyAsite analysis pipeline); requires at least 8 columns, including 'start', 'end', 'strand', 'rep' and a column that includes information on whether a site might be an internal priming site (currently 'repSite_signals'; value is hardcoded below, as well as the actual value that indicates that a site might be such a site, currently 'IPCandidate'); also requires chromosome information in the first and some sort of score (currently 'total_tpm') in the sixth column
# 2. Matching gene annotations file (GTF), filtered for gene entries only
# 3. Sample info / design table file; requires columns 'sample_id', 'protocol', 'source', 'title' and 'treatment'
# 4. Output filename (output is *not* compressed!)

# Processing steps carried out (in indicated order):
# 1. Potential internal priming sites are removed
# 2. Columns are re-arranged to yield a BED-like format (inclusion of a column indicating unique names for sites/clusters, from reference sequence, position of "representative" site and strand, e.g. 3:12345:+)
# 3. Two columns indicating names/symbols and Ensembl gene identifiers of any genes that the cluster/site is contained in; NA if it does not overlap any gene, multiple overlapped gene names/ids are separated by the pipe character, e.g. "GENE1|GENE2"); if multiple genes are overlapped, there *should* be the same number of gene names and Ensembl gene identifiers, and corresponding names *should* appear in the same order (not extensively tested), e.g. for "SiteXY" "GENE1|GENE2" "ID1|ID2", we should assume that "GENE1" and "ID1" are corresponding
# 4. Generic sample names (i.e., currently any column names after 7th) are replaced with actual sample names, pasted to other useful sample information via the pipe character, specifically: protocol, source, title and treatment; e.g,: SAMPLE_1|PROTOCOL_X|SOME_TISSUE|DESCRIPTIVE_TITLE|SOME_TREATMENT_CONDITION

# Output:
# Uncompressed file with the following columns (see info above for more info):
# chrom  chromStart  chromEnd  name  score  strand  rep  repSite_signals gene_name gene_id <sample_1_name> <sample_2_name> <...> <sample_n_name>

# User versions:
# R 3.6.0 (2019-04-26)
# BiocGenerics_0.30.0
# rtracklayer_1.44.0

# Load required packages
suppressMessages(library('rtracklayer'))

# Process CLI arguments
args = commandArgs(trailingOnly=TRUE)
path.cluster_file = args[1]
path.annotation_file = args[2]
path.design_file = args[3]
output_file = args[4]

# Hardcoded parameters
index_of_first_sample_column <- 14 # AFTER CONVERSION TO NEW FORMAT :)
index_of_column_containing_ip_info <- 10 # BEFORE CONVERSION TO NEW FORMAT
search_string_ip_info <- "IPCandidate"

# Read data
clusters <- read.delim(path.cluster_file, header=TRUE, stringsAsFactors=FALSE)
sample_info <- read.delim(path.design_file, header=TRUE, stringsAsFactors=FALSE)
anno <- import(path.annotation_file, format="gtf")

# Filter out internal priming sites
filter <- grep(search_string_ip_info, clusters[[index_of_column_containing_ip_info]])
if (! length(filter) ) {
    cat("No IP candidates found!\n")
} else {
	clusters <- clusters[-filter, ]
}

# Convert clusters table to BED-like format
nms <- paste(clusters[[1]], clusters[['rep']], clusters[['strand']], sep=":")
clusters <- cbind(
    chrom = clusters[[1]],
    chromStart = clusters[['start']],
    chromEnd = clusters[['end']],
    name = nms,
    score = clusters[['total_tpm']],
    strand = clusters[['strand']],
    rep = clusters[['rep']],
    frac_samples = clusters[['frac_samples']],
    nr_prots = clusters[['nr_prots']],
    annotation = clusters[['annotation']],
    clusters[index_of_column_containing_ip_info:ncol(clusters)]
)

# Add gene info
clusters.bed <- GRanges(seqnames=clusters[['chrom']], ranges=IRanges(clusters[['chromStart']], clusters[['chromEnd']]), strand=clusters[['strand']])
hits <- suppressWarnings(findOverlaps(clusters.bed, anno))
gene_names <- anno$gene_name[to(hits)]
gene_ids <- anno$gene_id[to(hits)]
dat.gene <- data.frame(cluster=from(hits), gene_name=gene_names, gene_id=gene_ids)
dat.aggr <- aggregate(dat.gene[2:3], by=list(dat.gene[[1]]), FUN=paste, collapse="|")
gene_info <- data.frame(gene_name=rep(NA, nrow(clusters)), gene_id=rep(NA, nrow(clusters)))
gene_info[['gene_name']][dat.aggr[[1]]] <- dat.aggr[['gene_name']]
gene_info[['gene_id']][dat.aggr[[1]]] <- dat.aggr[['gene_id']]
clusters <- cbind(clusters[1:10], gene_info, clusters[11:ncol(clusters)])

# Replace sample names
sample_names <- paste(sample_info[['sample_id']], sample_info[['protocol']], sample_info[['source']], sample_info[['title']], sample_info[['treatment']], sep="|")
colnames(clusters)[index_of_first_sample_column:ncol(clusters)] <- sample_names

# Write output
write.table(clusters, output_file, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
