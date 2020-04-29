# Segemehl mapping
## Include (Current setup)
The Snakefile can be called through the include directive in the main Snakefile of the polyASite_regular_updates.     
Samples-, script-, log-, and annotation directories have to match the corresponding directories in the main Snakefile!

## Standalone
Change directories, especially in function 'get_valid_reads' an input directory from where to get the samples from has to be specified!! Design file, config and cluster_config have to be created and specified accordingly.
scripts directory has to be provided with the following scripts:
- get_lines_w_pattern.sh
- gtf_exons_bed.1.1.2.R
- sam_remove_duplicates_inferior_alignments_multimappers.1_5.pl
- sam_trx_to_sam_gen.pl
- validation_fasta.py
- rs-bam2bed.py
