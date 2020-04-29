#configfile: "config.yaml"
# don't provide a default config file
# to ensure that the executing person knows what it is doing
# and provides the intended information

from snakemake.utils import makedirs
from snakemake.utils import listfiles

import pandas as pd
import numpy as np
import string
import random
import os

################################################################################
# Mapping subpipeline
################################################################################
include: "segemehl/Snakefile"

################################################################################
# Functions for exceptions/branching
################################################################################

#-------------------------------------------------------------------------------
# A-seq2 samples go through 5p-Adapter trimming (select_for_valid_5p_configuration)
#-------------------------------------------------------------------------------
def trim_5p_adapter_input(wildcards):
    if samples.loc[wildcards.sample, "protocol"] == "QuantSeq_REV":
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".fa.gz")
    else:
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".fa.gz")

def trim_adapter_input(wildcards):
    if samples.loc[wildcards.sample, "protocol"] == "A-seq2":
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".5ptrimmed.A-seq2.fa.gz")
    elif samples.loc[wildcards.sample, "protocol"] == "3'READS":
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".5ptrimmed.3READS.fa.gz")
    elif samples.loc[wildcards.sample, "protocol"] == "PAS-Seq":
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".leadingTs_trimmed.fa.gz")
    else:
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".fa.gz")

def get_reverse_compl_input(wildcards):
    if( samples.loc[wildcards.sample, "protocol"] == "PAS-Seq" or
        samples.loc[wildcards.sample, "protocol"] == "SAPAS"):
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".trimmed_add_nuc.fa.gz")
    elif( samples.loc[wildcards.sample, "protocol"] == "DRS" or
          samples.loc[wildcards.sample, "protocol"] == "3P-Seq"):
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".fa.gz")
    elif( samples.loc[wildcards.sample, "protocol"] == "QuantSeq_REV"):
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".5ptrimmed.fa.gz")
    else:
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".trimmed.fa.gz")

def get_valid_3p_3PSeq_file(wildcards):
    if samples.loc[wildcards.sample, "reverse_compl"]:
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".trimmed.rev_cmpl.fa.gz")
    else:
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".fa.gz")

def get_input_polyAtail_removal(wildcards):
    if samples.loc[wildcards.sample, "protocol"] == "3P-Seq":
        return os.path.join(config["samples_dir"],
                               wildcards.sample,
                               wildcards.sample + ".valid_3p_configuration.fa.gz")
    elif samples.loc[wildcards.sample, "protocol"] == "PAPERCLIP":
        if pd.isna( samples.loc[wildcards.sample, "fiveAdapter"]):
            return os.path.join(config["samples_dir"],
                               wildcards.sample,
                               wildcards.sample + ".fa.gz")
        else:
            return os.path.join(config["samples_dir"],
                               wildcards.sample,
                               wildcards.sample + ".5ptrimmed.fa.gz")
    elif samples.loc[wildcards.sample, "protocol"] == "3'-Seq (Mayr)":
        return os.path.join(config["samples_dir"],
                               wildcards.sample,
                               wildcards.sample + ".trimmed.fa.gz")
    else:
        return os.path.join(config["samples_dir"],
                                     wildcards.sample,
                                     wildcards.sample + ".trimmed.rev_cmpl.fa.gz")

#--------------------------------------------------------------------------------
# Reverse complement will NOT be applied for Aseq, Mayr and 3P-seq;
# those go from trim_3p_adapter directly to get_valid_reads
# 3READS,3P-Seq,QuantSeq_REV,PAPERCLIP, Mayr: additionally, the poly(A) tail is trimmed
#--------------------------------------------------------------------------------
def get_reads_after_trimming(wildcards):
    if(samples.loc[wildcards.sample, "protocol"] == "A-seq" or
       samples.loc[wildcards.sample, "protocol"] == "3P-seq"):
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".trimmed.fa.gz")
    elif samples.loc[wildcards.sample, "protocol"] == "3'-Seq (Mayr)":
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".trimmed_tail.fa.gz")
    elif( samples.loc[wildcards.sample, "protocol"] == "3'READS" or
          samples.loc[wildcards.sample, "protocol"] == "3P-Seq" or
          samples.loc[wildcards.sample, "protocol"] == "PAPERCLIP" or
          samples.loc[wildcards.sample, "protocol"] == "QuantSeq_REV"):
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".trimmed_tail.fa.gz")
    else:
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".trimmed.rev_cmpl.fa.gz")

#-------------------------------------------------------------------------------
# the valid reads change depending on whether reads are
# also length filtered
#-------------------------------------------------------------------------------
def get_valid_reads(wildcards):
    if(samples.loc[wildcards.sample, "protocol"] == "A-seq" or
       samples.loc[wildcards.sample, "protocol"] == "3'-Seq (Mayr)"):
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".valid_maxLength.fa.gz")
    else:
        return os.path.join(config["samples_dir"],
                            wildcards.sample,
                            wildcards.sample + ".valid.fa.gz")




def get_ds_patterns_for_ipAssignment(wildcards):
    if samples.loc[wildcards.sample, "protocol"] == "3'READS":
        return " ".join(["--ds_pattern=%s"
                         % pat for pat in config['IP.downstream_patterns.3READS'] ]) \
                         if config['IP.downstream_patterns.3READS'] is not None else ""
    else:
        return " ".join(["--ds_pattern=%s"
                         % pat for pat in config['IP.downstream_patterns'] ]) \
                         if config['IP.downstream_patterns'] is not None else ""

def get_excluded_chromosomes(wildcards):
    if config['excluded_chr'] is not None:
        excluded = {i:1 for i in config['excluded_chr']}
    else:
        excluded = {}
    if samples.loc[wildcards.sample, "sex"] == "F":
        excluded[ config['female_chr'] ] = 1
    else:
        if config['female_chr'] in excluded:
            del excluded[ config['female_chr'] ]
    excluded_list = list(excluded.keys())
    if len(excluded_list) == 0:
        return ""
    else:
        return("--exclude=" + ":".join(excluded_list))


#-------------------------------------------------------------------------------
# local rules
#-------------------------------------------------------------------------------
localrules: create_log_dir, create_log_dir_atlas, download_fastq_se, download_fastq_pe,
            change_fastq_to_fq_se, download_genome, download_annotation, fetch_chr_sizes_ucsc,
            make_track_info, complete_preprocessing, complete_clustering, complete_tracks, finish
samples = pd.read_table(config['atlas.samples_table'], index_col=0, comment='#')

atlas_outputdir = os.path.join(config['atlas_dir'],
                   config['organism'],
                   config['genome'],
                   config['atlas.release_name'])
################################################################################
# target rule
################################################################################
rule finish:
    ##LOCAL##
    ##No Singularity support required##
    input:
        #ip_cnt = expand( os.path.join(config["samples_dir"],
        #                              "counts",
        #                              "{sample}_" + config['genome'] + ".ip3pSites.out"),
        #                 sample = samples.index),
        #track_info = expand( os.path.join(config["samples_dir"],
        #                           "{sample}",
        #                           config['genome'],
        #                           "{sample}.track_info.txt"),
        #                sample = samples.index),
        prepro_cmplt = expand( os.path.join(config["samples_dir"],
                                   "{sample}",
                                   config['genome'],
                                   "{sample}.prepro_cmplt.txt"),
                        sample = samples.index),
        clst_cmplt = os.path.join(atlas_outputdir,
                                    "clst_cmplt.txt"),
        tracks_cmplt = os.path.join(atlas_outputdir,
                                    "tracks_cmplt.txt")
        #noBG_cnt = expand( os.path.join(config["samples_dir"],
        #                            "counts",
        #                            "{sample}_" + config['genome'] + ".noBG3pSites.out"),
        #                 sample = samples.index),
        #cluster_stats = os.path.join( atlas_outputdir,
        #                            "counts",
        #                            "clusters.stats.out" ),
        #track_info_cl = os.path.join(atlas_outputdir,
        #                            ("clusters."+ config['genome'] + "_"
        #                              + config['atlas.release_name']
        #                              + ".track_info.txt")),
        #final_atlas = os.path.join(atlas_outputdir,
        #                            "clusters.bed.gz")

################################################################################
# individual rules (if possible in chronological order)
################################################################################

#-------------------------------------------------------------------------------
# create dir for logfiles (samples)
#-------------------------------------------------------------------------------
rule create_log_dir:
    ##LOCAL##
    ##No Singularity support required##
    ''' This step creates the log directory, if necessary.
    This is required when jobs are submitted and the
    job output should be written to these files.
    '''
    params:
        cluster_samples_log = os.path.join(config["samples_dir"],
                                           "logs",
                                           "cluster_logs"),
        cluster_countings_log = os.path.join(config["samples_dir"],
                                             "logs",
                                             "cluster_logs",
                                             "counting")
    output:
        dirs_samples_created = touch(os.path.join(config["samples_dir"],
                                                  "logs",
                                                  "created_log_dir.out"))
    shell:
        '''
        mkdir -p {params.cluster_samples_log}
        mkdir -p {params.cluster_countings_log}
        '''
#-------------------------------------------------------------------------------
# create dir for logfiles (atlas)
#-------------------------------------------------------------------------------
rule create_log_dir_atlas:
    ##LOCAL##
    ##No Singularity support required##
    ''' This step creates the log directory, if necessary.
    This is required when jobs are submitted and the
    job output should be written to these files.
    '''
    params:
        cluster_atlas_log = os.path.join(atlas_outputdir,
                                         "logs",
                                         "cluster_logs")
    output:
        dirs_atlas_created = touch(os.path.join(atlas_outputdir,
                                                "logs",
                                                "created_log_dir.out")),
    shell:
        '''
        mkdir -p {params.cluster_atlas_log}
        '''

#-------------------------------------------------------------------------------
# download the genome sequence
#-------------------------------------------------------------------------------
rule download_genome:
    ##LOCAL##
    ##Singularity provided: zavolab_minimal:1, not tested##
    output:
        genome = os.path.join(config['genome_fasta_raw']),
        temp_genome = temp( "genome." + ''.join(random.choice(string.ascii_uppercase) for _ in range(6)) + ".fa.gz"),
        clean = os.path.join(config['genome_fasta'])
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        url = config['genome_fasta_url']
    resources:
        load = 20 # With "--resources load=100", max 5 instances of this rule are run in parallel!
    shell:
        '''
        wget -O {output.temp_genome} \
        {params.url} \
        &> /dev/null &&
        gzip -cd {output.temp_genome} \
        > {output.genome} &&
        sed 's/\s.*//' {output.genome} \
        > {output.clean}
        '''

#-------------------------------------------------------------------------------
# download the gene annotation file
#-------------------------------------------------------------------------------
rule download_annotation:
    ##LOCAL##
    ##Singularity provided: zavolab_minimal:1, not tested##
    output:
        anno = config['gene_annotation'],
        temp_anno = temp( "gene_anno." + ''.join(random.choice(string.ascii_uppercase) for _ in range(6)) + ".gtf.gz")
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        url = config['gene_anno_url']
    resources:
        load = 20 # With "--resources load=100", max 5 instances of this rule are run in parallel!
    shell:
        '''
        wget -O {output.temp_anno} \
        {params.url} \
        &> /dev/null &&
        gzip -cd {output.temp_anno} \
        > {output.anno}
        '''

#-------------------------------------------------------------------------------
# get filtered version of annotation
#-------------------------------------------------------------------------------
rule get_filtered_annotation:
    ##Singularity needed: perl, gzip##
    ## Singularity provided: zavolab_minimal:1, not tested ##
    input:
       anno = os.path.join(config['gene_annotation']),
       script = os.path.join( config["script_dir"],
                                      "rs-filter-gtf-by-type-and-support.pl")
    output:
        filtered_anno = config['gene_annotation_filtered']
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        type_id = config['transcript_biotype_id'],
        types = lambda wildcards: " ".join(["--type=" + i for i in config['transcript_type']]),
        tr_supp_level_id = "--support_level_id=" + config['transcript_support_level_id'] \
                                                 if config['transcript_support_level_id'] is not None else "",
        tr_supp_level = "--support_level=" + config['transcript_support_level'] \
                                           if config['transcript_support_level'] is not None else "",
        cluster_log = os.path.join(config['annotation_dir'],
                                   "filter_anno.log")
    shell:
        '''
        perl {input.script} \
        --type_id={params.type_id} \
        {params.types} \
        {params.tr_supp_level_id} {params.tr_supp_level} \
        {input.anno} \
        > {output.filtered_anno}
        '''

# ################################################################################
# # preprae mongoDB collection
# ################################################################################
# rule check_samples_in_mongoDB:
#     # check for each sample if it is already in the mongoDB
#     # (genome and organism specific)
#     input:
#         dirs_created = os.path.join(atlas_outputdir, "logs", "created_log.tmp" )
#     output:
#         checked_db = os.path.join(atlas_outputdir,
#                                   "logs",
#                                   "mongoDB_samples_checked.log")
#     params:
#         samples = samples.index,
#         organism = config['organism'],
#         genome = config['genome']


################################################################################
# samples preprocessing
################################################################################

#-------------------------------------------------------------------------------
# download fastq files (paired-end data)
#-------------------------------------------------------------------------------
rule download_fastq_pe:
    ##LOCAL##
    input:
        script = os.path.join( config["script_dir"],
                               "rs-download-fastq-files-from-ena-via-ascp.py")
    output:
        sample_fq = expand(os.path.join(config["samples_dir"],
                                        "{{sample_name}}",
                                        "{{sample_id}}_{read}.fastq.gz"),
                           read = [1,2])
    singularity:
        "docker://cjh4zavolab/aspera:5"
    params:
        outdir = os.path.join(config["samples_dir"],
                              "{sample_name}"),
        srr_id = "{sample_id}"
    resources:
        load = 20 # With "--resources load=100", max 5 instances of this rule are run in parallel!
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "download_fq",
                     "{sample_name}.log")
    shell:
        '''
        python3 {input.script} \
        --srr_id {params.srr_id} \
        --outdir {params.outdir}
        --paired \
        2> {log}
        '''

#-------------------------------------------------------------------------------
# download fastq file (single_end data)
#-------------------------------------------------------------------------------
rule download_fastq_se:
    ##LOCAL##
    input:
        script = os.path.join( config["script_dir"],
                               "rs-download-fastq-files-from-ena-via-ascp.py")
    output:
        sample_fq = os.path.join(config["samples_dir"],
                                 "{sample_name}",
                                 "{sample_id}.fastq.gz")
    singularity:
        "docker://cjh4zavolab/aspera:5"
    params:
        outdir = os.path.join(config["samples_dir"],
                              "{sample_name}"),
        srr_id = "{sample_id}"
    resources:
        load = 20 # With "--resources load=100", max 5 instances of this rule are run in parallel!
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "download_fq",
                     "{sample_name}.log")
    shell:
        '''
        python3 {input.script} \
        --srr_id {params.srr_id} \
        --outdir {params.outdir} \
        2> {log}
        '''

#-------------------------------------------------------------------------------
# convert file names:
# - from SRR id to GSM/SRA id
# - from fastq.gz to fq.gz
#-------------------------------------------------------------------------------
rule change_fastq_to_fq_se:
    ##LOCAL##
    #ATTENTION: For some samples, multiple SRR run files (aka fastq files) belong to
    #           a single sample. If this is the case, they are concatenated here
    #           before the softlink is established
    ##No Singularity support required##
    input:
        sample_fq = lambda wildcards: expand(os.path.join(config["samples_dir"],
					                  wildcards.sample,
                                                          "{srr}.fastq.gz"),
                                             srr = samples.loc[wildcards.sample, "SRR"].split(","))
    output:
        sample_fq = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 "{sample}.fq.gz")
    params:
        file_dir = os.path.join(config["samples_dir"],
                                "{sample}"),
        sample_srr = lambda wildcards: samples.loc[wildcards.sample, "SRR"],
        first_srr = lambda wildcards: samples.loc[wildcards.sample, "SRR"].split(",")[0],
        sample_id = "{sample}"
    shell:
        '''
        cd {params.file_dir}
        IFS=',' read -ra SRR <<< "{params.sample_srr}"
        if [[ "${{#SRR[@]}}" > "1" ]];then
        first_file="${{SRR[0]}}.fastq.gz"
        for i in $(seq 1 $((${{#SRR[@]}}-1))); do curr_file="${{SRR[$i]}}.fastq.gz"; cat ${{curr_file}} >> ${{first_file}};done
        fi
        ln -fs {params.first_srr}.fastq.gz {params.sample_id}.fq.gz
        cd -
        '''

#-------------------------------------------------------------------------------
# convert fq to fasta
# hint: I (Ralf) do not use fastq_to_fasta anymore because I had issues with the
# number of output reads and I believe it has to do with the -Q33 option
# I use to indicate the correct offset for the quality scores
# (see here: https://www.biostars.org/p/120311/ )
#-------------------------------------------------------------------------------
rule fq2fasta_se:
    input:
        dirs_samples_created = os.path.join(config["samples_dir"],
                                            "logs",
                                            "created_log_dir.out"),
        sample_fq = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 "{sample}.fq.gz"),
        script = os.path.join( config["script_dir"],
                               "rs-fastq_to_fasta_awk.sh")
    output:
        sample_fa = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 "{sample}.fa.gz")
    singularity:
        "docker://cjh4zavolab/fastx:0.0.14"
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "fq2fa",
                     "{sample}.log")
    shell:
        '''
        (zcat {input.sample_fq} \
        | {input.script} \
        | fastx_renamer -n COUNT -z \
        > {output.sample_fa})
        2> {log}
        '''

#-------------------------------------------------------------------------------
# get number of raw reads
#-------------------------------------------------------------------------------
rule raw_read_cnt_se:
    ##No Singularity support required##
    input:
        sample_fa = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 "{sample}.fa.gz")
    output:
        raw_cnt = temp(os.path.join(config["samples_dir"],
                                    "counts",
                                    "{sample}.raw.nr.out" ) )
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "counting",
                                   "{sample}.log")
    run:
        import gzip
        n = 0
        with gzip.open(input.sample_fa, "rt") as infile:
            n = sum([1 for line in infile if line.startswith(">")])
        with open(output.raw_cnt, "w") as out:
            out.write("reads.raw.nr\t%i\n" % n)

#-------------------------------------------------------------------------------
# get length of raw reads
#-------------------------------------------------------------------------------
rule raw_read_length_se:
    ##No Singularity support required##
    input:
        raw_cnt = os.path.join(config["samples_dir"],
                               "counts",
                               "{sample}.raw.nr.out" ),
        input_fa = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 "{sample}.fa.gz")
    output:
        raw_len = temp(os.path.join(config["samples_dir"],
                                    "counts",
                                    "{sample}.raw.len.out" ) )
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "counting",
                                   "{sample}.log")
    script:
        os.path.join(config['snakemake_script_dir'],
                     "raw-read-length.py")

#-------------------------------------------------------------------------------
# filter reads without expected 5' start
#-------------------------------------------------------------------------------
rule select_for_valid_5p_configuration:
    input:
        sample_fa = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 "{sample}.fa.gz"),
        script = os.path.join( config["script_dir"],
                               "rs-filter-by-5p-adapter.pl")
    output:
        selected_5p = os.path.join(config["samples_dir"],
                                   "{sample}",
                                   "{sample}.5ptrimmed.A-seq2.fa.gz")
    singularity:
        "docker://cjh4zavolab/select_valid_5p:3"
    params:
        adapt = config['to_trim_from_5p_Aseq2'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "trim_5p_part",
                     "{sample}.log")
    shell:
        '''
        (zcat {input.sample_fa} \
        | perl {input.script} \
        --adapter={params.adapt} \
        | gzip > {output.selected_5p}) 2> {log}
        '''

#-------------------------------------------------------------------------------
# trim 4 nucleotides from the read start (used for 3' READS)
# 3' READS: according to doi:10.1038/nmeth.2288
# each valid read (from rev sequencing which was applied in these samples;
# BUT: new samples have to be checked whether they were still reverse sequenced)
# should have 4 random nt at the 5' end followed by remaining Ts from the
# reverse transcription of the poly(A) tail;
# according to the published protocol, only reads with at least 2 nongenomic
# As were considered valid
# hence, here we select valid 5' configuration as: ....TT and the remaining
# part of the poly(A) tail is trimmed later after reverse complementation
#-------------------------------------------------------------------------------
rule select_for_valid_5p_configuration_3READS:
    input:
        sample_fa = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 "{sample}.fa.gz"),
        script = os.path.join( config["script_dir"],
                               "rs-filter-by-5p-adapter.pl")
    output:
        trimmed_5p = os.path.join(config["samples_dir"],
                                  "{sample}",
                                  "{sample}.5ptrimmed.3READS.fa.gz")
    singularity:
        "docker://cjh4zavolab/select_valid_5p:3"
    params:
        adapt = lambda wildcards: samples.loc[ wildcards.sample, "fiveAdapter"],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "trim_5p_part",
                     "{sample}.log")
    shell:
        '''
        (zcat {input.sample_fa} \
        | perl {input.script} \
        --adapter={params.adapt} \
        | gzip > {output.trimmed_5p}) 2> {log}
        '''

#-------------------------------------------------------------------------------
# trim leading Ts
# (used for samples from PAS-Seq)
#-------------------------------------------------------------------------------
rule trim_leading_Ts:
    input:
        sample_fa = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 "{sample}.fa.gz"),
        script = os.path.join( config["script_dir"],
                               "rs-trim-5p-T.pl")
    output:
        nuc_trimmed = os.path.join(config["samples_dir"],
                                   "{sample}",
                                   "{sample}.leadingTs_trimmed.fa.gz")
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        adapt = "T",
        minLen=config['min_length'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "trim_5p_T",
                     "{sample}.log")
    shell:
        '''
        zcat {input.sample_fa} \
        | perl {input.script} \
        --minLen={params.minLen} \
        --nuc={params.adapt} \
        | gzip > {output.nuc_trimmed}
        '''

#-------------------------------------------------------------------------------
# trim the 5' adapter
# Currently only QuantSeq_REV
#-------------------------------------------------------------------------------
rule trim_5p_adapter_se:
    input:
        in_fa = trim_5p_adapter_input
    output:
        no_5p_adapter = os.path.join(config["samples_dir"],
                                     "{sample}",
                                     "{sample}.5ptrimmed.fa.gz")
    singularity:
        "docker://zavolab/cutadapt:1.16"
    params:
        adapt = lambda wildcards: samples.loc[ wildcards.sample, "fiveAdapter"] if not "*" in samples.loc[ wildcards.sample, "fiveAdapter"] else samples.loc[ wildcards.sample, "fiveAdapter"].split("*")[1],
        minLen=config['min_length'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    resources:
        time = 6
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "trim_5p_adapter",
                     "{sample}.log")
    shell:
        '''
        cutadapt \
        -g {params.adapt} \
        --minimum-length {params.minLen} \
        -o {output.no_5p_adapter} \
        {input.in_fa} \
        &> {log}
        '''


#-------------------------------------------------------------------------------
# trim the 3' adapter
#-------------------------------------------------------------------------------
rule trim_adapter_se:
    input:
        in_fa = trim_adapter_input
    output:
        no_3p_adapter = os.path.join(config["samples_dir"],
                                     "{sample}",
                                     "{sample}.trimmed.fa.gz")
    singularity:
        "docker://zavolab/cutadapt:1.16"
    params:
        adapt = lambda wildcards: samples.loc[ wildcards.sample, "threeAdapter"] if not "*" in samples.loc[ wildcards.sample, "threeAdapter"] else samples.loc[ wildcards.sample, "threeAdapter"].split("*")[1],
        five_p_adapt = lambda wildcards: "" if (pd.isna( samples.loc[ wildcards.sample, "fiveAdapter"]) or (samples.loc[ wildcards.sample, "protocol"] == "3'READS")) else "-g " + samples.loc[ wildcards.sample, "fiveAdapter"],
        minLen=config['min_length'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    resources:
        time = 6
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "trim_3p_adapter",
                     "{sample}.log")
    shell:
        '''
        cutadapt \
        -a {params.adapt} \
        {params.five_p_adapt} \
        --minimum-length {params.minLen} \
        -o {output.no_3p_adapter} \
        {input.in_fa} \
        &> {log}
        '''

#-------------------------------------------------------------------------------
# trim additional nucleotides that might occur between the
# 3' end and the 3' adapter
#-------------------------------------------------------------------------------
rule trim_additional_3p_nuc:
    ##Singularity required: perl##
    input:
        no_3p_adapter = os.path.join(config["samples_dir"],
                                     "{sample}",
                                     "{sample}.trimmed.fa.gz"),
        script = os.path.join( config["script_dir"],
                               "ag-trimm-3p-end.pl")
    output:
        nuc_trimmed = os.path.join(config["samples_dir"],
                                     "{sample}",
                                     "{sample}.trimmed_add_nuc.fa.gz")
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        adapt = lambda wildcards: samples.loc[ wildcards.sample, "threeAdapter"].split("*")[0].rstrip("]").lstrip("["),
        minLen=config['min_length'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "trim_3p_nuc",
                     "{sample}.log")
    shell:
        '''
        zcat {input.no_3p_adapter} \
        | perl {input.script} \
        --minLen={params.minLen} \
        --nuc={params.adapt} \
        | gzip > {output.nuc_trimmed}
        '''

#-------------------------------------------------------------------------------
# reverse complement
#-------------------------------------------------------------------------------
rule reverse_complement:
    input:
       input_seqs = get_reverse_compl_input
    output:
        rev_cmpl = os.path.join(config["samples_dir"],
                                "{sample}",
                                "{sample}.trimmed.rev_cmpl.fa.gz")
    singularity:
        "docker://cjh4zavolab/fastx:0.0.14"
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "rev_cmpl",
                     "{sample}.log")
    shell:
        '''
        zcat {input.input_seqs} \
        | fastx_reverse_complement -z \
        -o {output.rev_cmpl} \
        &> {log}
        '''

#-------------------------------------------------------------------------------
# select valid 3' configuration for 3P-Seq samples
#-------------------------------------------------------------------------------
rule select_for_valid_3p_configuration_3PSeq:
    input:
        in_fa = get_valid_3p_3PSeq_file,
        script = os.path.join( config["script_dir"],
                               "rs-filter-by-3p-adapter.pl")
    output:
        selected_3p = os.path.join(config["samples_dir"],
                               "{sample}",
                               "{sample}.valid_3p_configuration.fa.gz")
    singularity:
        "docker://cjh4zavolab/select_valid_5p:3"
    params:
        adapt = config['to_trim_from_3p_3PSeq'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "trim_3p_part",
                     "{sample}.log")
    shell:
        '''
        (zcat {input.in_fa} \
        | perl {input.script} \
        --adapter={params.adapt} \
        | gzip > {output.selected_3p}) 2> {log}
        '''

#-------------------------------------------------------------------------------
# remove putative leftOvers of the poly(A) tail from the read 3' ends
# 3'READS: 2 As were cleaved already initially as Ts from the 5' end
# 3P-Seq: 2 As were cleaved already initially
#-------------------------------------------------------------------------------
rule remove_polyAtail:
    input:
        no_3p_adapter = get_input_polyAtail_removal
    output:
        no_polyAtail = os.path.join(config["samples_dir"],
                                    "{sample}",
                                    "{sample}.trimmed_tail.fa.gz")
    singularity:
        "docker://zavolab/cutadapt:1.16"
    params:
        adapt = "AAAAAAAAAAAAAA",
        error_rate = config['polyA_trimming_errorRate'],
        minLen=config['min_length'],
        min_overlap = config['polyA_minOverlap'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    resources:
        time = 6
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "trim_polyAtail_adapter",
                     "{sample}.log")
    shell:
        '''
        cutadapt \
        --adapter {params.adapt} \
        --minimum-length {params.minLen} \
        --overlap {params.min_overlap} \
        -e {params.error_rate} \
        -o {output.no_polyAtail} \
        {input.no_3p_adapter} \
        &> {log}
        '''

#-------------------------------------------------------------------------------
# get number of reads after 3' adapter trimming
#-------------------------------------------------------------------------------
rule no_3pAdapter_read_cnt_se:
    ##No Singularity support required##
    input:
        prev_cnt = os.path.join(config["samples_dir"],
                                "counts",
                                "{sample}.raw.len.out" ),
        in_fa = get_reads_after_trimming
    output:
        trimmed_cnt = temp(os.path.join(config["samples_dir"],
                                        "counts",
                                        "{sample}.after_trim.out" ))
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "counting",
                                   "{sample}.log")
    run:
        import gzip
        n = 0
        with gzip.open(input.in_fa, "rt") as infile:
            n = sum([1 for line in infile if line.startswith(">")])
        with open(output.trimmed_cnt, "w") as out:
            with open(input.prev_cnt, "r") as cnt:
                out.write("%s" % cnt.read() )
            out.write("reads.trim.out\t%i\n" % n)

#-------------------------------------------------------------------------------
# collect high confident reads
#-------------------------------------------------------------------------------
rule get_valid_reads:
    ##Singularity needed: perl, zcat, gzip##
    ##Singularity provided: zavolab_minimal:1, not tested##
    '''
    valid reads have:
    not more than 2 Ns
    A-content: maximum 80%
    a 3' nucleotide other than A
    '''
    input:
        valid_rds_in = get_reads_after_trimming,
        script_filter = os.path.join( config["script_dir"],
                                      "ag-filter-seqs-by-nucleotide-composition.pl"),
        script_last = os.path.join( config["script_dir"],
                                    "ag-filter-seqs-by-last-nuc.pl")
    output:
        valid_reads = os.path.join(config["samples_dir"],
                                   "{sample}",
                                   "{sample}.valid.fa.gz")
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        maxN = config['maxN'],
        maxAcontent = config['maxAcontent'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "get_valid_reads",
                     "{sample}.log")
    shell:
        '''
        (zcat {input.valid_rds_in} \
        | perl {input.script_filter} \
        --max {params.maxN} --nuc N \
        | perl {input.script_filter} \
        --max {params.maxAcontent} --nuc A  \
        | perl {input.script_last} \
        | gzip > {output.valid_reads}) 2> {log}
        '''

#-------------------------------------------------------------------------------
# collect high confident reads
#-------------------------------------------------------------------------------
rule get_valid_reads_with_maxLength:
    ##Singularity needed: perl, zcat, gzip##
    ##Singularity provided: zavolab_minimal:1, not tested##
    '''
    valid reads have:
    not more than 2 Ns
    A-content: maximum 80%
    a 3' nucleotide other than A
    a length shorter than a given maximum
    '''
    input:
        valid_rds_in = get_reads_after_trimming,
        script_filter = os.path.join( config["script_dir"],
                                      "ag-filter-seqs-by-nucleotide-composition.pl"),
        script_len_filter = os.path.join( config["script_dir"],
                                          "ag-filter-seqs-by-length.pl"),
        script_last = os.path.join( config["script_dir"],
                                    "ag-filter-seqs-by-last-nuc.pl")
    output:
        valid_reads = os.path.join(config["samples_dir"],
                                   "{sample}",
                                   "{sample}.valid_maxLength.fa.gz")
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        maxLen = lambda wildcards: int(samples.loc[wildcards.sample, "readlen"]) - int(config['min_sense_strand_shortening']),
        maxN = config['maxN'],
        maxAcontent = config['maxAcontent'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "get_valid_reads",
                     "{sample}.log")
    shell:
        '''
        (zcat {input.valid_rds_in} \
        | perl {input.script_len_filter} --max={params.maxLen} \
        | perl {input.script_filter} \
        --max={params.maxN} --nuc=N \
        | perl {input.script_filter} \
        --max={params.maxAcontent} --nuc=A  \
        | perl {input.script_last} \
        | gzip > {output.valid_reads}) 2> {log}
        '''

#-------------------------------------------------------------------------------
# count valid reads
#-------------------------------------------------------------------------------
rule valid_read_cnt:
    ##No Singularity support required##
    '''
    count the reads after filtering for valid read configuration
    '''
    input:
        prev_cnt = os.path.join(config["samples_dir"],
                                "counts",
                                "{sample}.after_trim.out" ),
        in_fa = get_valid_reads
    output:
        valid_cnt = temp(os.path.join(config["samples_dir"],
                                      "counts",
                                      "{sample}.valid.out" ))
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "counting",
                                   "{sample}.log")
    run:
        import gzip
        n = 0
        with gzip.open(input.in_fa, "rt") as infile:
            n = sum([1 for line in infile if line.startswith(">")])
        with open(output.valid_cnt, "w") as out, open(input.prev_cnt, "r") as cnt:
            out.write("%s" % cnt.read() )
            out.write("reads.valid.nr\t%i\n" % n)




################################################################################
# Mapping subpipeline will be called here.
# Currently, unique and multi mappers are returned. This behaviour can be changed
# by setting/removing flag --keep-mm in rule remove_inferiors in the subpipeline.
#
# SEGEMEHL MAPPING
# #
# # Author: adapted from mir-map by Paula Iborra de Toledo
# # Maintainer: christina.herrmann@unibas.ch
# # Date: 2019-05-01
# #
# # This workflow processes appropriate genome and annotation files,
# # performs mapping to genome and transcriptome separately,
# # and finally selects the best mappers.
# #
# # INPUT: transcriptome and genome fasta files, gtf annotation, filtered reads (.fa.gz)
# # OUTPUT: bed.gz of mapped reads, sorted by position
#
# # If used as subworkflow via 'include', don't provide config file!
# #  Configs are specified in config.yaml of main Snakefile!
# # configfile: "segemehl_config.yaml"
# ##################################################################################



#-------------------------------------------------------------------------------
# only consider unique mappers
#-------------------------------------------------------------------------------
rule select_unique_mappers:
    ##Singularity needed: python2##
    ##packages needed: argpase, gzip##
    ## Singularity provided: python:2.7-slim, not tested ##
    input:
        reads_bed = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 config['genome'],
                                 "{sample}.reads.bed.gz"),
        script = os.path.join( config["script_dir"],
                               "rs-select-unique-mappers.py")
    output:
        unique_bed = os.path.join(config["samples_dir"],
                                     "{sample}",
                                     config['genome'],
                                     "{sample}.reads.unique.bed.gz")
    singularity:
        "docker://python:2.7-slim"
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    resources:
        mem = config['unique_mappers.total_RAM']
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "select_unique",
                     "{sample}.log")
    shell:
        '''
        (python {input.script} \
        --bed {input.reads_bed} \
        | gzip > {output.unique_bed}) 2>> {log}
        '''


#-------------------------------------------------------------------------------
# count mapped reads
# (write them in the appropriate file only in the next rule
#-------------------------------------------------------------------------------
rule mapped_read_cnt:
    ##No Singularity needed##
    input:
        prev_cnt = os.path.join(config["samples_dir"],
                                 "counts",
                                 "{sample}.valid.out" ),
        unique_bed = os.path.join(config["samples_dir"],
                                     "{sample}",
                                     config['genome'],
                                     "{sample}.reads.unique.bed.gz"),
        reads_bed = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 config['genome'],
                                 "{sample}.reads.bed.gz")
    output:
        mapped_cnt = temp( os.path.join(config["samples_dir"],
                                  "counts",
                                  "{sample}_" + config['genome'] + ".mapped.out") )
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "counting",
                                   "{sample}.log")
    resources:
        mem = config['mapped_read_nt.total_RAM']
    run:
        import gzip
        unique = 0
        mapped = {}
        with gzip.open(input.reads_bed, "rt") as in_all:
            total_mapped = {line.split("\t")[3]:1 for line in in_all.readlines()}
        with gzip.open(input.unique_bed, "rt") as in_bed:
            unique = sum([1 for line in in_bed])
        multi = len(total_mapped) - unique
        with open(output.mapped_cnt, "w") as out, open(input.prev_cnt, "r") as cnt:
            out.write("%s" % cnt.read() )
            out.write("reads.mapped.uniqueMappers.nr\t%i\n" % unique)
            out.write("reads.mapped.multiMappers.nr\t%i\n" % multi)

#-------------------------------------------------------------------------------
# get 3' ends
#-------------------------------------------------------------------------------
rule get_3p_ends:
    ##Singularity needed: perl##
    ## Singularity provided: zavolab_minimal:1, not tested ##
    '''Only 3' ends with the following characteristics are reported:
    minimum the last 4 nt map perfectly to the genome
    the read was found to be valid before
    '''
    input:
        unique_bed = os.path.join(config["samples_dir"],
                                     "{sample}",
                                     config['genome'],
                                     "{sample}.reads.unique.bed.gz"),
        script = os.path.join( config["script_dir"],
                               "cjh-get-3pEnds-from-bed.pl")
    output:
        end_sites = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 config['genome'],
                                 "{sample}.3pSites.bed.gz")
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        correction = lambda wildcards: "1" if samples.loc[wildcards.sample, "protocol"] == "DRS" else "0",
        exclude_chr = get_excluded_chromosomes,
        min_align = config['min_3p_align'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    resources:
        mem = config['get_3p_ends.total_RAM']
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "get_3p_ends",
                     "{sample}.log")
    shell:
        '''
        (perl {input.script} \
        {params.exclude_chr} \
        --correction={params.correction} \
        --strict \
        --min_align={params.min_align} \
        {input.unique_bed} \
        | gzip > {output.end_sites}) 2>> {log}
        '''

#-------------------------------------------------------------------------------
# count the number of single 3' ends
# The difference to unique mappers are reads
# that don't map perfectly in the last 4 nucleotides
#-------------------------------------------------------------------------------
rule raw_3pSites_cnt:
    ##No Singularity needed##
    input:
        prev_cnt = os.path.join(config["samples_dir"],
                                "counts",
                                "{sample}_" + config['genome'] + ".mapped.out"),
        end_sites = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 config['genome'],
                                 "{sample}.3pSites.bed.gz")
    output:
        sites_cnt = temp( os.path.join(config["samples_dir"],
                                  "counts",
                                  "{sample}_" + config['genome'] + ".raw3pSites.out") )
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "counting",
                                   "{sample}.log")
    run:
        import gzip
        plus = 0
        plus_reads = 0
        minus = 0
        minus_reads = 0
        with gzip.open(input.end_sites, "rt") as in_bed:
            for line in in_bed:
                F = line.rstrip().split("\t")
                if F[5] == "+":
                    plus += 1
                    plus_reads += float(F[4])
                else:
                    minus += 1
                    minus_reads += float(F[4])
        with open(output.sites_cnt, "w") as out, open(input.prev_cnt, "r") as cnt:
            out.write("%s" % cnt.read() )
            out.write("sites.highconfidence.number.plus\t%i\n" % plus)
            out.write("sites.highconfidence.number.minus\t%i\n" % minus)
            out.write("sites.highconfidence.reads.plus\t%i\n" % plus_reads)
            out.write("sites.highconfidence.reads.minus\t%i\n" % minus_reads)

#-------------------------------------------------------------------------------
# extract the sequences that surround the 3' ends
#-------------------------------------------------------------------------------
rule fetch_flanking_seqs:
    ## Singularity available, not tested##
    ## Singularity needed: bedtools 2.27, perl##
    input:
        genome = config["genome_fasta"],
        ends = os.path.join(config["samples_dir"],
                            "{sample}",
                            config['genome'],
                            "{sample}.3pSites.bed.gz"),
        script = os.path.join( config["script_dir"],
                               "rs-fetch-flanking-fasta.pl")
    output:
        seqs = os.path.join(config["samples_dir"],
                            "{sample}",
                            config['genome'],
                            "{sample}.3pSites.bed.seqs.gz")
    singularity:
        "docker://cjh4zavolab/bedtools:2.27"
    params:
        upstream_ext = config['IP.upstream_region'],
        downstream_ext = config['IP.downstream_region'],
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    resources:
        mem = config['fetch_seqs.total_RAM']
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "fetch_flanking_region_seqs",
                     "{sample}.log")
    shell:
        '''
        (perl {input.script} \
        --genome={input.genome} \
        --upstream={params.upstream_ext} \
        --downstream={params.downstream_ext} \
        {input.ends} \
        | gzip > {output.seqs}) 2>> {log}
        '''

#-------------------------------------------------------------------------------
# assign internal priming sites
#-------------------------------------------------------------------------------
rule assign_IP_sites:
    ## Singularity needed: perl##
    ## Singularity provided: zavolab_minimal:1, not tested ##
    input:
        seqs = os.path.join(config["samples_dir"],
                            "{sample}",
                            config['genome'],
                            "{sample}.3pSites.bed.seqs.gz"),
        script = os.path.join( config["script_dir"],
                               "ag-assign-internal-priming-sites.pl")
    output:
        ip_assigned = os.path.join(config["samples_dir"],
                                   "{sample}",
                                   config['genome'],
                                   "{sample}.3pSites.ip.bed.gz")
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        upstream_ext = config['IP.upstream_region'],
        downstream_ext = config['IP.downstream_region'],
        tot_As = config['IP.total_As'],
        consec_As = config['IP.consecutive_As'],
        ds_patterns = get_ds_patterns_for_ipAssignment,
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(config["samples_dir"],
                     "logs",
                     "assign_IP_sites",
                     "{sample}.log")
    shell:
        '''
        (perl {input.script} \
        --upstream_len={params.upstream_ext} \
        --downstream_len={params.downstream_ext} \
        --consecutive_As={params.consec_As} \
        --total_As={params.tot_As} \
        {params.ds_patterns} \
        {input.seqs} \
        | gzip > {output.ip_assigned}) 2>> {log}
        '''

#-------------------------------------------------------------------------------
# count number of IP sites
#-------------------------------------------------------------------------------
rule IP_3pSites_cnt:
    ##No Singularity needed##
    input:
        prev_cnt = os.path.join(config["samples_dir"],
                                 "counts",
                                 "{sample}_" + config['genome'] + ".raw3pSites.out"),
        end_sites = os.path.join(config["samples_dir"],
                                 "{sample}",
                                 config['genome'],
                                 "{sample}.3pSites.ip.bed.gz")
    output:
        ip_cnt = os.path.join(config["samples_dir"],
                              "counts",
                              "{sample}_" + config['genome'] + ".ip3pSites.out")
    params:
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "counting",
                                   "{sample}.log")
    run:
        import gzip
        plus = 0
        plus_reads = 0
        minus = 0
        minus_reads = 0
        with gzip.open(input.end_sites, "rt") as in_bed:
            for line in in_bed:
                F = line.rstrip().split("\t")
                if F[3] == "IP":
                    if F[5] == "+":
                        plus += 1
                        plus_reads += float(F[4])
                    else:
                        minus += 1
                        minus_reads += float(F[4])
        with open(output.ip_cnt, "w") as out, open(input.prev_cnt, "r") as cnt:
            out.write("%s" % cnt.read() )
            out.write("sites.highconfidence.internalpriming.number.plus\t%i\n" % plus)
            out.write("sites.highconfidence.internalpriming.number.minus\t%i\n" % minus)
            out.write("sites.highconfidence.internalpriming.reads.plus\t%i\n" % plus_reads)
            out.write("sites.highconfidence.internalpriming.reads.minus\t%i\n" % minus_reads)


#-------------------------------------------------------------------------------
# Target rule for pre-processing
#-------------------------------------------------------------------------------
rule complete_preprocessing:
    ## LOCAL ##
    input:
        counts = os.path.join(config["samples_dir"],
                              "counts",
                              "{sample}_" + config['genome'] + ".ip3pSites.out")
    output:
        prepro_cmplt = os.path.join(config["samples_dir"],
                                   "{sample}",
                                   config['genome'],
                                   "{sample}.prepro_cmplt.txt")
    shell:
        '''
        echo '#########################\n \
              Pre-processing completed.\n#########################\n \
              Created "{input.counts}"' \
              > {output.prepro_cmplt}
        '''


################################################################################
# combining all samples into the full atlas
################################################################################

#-------------------------------------------------------------------------------
# merge all samples to a full set of 3' end sites
#-------------------------------------------------------------------------------
rule pool_samples:
    ## Singularity needed: perl##
    ## Singularity provided: zavolab_minimal:1, not tested ##
    input:
        dirs_atlas_created = touch(os.path.join(atlas_outputdir,
                                                "logs",
                                                "created_log_dir.out")),
        files = expand( os.path.join(config["samples_dir"],
                                     "{sample}",
                                     config['genome'],
                                     "{sample}.3pSites.ip.bed.gz"),
                        sample = samples.index),
        counts = expand( os.path.join(config["samples_dir"],
                                      "counts",
                                      "{sample}_" + config['genome'] + ".ip3pSites.out"), sample = samples.index),
        script = os.path.join( config["script_dir"],
                               "ag-pool-sites.pl")
    output:
        pooled_sites = os.path.join( atlas_outputdir,
                                     "3pSites.tsv.gz" )
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "pool_samples.log")
    resources:
        mem = lambda wildcards: int( 1.8 * len(samples.index) ),
        time = config['pool_samples.time']
    log:
        os.path.join( atlas_outputdir,
                      "logs",
                      "pool_samples.log")
    shell:
        '''
        (perl {input.script} \
        --noip \
        {input.files} \
        | gzip > {output.pooled_sites}) 2>> {log}
        '''

#-------------------------------------------------------------------------------
# get overall number of unique 3' end sites (without IP sites)
#-------------------------------------------------------------------------------
rule get_unique_3pSites_cnt:
    ##No Singularity needed##
    input:
        pooled_sites = os.path.join( atlas_outputdir,
                                     "3pSites.tsv.gz" )
    output:
        pooled_sites_cnt = os.path.join( atlas_outputdir,
                                         "counts",
                                         "pooled_3p_ends.nr.out" )
    params:
        cluster_log = os.path.join( atlas_outputdir,
                                    "logs",
                                    "cluster_logs",
                                    "counting_sites.log")
    run:
        import gzip
        n = 0
        with gzip.open(input.pooled_sites, "rt") as infile:
            n = sum([1 for line in infile if not line.startswith("#")])
        with open(output.pooled_sites_cnt, "w") as out:
            out.write("3pSites.pooled:\t%i\n" % n)

#-------------------------------------------------------------------------------
# assign poly(A) signals
#-------------------------------------------------------------------------------
rule assign_polyA_signals:
    ##Singularity needed: perl, bedtools2.27##
    ## Singularity provided: bedtools:2.27 ##
    '''
    Assign poly(A) signals to the 3' end sites. Check for signals in the region
    of -60 to +10 around each 3' end site. This region is hardcoded in "ag-assign-polyAsignals.pl".
    NOTE: Order of PAS in column 82 of output file might not be preserved when repeating the run.
    '''
    input:
        pooled_sites = os.path.join( atlas_outputdir,
                                     "3pSites.tsv.gz"),
        script = os.path.join( config["script_dir"],
                               "ag-assign-polyAsignals.pl")
    output:
        sites_with_pas = os.path.join( atlas_outputdir,
                                       "3pSites.PAS.tsv.gz")
    singularity:
        "docker://cjh4zavolab/bedtools:2.27"
    params:
        signals = " ".join(["--motif=%s" % sig for sig in config['polyA_signals'] ]),
        genome = config['genome_fasta'],
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "assign_polyA_signals.log")
    log:
        os.path.join( atlas_outputdir,
                      "logs",
                      "assign_polyA_signals.log")
    shell:
        '''
        (perl {input.script} \
        {params.signals} \
        --genome={params.genome} \
        {input.pooled_sites} \
        | gzip > {output.sites_with_pas}) 2>> {log}
        '''

#-------------------------------------------------------------------------------
# define sample-specific backgrounds
#-------------------------------------------------------------------------------
rule sample_specific_bg:
    ## Singularity needed: perl##
    ## Singularity provided: zavolab_minimal:1, not tested ##
    '''Based on the annotated poly(A) signals,
    iterate over the 3' ends from highest to lowest supported end
    determine the minimum number of reads per 3' end such that
    among all 3' end sites with at least this minimum number of reads
    x % have at least one annotated poly(A) signal
    '''
    input:
        sites_with_pas = os.path.join( atlas_outputdir,
                                       "3pSites.PAS.tsv.gz"),
        script = os.path.join( config["script_dir"],
                               "rs-find-sample-specific-cutoff.pl")
    output:
        sites_filtered = os.path.join( atlas_outputdir,
                                             "filteredSites",
                                             "{sample}.filtered.tsv" )
        #sites_filtered = temp( os.path.join( atlas_outputdir,
        #                                    "filteredSites",
        #                                     "{sample}.filtered.tsv" ) )
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        cutoff = config['sample.BG_polyAsignal_cutoff'],
        upstream_reg = config['sample.BG_upstream_clustering'],
        downstream_reg = config['sample.BG_downstream_clustering'],
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    resources:
        mem = config['sample.BG_total_RAM']
    log:
        os.path.join( atlas_outputdir,
                      "logs",
                      "sample_specific_bg",
                      "{sample}.log")
    shell:
        '''
        perl {input.script} \
        --cutoff={params.cutoff} \
        --upstream={params.upstream_reg} \
        --downstream={params.downstream_reg} \
        --sample={wildcards.sample} \
        {input.sites_with_pas} \
        > {output.sites_filtered} 2>> {log}
        '''

#-------------------------------------------------------------------------------
# merge the sample-specific results to a new overall table of 3' end sites
#-------------------------------------------------------------------------------
rule create_noBG_3pSites_table:
    ##No Singularity needed##
    input:
        filtered = expand( os.path.join( atlas_outputdir,
                                         "filteredSites",
                                         "{sample}.filtered.tsv" ),
                           sample = samples.index),
        raw_table = os.path.join( atlas_outputdir,
                                  "3pSites.PAS.tsv.gz")
    output:
        table_adjusted =  os.path.join(atlas_outputdir,
                                            "3pSites.PAS.filtered.tsv.gz")
        #table_adjusted = temp( os.path.join(atlas_outputdir,
        #                                    "3pSites.PAS.filtered.tsv.gz") )
    params:
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "create_noBG_3pSites_table.log")
    resources:
        mem = lambda wildcards: ( int( len(samples.index) / 50 ) + 1) * 12,
        time = config["noBG_table.time"]
    script:
        os.path.join( config["snakemake_script_dir"],
                      "merge-sample-bg-files-stable.py")

#-------------------------------------------------------------------------------
# delete 3' end sites without cutoff-corrected read support from any sample
#-------------------------------------------------------------------------------
rule delete_noReadSupport_rows:
    ##No Singularity needed##
    input:
        table_adjusted = os.path.join(atlas_outputdir,
                                            "3pSites.PAS.filtered.tsv.gz")
        #table_adjusted = temp( os.path.join(atlas_outputdir,
        #                                    "3pSites.PAS.filtered.tsv.gz") )
    output:
        table_filtered = os.path.join(atlas_outputdir,
                                      "3pSites.PAS.noBG.tsv.gz")
    params:
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "delete_noReadSupport_rows_3pSites_table.log")
    run:
        import gzip
        with gzip.open(output.table_filtered, "wt") as out_file, gzip.open(input.table_adjusted, "rt") as infile:
            for line in infile:
                if line.startswith("#"):
                    out_file.write(line)
                    continue
                line_list = line.rstrip().split("\t")
                read_sum = sum( [1 for i in line_list[3:-2] if float(i) > 0] )
                if read_sum > 0:
                    # this site has still read support
                    out_file.write(line)

#-------------------------------------------------------------------------------
# For each SAMPLE
# get background-corrected number of 3' sites
# get number of sites with PAS
#-------------------------------------------------------------------------------
rule get_noBG_3pSites_per_sample:
    input:
        noBG_sites = os.path.join( atlas_outputdir,
                                     "3pSites.PAS.noBG.tsv.gz" ),
        prev_cnt = os.path.join(config["samples_dir"],
                              "counts",
                              "{sample}_" + config['genome'] + ".ip3pSites.out")
    output:
        noBG_cnt = os.path.join(config["samples_dir"],
                              "counts",
                              "{sample}_" + config['genome'] + ".noBG3pSites.out")
    params:
        sample = "{sample}",
        cluster_log = os.path.join(config["samples_dir"],
                                   "logs",
                                   "cluster_logs",
                                   "counting",
                                   "{sample}.log")
    run:
        import gzip
        sites = 0
        reads = 0
        pas = 0
        pas_reads = 0
        col = 0
        with gzip.open(input.noBG_sites,"rt") as all_sites:
            for line in all_sites:
                if line.startswith("#"):
                    if params.sample in line:
                        F = line.rstrip().split(";")
                    	col = int(F[0].lstrip("#"))
                else:
                    if col == 0:
                        print("Column for sample could not be identified!")
                        print(params.sample)
                        exit()
                    else:
                        line_list = line.rstrip().split("\t")
                        if line_list[col] != "0":
                            sites += 1
                            reads += int(line_list[col])
                            if line_list[-2] != "NA":
                                pas += 1
				pas_reads += int(line_list[col])
        with open(output.noBG_cnt, "w") as out, open(input.prev_cnt, "r") as cnt:
            out.write("%s" % cnt.read() )
            out.write("sites.noBG.all.reads\t%i\n" % reads)
            out.write("sites.noBG.all.number\t%i\n" % sites)
            out.write("sites.noBG.withPAS.reads\t%i\n" % pas_reads)
            out.write("sites.noBG.withPAS.number\t%i\n" % pas)
            if sites != 0:
                out.write("sites.noBG.withPAS.percent\t%i\n" % (pas/sites*100)) # For put in mongo we need int
            else:
                out.write("sites.noBG.withPAS.percent\t%i\n" % sites)
#-------------------------------------------------------------------------------
# For the ATLAS
# get background-corrected number of unique 3' end sites (without IP sites)
#-------------------------------------------------------------------------------
rule get_noBG_3pSites_cnt:
    ##No Singularity needed##
    input:
        prev_cnt = os.path.join( atlas_outputdir,
                                         "counts",
                                         "pooled_3p_ends.nr.out"),
        noBG_sites = os.path.join( atlas_outputdir,
                                     "3pSites.PAS.noBG.tsv.gz" )
    output:
        noBG_sites_cnt = os.path.join( atlas_outputdir,
                                         "counts",
                                         "noBG_3p_ends.nr.out" )
    params:
        cluster_log = os.path.join( atlas_outputdir,
                                    "logs",
                                    "cluster_logs",
                                    "counting_noBG_sites.log")
    run:
        import gzip
        n = 0
        with gzip.open(input.noBG_sites, "rt") as infile:
            n = sum([1 for line in infile if not line.startswith("#")])
        with open(output.noBG_sites_cnt, "w") as out, open(input.prev_cnt, "r") as cnt:
            out.write("%s" % cnt.read() )
            out.write("3pSites.noBG:\t%i\n" % n)

#-------------------------------------------------------------------------------
# cluster individual closely spaced 3' end sites
#-------------------------------------------------------------------------------
rule cluster_sites:
    ##Singularity needed: perl##
    ## Singularity provided: zavolab_minimal:1, not tested ##
    input:
        table_filtered = os.path.join(atlas_outputdir,
                                      "3pSites.PAS.noBG.tsv.gz"),
        script = os.path.join( config["script_dir"],
                                      "ag-generate-clusters.pl")
    output:
        primary_clusters = os.path.join( atlas_outputdir,
                                         "clusters.primary.tsv.gz" )
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        upstream_ext = config['CL.upstream_clustering'],
        downstream_ext = config['CL.downstream_clustering'],
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "cluster_sites.log")
    resources:
        mem = config['CL.total_RAM'],
        time = config['CL.time']
    log:
        os.path.join( atlas_outputdir,
                      "logs",
                      "cluster_sites.log")
    shell:
        '''
        (perl {input.script} \
        --upstream={params.upstream_ext} \
        --downstream={params.downstream_ext} \
        {input.table_filtered} \
        | gzip > {output.primary_clusters}) 2> {log}
        '''

#-------------------------------------------------------------------------------
# get number of primary clusters
#-------------------------------------------------------------------------------
rule get_prim_clusters_cnt:
    ##No Singularity needed##
    input:
        prev_cnt = os.path.join( atlas_outputdir,
                                         "counts",
                                         "noBG_3p_ends.nr.out"),
        clusters = os.path.join( atlas_outputdir,
                                     "clusters.primary.tsv.gz" )
    output:
        clusters_cnt = os.path.join( atlas_outputdir,
                                         "counts",
                                         "prim_clusters.nr.out" )
    params:
        cluster_log = os.path.join( atlas_outputdir,
                                    "logs",
                                    "cluster_logs",
                                    "counting_prim_clusters.log")
    run:
        import gzip
        n = 0
        with gzip.open(input.clusters, "rt") as infile:
            n = sum([1 for line in infile if not line.startswith("#")])
        with open(output.clusters_cnt, "w") as out, open(input.prev_cnt, "r") as cnt:
            out.write("%s" % cnt.read() )
            out.write("clusters.primary:\t%i\n" % n)

#-------------------------------------------------------------------------------
# merge closely spaced clusters
#-------------------------------------------------------------------------------
rule merge_clusters:
    ##Singularity needed: perl##
    ## Singularity provided: zavolab_minimal:1, not tested ##
    ''' ATTENTION:
    The script expects the input file to be formatted
    according to ag-generate-clusters.pl from the A-seq-processing pipeline
    -> all data is accessed with hard coded indices

    cluster are further merged if:
    - an IP candidate has another downstream cluster that shares all its
    PAS with the IP candidate
    - a cluster shares all its PAS with the next cluster upstream
    - two clusters with (?) independent PAS have a combined length smaller the maxsize
    Keep all un-merged clusters without PAS
    '''
    input:
        primary_clusters = os.path.join( atlas_outputdir,
                                         "clusters.primary.tsv.gz"),
        script = os.path.join( config["script_dir"],
                                      "rs-merge-clusters.pl")
    output:
        merged_clusters = os.path.join( atlas_outputdir,
                                         "clusters.merged.tsv.gz")
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    params:
        maxsize = config['CL.max_cluster_size'],
        minDistToPAS = config['CL.min_dist_to_PAS'],
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "merge_clusters.log")
    resources:
        mem = config['CL.merge_RAM'],
        time = config['CL.time']
    log:
        os.path.join( atlas_outputdir,
                      "logs",
                      "merge_clusters.log")
    shell:
        '''
        (perl {input.script} \
        --minDistToPAS={params.minDistToPAS} \
        --maxsize={params.maxsize} \
        {input.primary_clusters} \
        | gzip > {output.merged_clusters}) 2> {log}
        '''

#-------------------------------------------------------------------------------
# get number of merged clusters
#-------------------------------------------------------------------------------
rule get_merged_clusters_cnt:
    ##No Singularity needed##
    input:
        prev_cnt = os.path.join( atlas_outputdir,
                                         "counts",
                                         "prim_clusters.nr.out"),
        clusters = os.path.join( atlas_outputdir,
                                     "clusters.merged.tsv.gz" )
    output:
        clusters_cnt = os.path.join( atlas_outputdir,
                                         "counts",
                                         "merged_clusters.nr.out" )
    params:
        cluster_log = os.path.join( atlas_outputdir,
                                    "logs",
                                    "cluster_logs",
                                    "counting_merged_clusters.log")
    run:
        import gzip
        n = 0
        with gzip.open(input.clusters, "rt") as infile:
            n = sum([1 for line in infile if not line.startswith("#")])
        with open(output.clusters_cnt, "w") as out, open(input.prev_cnt, "r") as cnt:
            out.write("%s" % cnt.read() )
            out.write("clusters.merged:\t%i\n" % n)

#-------------------------------------------------------------------------------
# annotate the location of the clusters with respect to the
# filtered annotation
# annotated are (with the first and last having the highest and lowest
# priority, respectively: TE-terminal exon, EX-exon,
# IN-intron, DS-up to n nt downstream of TE,
# AE-antisense exon, AI-antisense intron, AU-antisense upstream,
# IG-intergenic
#-------------------------------------------------------------------------------
rule annotate_gene_features:
    ##No Singularity needed##
    input:
        merged_clusters = os.path.join( atlas_outputdir,
                                         "clusters.merged.tsv.gz"),
        script = os.path.join( config["script_dir"],
                                      "rs-annotate-gene-features-tsv.py"),
        #anno_filtered = config['gene_annotation_filtered']
        anno = config['gene_annotation']
    output:
        clusters_annotated = os.path.join( atlas_outputdir,
                                               "clusters.anno.tsv.gz")
    params:
        downstream_region = config['ds_range'],
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "annotate_clusters.log")
    log:
        os.path.join( atlas_outputdir,
                      "logs",
                      "annotate_clusters.log")
    shell:
        '''
        python {input.script} \
        --verbose \
        --gtf {input.anno} \
        --ds-range {params.downstream_region} \
        --input {input.merged_clusters} \
        | gzip > {output.clusters_annotated} \
        2> {log}
        '''

#-------------------------------------------------------------------------------
# calculate cluster support measures
#-------------------------------------------------------------------------------
rule cluster_support:
    ''' ATTENTION: This script requires the original design file.
    The correct order of samples must be given.
    '''
    input:
        clusters_annotated = os.path.join( atlas_outputdir,
                                               "clusters.anno.tsv.gz"),
        script = os.path.join( config["script_dir"],
                                      "cjh-get-clusters-support.py")
    output:
        clusters_support = os.path.join( atlas_outputdir,
                                         "clusters.support.tsv.gz"),
        clusters_temp = temp(os.path.join( atlas_outputdir,
                                               "clusters.support.tsv"))
    singularity:
        "docker://python:3.6.9-slim-stretch"
    params:
        design_file = config['atlas.samples_table'],
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "clusters_support.log")
    log:
        os.path.join( atlas_outputdir,
                      "logs",
                      "clusters_support.log")
    shell:
        '''
        python {input.script} \
        --verbose \
        --design={params.design_file} \
        --in {input.clusters_annotated} \
        --out {output.clusters_temp} \
        2> {log} &&
        gzip -c {output.clusters_temp} \
        > {output.clusters_support}
        '''


#-------------------------------------------------------------------------------
# make a bed file with cluster tpms for each sample
#-------------------------------------------------------------------------------
rule make_bed:
    input:
        clusters = os.path.join( atlas_outputdir,
                                               "clusters.support.tsv.gz"),
        script = os.path.join( config["script_dir"],
                                 "cjh-bed-per-sample-from-clusters.py")
    output:
        samples_bed = os.path.join("{path}",
                           "{sample}.clusters.bed.gz"),
        samples_temp = temp(os.path.join("{path}",
                           "{sample}.clusters.bed"))
    singularity:
        "docker://python:3.6.9-slim-stretch"
    params:
        id = "{sample}",
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(atlas_outputdir,
                     "logs",
                     "{sample}.log")

    shell:
        '''
        python {input.script} \
        -i {input.clusters} \
        -s {params.id} \
        -o {output.samples_temp} \
        2> {log} &&
        gzip -c {output.samples_temp} \
        > {output.samples_bed}
        '''

#-------------------------------------------------------------------------------
# sort bed 
#-------------------------------------------------------------------------------
rule sort_bed:
    ## Singularity: bedtools ##
    input:
        bed = os.path.join("{path}",
                           "{sample}.clusters.bed.gz")
    output:
        sorted_bed = os.path.join("{path}",
                                  "{sample}.clusters." + config['ucsc_db'] + "."
                                     + config['atlas.release_name'] + ".bed.gz")
    singularity:
        "docker://cjh4zavolab/bedtools:2.27"
    params:
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(atlas_outputdir,
                     "logs",
                     "{sample}.log")
    shell:
        '''
        sortBed \
        -i {input.bed} \
        | gzip \
        > {output.sorted_bed}
        '''
#-------------------------------------------------------------------------------
# get number of final clusters
# count PAS and annotations
#-------------------------------------------------------------------------------
rule get_final_clusters_stats:
    ##No Singularity needed##
    input:
        prev_cnt = os.path.join( atlas_outputdir,
                                         "counts",
                                         "merged_clusters.nr.out"),
        clusters_bed = os.path.join( atlas_outputdir,
                                         "atlas.clusters." + config['ucsc_db'] + "."
                                         + config['atlas.release_name'] + ".bed.gz")
    output:
        clusters_cnt = os.path.join( atlas_outputdir,
                                         "counts",
                                         "clusters.stats.out" )
    params:
        cluster_log = os.path.join( atlas_outputdir,
                                    "logs",
                                    "cluster_logs",
                                    "clusters_stats.log")
    run:
        import gzip
        n = 0
        p = 0
        annos = {'TE': 0,
                'EX': 0,
                'IN': 0,
                'DS': 0,
                'AE': 0,
                'AI': 0,
                'AU': 0,
                'IG': 0}

        with gzip.open(input.clusters_bed, "rt") as infile:
            for line in infile:
                # Count clusters
                n += 1
                # Count clusters with PAS
                if not "NA" in line:
                    p += 1
                # For each cluster get annotation
                a = line.split('\t')[9]
                annos[a] += 1
        with open(output.clusters_cnt, "w") as out, open(input.prev_cnt, "r") as cnt:
            out.write("{}".format(cnt.read() ))
            out.write("clusters.all:\t{:d}\n".format(n))
            out.write("clusters.PAS.nr:\t{:d}\n".format(p))
            out.write("clusters.PAS.percent:\t{:d}\n".format(int(p/n*100))) # For put in mongo we need int
            for k in annos.keys():
                out.write("clusters.annos.%s:\t%s\n" % (k, annos[k]))



#-------------------------------------------------------------------------------
# Target rule clustering
#-------------------------------------------------------------------------------
rule complete_clustering:
    input:
        noBG_cnt = expand( os.path.join(config["samples_dir"],
                                    "counts",
                                    "{sample}_" + config['genome'] + ".noBG3pSites.out"),
                         sample = samples.index),
        cluster_stats = os.path.join( atlas_outputdir,
                                    "counts",
                                    "clusters.stats.out" ),
        clusters_bed = os.path.join( atlas_outputdir,
                                         "atlas.clusters." + config['ucsc_db'] + "."
                                         + config['atlas.release_name'] + ".bed.gz"),
        samples_bed = expand( os.path.join(config["samples_dir"],
                                     "{sample}",
                                     config['genome'],
                                     "{sample}.clusters." + config['ucsc_db'] + "."
                                     + config['atlas.release_name'] + ".bed.gz"),
                        sample = samples.index)
    output:
        clst_cmplt = os.path.join(atlas_outputdir,
                                    "clst_cmplt.txt")
    shell:
        '''
        echo '#########################\n \
        Clustering completed.\n \
        #########################\n \
        Created "{input.noBG_cnt}"\n \
        "{input.cluster_stats}"\n \
        "{input.clusters_bed}"\n \
        "{input.samples_bed}"\n' \
        > {output.clst_cmplt}
        '''

################################################################################
# Make files for visualization of custom tracks in UCSC genome browser
################################################################################

#-------------------------------------------------------------------------------
# get the UCSC chromosome sizes file
#-------------------------------------------------------------------------------
rule fetch_chr_sizes_ucsc:
    ##LOCAL##
    ##Singularity needed: wget##
    params:
        url = config['ucsc_chromSizes_URL']
    output:
        chr_sizes_ucsc = config["ucsc_chromSizes_file"]
    singularity:
        "docker://cjh4zavolab/zavolab_minimal:1"
    shell:
        '''
        wget -O {output.chr_sizes_ucsc} \
        {params.url} \
        &> /dev/null
        '''

#-------------------------------------------------------------------------------
# prepare bedGraph for bigWig
#-------------------------------------------------------------------------------
rule clusters_bedGraph:
    input:
        clusters = os.path.join( atlas_outputdir,
                                               "clusters.anno.tsv.gz"),
        script = os.path.join( config["script_dir"],
                                 "cjh-bedGraph-from-tsv.py")
    output:
        plus =  os.path.join(atlas_outputdir,
                                     "tracks",
                                     "{sample}_plus.bedGraph"),
        minus =  os.path.join(atlas_outputdir,
                                     "tracks",
                                     "{sample}_minus.bedGraph")
    singularity:
        "docker://python:2.7-slim"
    params:
        id = "{sample}",
        chr_names = lambda wildcards: " ".join([(str(c) +
                                                 ":" +
                                                 config['chromosome_names'][c]) for c in  config['chromosome_names']]),
        cluster_log = os.path.join( atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(atlas_outputdir,
                     "logs",
                     "{sample}.log")
    shell:
        '''
        python {input.script} \
        -i {input.clusters} \
        -s {params.id} \
        --chr-names {params.chr_names} \
        -p {output.plus} \
        -m {output.minus} \
        2> {log}
        '''


#-------------------------------------------------------------------------------
# sort bedGraphs 
#-------------------------------------------------------------------------------
rule sort_bed_4_big:
    ## Singularity: bedtools ##
    input:
        ucsc_bed = os.path.join(atlas_outputdir,
                                     "tracks",
                                     "{sample}_{strand}.bedGraph")
    output:
        sorted_bed = os.path.join(atlas_outputdir,
                                  "tracks",
                                  "{sample}_{strand}.sorted.bedGraph")
    singularity:
        "docker://cjh4zavolab/bedtools:2.27"
    params:
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "{sample}.log")
    log:
        os.path.join(atlas_outputdir,
                     "logs",
                     "{sample}.log")
    shell:
        '''
        sortBed \
        -i {input.ucsc_bed} \
        > {output.sorted_bed}
        '''

#-------------------------------------------------------------------------------
# prepare bigWig from bedGraph
#-------------------------------------------------------------------------------
rule prepare_bigWig:
    ##Singularity needed: bedToBigWig##
    input:
        ucsc_bed = os.path.join(atlas_outputdir,
                                "tracks",
                                "{sample}_{strand}.sorted.bedGraph"),
        chr_sizes = config['ucsc_chromSizes_file']
    output:
        bigWig = os.path.join(atlas_outputdir,
                    "tracks",
                    "{sample}_{strand}." + config['ucsc_db'] + "." + config['atlas.release_name'] + ".bw")
    singularity:
        "docker://kerstenbreuer/kent-bedgraphtobigwig"
    params:
        cluster_log = os.path.join(atlas_outputdir,
                                   "logs",
                                   "cluster_logs",
                                   "bed2bigWig.log")
    log:
        os.path.join(atlas_outputdir,
                     "logs",
                     "bed2bigWig.log")
    shell:
        '''
        bedGraphToBigWig \
        {input.ucsc_bed} \
        {input.chr_sizes} \
        {output.bigWig}
        '''

#-------------------------------------------------------------------------------
# Make accessory files containing track line info
#-------------------------------------------------------------------------------
rule make_track_info:
    ## LOCAL ##
    ##No Singularity needed##
    input:
        bigplus = os.path.join(atlas_outputdir,
                            "tracks",
                            "{sample}_plus." + config['ucsc_db'] + "." + config['atlas.release_name'] + ".bw"),
        bigminus = os.path.join(atlas_outputdir,
                            "tracks", 
                            "{sample}_minus." + config['ucsc_db'] + "." + config['atlas.release_name'] + ".bw")
    params:
        name = "{sample}",
        atlas_public_name = config['ucsc_db'] + "." + config['atlas.release_name'],
        url = os.path.join(config["polyasite_download_url"],
                                   "tracks",
                                   config["genome"]),
        plus = "{sample}_plus." + config['ucsc_db'] + "." + config['atlas.release_name'] + ".bw",
        minus = "{sample}_minus." + config['ucsc_db'] + "." + config['atlas.release_name'] + ".bw",
    output:
        track_info = os.path.join(atlas_outputdir,
                       "tracks",
                      ("{sample}." + config['ucsc_db'] + "." + config['atlas.release_name'] + ".track_info.txt"))

    run:
        with open(output.track_info, "wt") as out:
            out.write('track type=bigWig name="%s: poly(A) clusters plus strand %s" \
                       visibility="full" color="4,177,216" maxHeightPixels="128:60:8"\
                       bigDataUrl="%s/%s"\n'\
                       % (params.name, params.atlas_public_name, params.url, params.plus))
            out.write('track type=bigWig name="%s: poly(A) clusters minus strand %s" \
                       visibility="full" color="241,78,50" maxHeightPixels="128:60:8"\
                       bigDataUrl="%s/%s"\n' \
                       % (params.name, params.atlas_public_name, params.url, params.minus))


#-------------------------------------------------------------------------------
# Target rule tracks
#-------------------------------------------------------------------------------
rule complete_tracks:
    input:
        atlas_track = os.path.join(atlas_outputdir,
                       "tracks",
                      ("atlas." + config['ucsc_db'] + "." + config['atlas.release_name'] + ".track_info.txt")),
        sample_tracks = expand( os.path.join(atlas_outputdir,
                       "tracks",
                      ("{sample}." + config['ucsc_db'] + "." + config['atlas.release_name'] + ".track_info.txt")),
                         sample = samples.index),
    output:
        tracks_cmplt = os.path.join(atlas_outputdir,
                                    "tracks_cmplt.txt")
    shell:
        '''
        echo '#########################\n \
        Track files completed.\n \
        #########################\n \
        Created "{input.atlas_track}"\n \
        "{input.sample_tracks}"\n' \
        > {output.tracks_cmplt}
        '''
#-------------------------------------------------------------------------------
# How did it go?
#-------------------------------------------------------------------------------
onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred, check log at %s." % {log})
