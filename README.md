# PolyASite Workflow

Pipeline to infer poly(A) site clusters through processing of 3' end sequencing libraries prepared according
to various protocols. The pipeline was used for the generation of the
[PolyASite](http://polyasite.unibas.ch/) atlas.

## Pipeline schematic
Further information on the implemented processing can be found [below](#protocol-specific-notes) and on [PolyASite](http://polyasite.unibas.ch/).
![Pipeline DAG]( ./pipeline_dag.png )


## Requirements

* [Snakemake v4.4.0+](https://snakemake.readthedocs.io/en/stable/)
* [Singularity v2.4.4+](https://singularity.lbl.gov/)

> The pipeline was tested on an HPC environment managed by
> [Slurm](https://slurm.schedmd.com/).

### Conda environment

We recommend that to use a Conda environment that contains the necessary
software.

The environment was created with:

```bash
conda env create \
      --name polyA_atlas_pipeline \
      --file snakemake_run_env_requirements.yaml
```

Activate the environment with:

```bash
source activate polyA_atlas_pipeline
```

Deactivate the environment with:

```bash
source deactivate
```

## Prerequisites

The pipeline has three central elements: the **Snakefile**, a **config file**,
and a **sample table**:

* The **Snakefile** does not need to be changed, modified or updated unless
  bugs, intended updates etc., require changes.
* The config file (called `config.yaml`) needs to be adjusted according to your
  needs. It requires a set of **samples**, an **organism**, a **genome version**
  and an **annotation version**. Organism and genome/annotation versions are
  edited directly in the config file. Samples are provided indirectly by
  indicating the path to the **sample table**. An example for a config file can
  be found here: `tests/EXAMPLE_config.yaml`.  
  Please see below in the **Run section** part how you can also run Snakemake
  without modifying the config file, instead specifying the required information
  as arguments to command line parameters.
* The sample table lists sample-specific information required for the successful
  run of the pipeline. Follow the format in `tests/EXAMPLE_samples.tsv` and
  include one line for each sample to be processed.

## Run the pipeline

### Local

Local execution is **not recommended** and should only be used for testing
purposes.

```bash
snakemake \
    -p \
    --use-singularity \
    --singularity-args "--bind ${PWD}" \
    --configfile config.yaml \
    &>> run_update.Organism_genomeVersion_annotationVersion.log
```

### Slurm

```bash
snakemake \
    -p \
    --use-singularity \
    --singularity-args "--bind ${PWD}" \
    --configfile config.yaml \
    --cluster-config cluster_config.json \
    --jobscript jobscript.sh \
    --cores 500 \
    --local-cores 10 \
    --cluster "sbatch --cpus-per-task {cluster.threads} \
              --mem {cluster.mem} --qos {cluster.queue} \
              --time {cluster.time} -o {params.cluster_log} \
              -p [NODE_ID] --export=JOB_NAME={rule} \
              --open-mode=append" \
    &>> run_update.Organism_genomeVersion_annotationVersion.log
```

### General notes on running the pipeline

#### Instant config changes

If you do not want to change the base config file, you can also specify the
appropriate values in the Snakemake command itself, e.g.:

```bash
snakemake \
    -p \
    --use-singularity \
    --singularity-args "--bind ${PWD}" \
    --configfile config.yaml \
    --config organism=HomoSapiens genome=GRCh38.90 atlas.release_name:r2.0 -- \
    &>> run_update.Organism_genomeVersion_annotationVersion.log
```
#### Running only parts of the pipeline
With the Snakemake option `--until` you can specify a target rule for running the pipeline. 
This is useful if you only want to run a subset of rules. For example, the pipeline includes 
the creation of bigWig- and track-info files for display of data in the UCSC genome browser. 
If you don't need these files, run with `--until complete_clustering`. For examples see 
[here](./tests/EXAMPLE_run_clustering.sh).

#### Generate graph for preprocessing of individual sample of specific protocol

When a target file is provided in the Snakemake command, Snakemake will run the
pipeline only until the point at which the desired file is created. This can be
used to generate overview graphs on the processing of samples from specific
protocols to check their processing steps. For example:

```bash
snakemake \
    -p \
    --configfile config.yaml \
    --dag \
    samples/counts/SRX517313_GRCh38.90.ip3pSites.out \
    | dot \
        -T png \
    > graph.SAPAS.png
```

### Protocol-specific notes

#### 3' READS

Pre-processing involves:

1. Filtering of reads based on the 5' configuration  
2. 3' adapter trimming
3. Reverse complementing
4. 3' trimming of potentially remaining `A`s from the poly(A) tail

> Only reads that start with a specified number of random nucleotides and two
> `T`s are considered, with the number of random nucleotides having to be
> extracted from the corresponding GEO/SRA entries or publication. See
> [here](https://doi.org/10.1371/journal.pgen.1005166) for more information.
>  
> According to [this paper](https://doi.org/10.1038/nmeth.2288), sequencing can
> be done in sense and antisense direction. The samples that are currently
> processed here were sequenced in antisense direction. Future samples should be
> checked carefully in order to decide whether current settings are appropriate.

#### SAPAS

Pre-processing involves:

1. Combined 5' and 3' adapter trimming
2. Trimming of remaining Cs at the 3' end (they are a result of template
   switching reverse transcription)
3. Reverse complementing

> Sequencing libraries are prepared such that sequencing can be done either in
> antisense direction (Illumina) or in sense direction (454). So far, only
> samples from Illumina sequencing have been processed. If "sense direction
> samples" need to be processed, the pipeline **must** be adapted accordingly.
>  
> In the supplementary material of the
> [APASdb paper](http://dx.doi.org/10.1093/nar/gku1076), the authors state that
> they only consider reads that have the expected 5' linker sequence
> `5'-TTTTCTTTTTTCTTTTTT-3'`. However, manual comparison of the first reads from
> an old sample (`SRX026584`) with one from the mentioned publication revealed
> that the abundance of this linker is very low in the samples from the APASdb
> publication. Therefore, we have decided not to be too strict with the 5' linker.
> Very often, only poly(T)s are found: therefore, for the newer samples only a
> stretch of Ts is given as 5' adapter and trimmed together with the 3' adapter.

#### A-seq2

Processing is done according to
[this protocol](http://dx.doi.org/10.3791/56129), but without using the first
nucleotide as barcode information.

#### A-seq

Pre-processing involves:

1. 3' adapter trimming
2. Valid reads filtering with additional consideration of a maximum read length
   (see below).

#### 3'-seq (Mayr)

Pre-processing involves:

1. 3' adapter trimming
2. Additional trimming of As and Ns at the 3' end
3. Valid reads filtering with additional consideration of a maximum read length
   (see below).

#### PAS-Seq

Pre-processing involves:

1. Trimming of poly(T) at the 5' end
2. 3' adapter trimming
3. Trimming of additional Cs at the 3' end
4. Reverse complementing

#### DRS

Current pre-processing involves:

1. Reverse complementing
2. Correction by 1 nt to obtain the true 3' end position (note that this is
   directly encoded in the Snakefile, *not* the config file.)

> The protocol facilitates direct sequencing of the RNA 3' end. Due to an
> initial `T`-fill step that involves the incorporation of a blocking
> nucleotide (anything except for a `T`), sequencing begins actually one
> nucleotide upstream of the RNA 3'-most nucleotide. Therefore, a correction of
> 1 in the downstream direction of the read's reverse complement is necessary to
> obtain the 3' end (See "Extended Experimental Procedures" in
> [Ozsolak et al](http://dx.doi.org/10.1016/j.cell.2010.11.020) for more info.

#### PolyA-seq

Pre-processing involves:

1. 3' adapter trimming
2. Reverse complementing

#### 3P-Seq

Pre-processing involves:

1. Reverse complementing if necessary
2. Filtering reads: only proceed with reads with at least 2 As at the 3' end
3. Remove additional `A`s at the 3' end that might remain from the poly(A) tail

> Note that processing for samples of this protocol is **sample-specific**. In
> particular, only a subset of samples requires reverse complementing, and
> hence, each sample has to be checked manually to infer whether reverse
> complementing is required or not. Once, this information is provided in the
> design file, the pipeline will process samples accordingly.
>  
> An easy way to check whether a file needs to be revese complemented is to
> count the occurence of the poly(A) signal in the first reads and their revese
> complements:
>  
> ```bash
> zcat samples/GSM1268942/GSM1268942.fa.gz \
>     | head -n10000 \
>     | tail -n1000 \
>     | grep TTTATT | wc -l
> zcat samples/GSM1268942/GSM1268942.fa.gz \
>     | head -n10000 \
>     | tail -n1000 \
>     | grep AATAAA | wc -l
> ```
>  
> Comparing these numbers should give a clear prefernce for one of the two signals.

#### 2P-Seq

Pre-processing involves:

1. 3' adapter trimming
2. Reverse complementing

#### Maximum read length

Samples prepared with protocols **3'-seq (Mayr)** and **A-seq** have a
restriction on the maximum read length for processed reads to count as _valid_.
As these protocols require sequencing in the sense direction, the length
restriction ensures that the 3' end of the transcript is reached.
