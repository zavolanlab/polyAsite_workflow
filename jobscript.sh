#!/bin/bash
# properties = {properties}

#PATH=/scicore/home/zavolan/schmiral/soft/miniconda3/bin:${{PATH}}
#export PATH

echo -e "JOB ID\t$SLURM_JOBID"
echo "=============================="
echo -e "rule\t$JOB_NAME"
echo -e "=============================="

ml purge
ml Singularity
{exec_job}

echo -e "==============================\n"
