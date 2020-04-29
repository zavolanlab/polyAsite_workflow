# set -e

snakemake \
-np \
-s ../Snakefile \
--use-singularity \
--singularity-args "--bind '${PWD}/../'" \
--configfile EXAMPLE_config.yaml \
--cluster-config ../cluster_config.json \
--jobscript ../jobscript.sh \
--cores 500 \
--local-cores 10 \
--cluster "sbatch \
	--cpus-per-task {cluster.threads} \
	--mem {cluster.mem} \
	--qos {cluster.queue} \
	--time {cluster.time} \
	-o {params.cluster_log} \
	-p scicore \
	--export=JOB_NAME={rule} \
	--open-mode=append" \
--resources load=100 \
--reason \
--until complete_preprocessing \
#&>> preprocessing.log
