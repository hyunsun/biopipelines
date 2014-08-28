#!/bin/bash
# $1 - file_name-sample_id mapping file path

source ./pipeline.cfg

cluster_name=$2
while read sample_id; do
  echo "[$(date "+%F %T")] Running sample ${sample_id}..."
  echo ""
  
  qsub -pe make 8 -cwd -j y ./pipeline_bam.sh ${sample_id} ${cluster_name}
done < $1
