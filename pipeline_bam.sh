#!/bin/bash
#$ -V
#$ -S /bin/bash

source ./pipeline.cfg

# arguments
category=$1
file_name=$2
sample_id=$3

# log function
time_start=$(date +%s)
time_now=${time_start}
time_end=${time_start}
LOG_FILE="./logs/${sample_id}.log"
LOG() {
  time_now=$(date +%s)
  echo "[${sample_id}][$(date "+%F %T")][time: $((${time_now} - ${time_end})) sec] $*" >> ${LOG_FILE}
  time_end=${time_now}
}

# functions 
check_error() {
if [ $1 -ne 0 ]; then
  echo "ERROR: pipeline.sh"
  echo "ERROR CODE: $1"
fi
}

check_line_size_bam() {
  fastqR1linesize=`zcat $1 | wc -l | cut -d' ' -f1`
  check_error $?

  fastqhalfsize=`expr ${fastqR1linesize} / 2`
  check_error $?

  bamlinesize=`${samtools} view $2 | wc -l ${FASTQ} | cut -d' ' -f1`
  check_error $?

  LOG "$1 half line size: ${fastqhalfsize}"
  LOG "$2 line size: ${bamlinesize}"
  if [ ${fastqhalfsize} -ne ${bamlinesize} ]; then
    LOG "check_line_size_bam failed!"
  fi
}

# variables
fastq_r1="${work_dir}fastq/disk5/${file_name}_1.fastq"
fastq_r2="${work_dir}fastq/disk5/${file_name}_2.fastq"
output="${work_dir}results/${category}/bam/${sample_id}"
stat="${work_dir}results/${category}/stat/${sample_id}"
fastqc="${work_dir}results/${category}/fastqc/"

LOG "Start processing sample ${sample_id}..."

# 0 Stat
(/SNUH/app/FastQC/fastqc -t 4 -o ${fastqc} -f fastq ${fastq_r1})&
(/SNUH/app/FastQC/fastqc -t 4 -o ${fastqc} -f fastq ${fastq_r2})&
wait
LOG "#0 FastQC done."

mv ${fastqc}${file_name}"_1_fastqc.zip" ${fastqc}${sample_id}"_1_fastqc.zip"
mv ${fastqc}${file_name}"_2_fastqc.zip" ${fastqc}${sample_id}"_2_fastqc.zip"

# 1
echo "#1"
echo "${bwa} mem -M -t 8 ${ref} ${fastq_r1} ${fastq_r2} > ${output}.sam"
${bwa} mem -M -t 8 ${ref} ${fastq_r1} ${fastq_r2} > ${output}.sam
LOG "#1 BWA done."

# 2
echo "#2"
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}AddOrReplaceReadGroups.jar INPUT=${output}.sam OUTPUT=${output}.sorted.bam SORT_ORDER=coordinate RGID=${output} RGLB=${output} RGPL=illumina RGPU=SureSelectAllExon RGSM=${output} VALIDATION_STRINGENCY=LENIENT"
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}AddOrReplaceReadGroups.jar INPUT=${output}.sam OUTPUT=${output}.sorted.bam SORT_ORDER=coordinate RGID=${sample_id} RGLB=${sample_id} RGPL=illumina RGPU=SureSelectAllExon RGSM=${sample_id} VALIDATION_STRINGENCY=LENIENT
LOG "#2 AddOrReplaceReadGroups done."
#check_line_size_bam ${fastq_r1} ${output}.sorted.bam

# 3
echo "#3"
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}MarkDuplicates.jar INPUT=${output}.sorted.bam OUTPUT=${output}.sorted.dp.bam METRICS_FILE=${output}.sorted.dp.metrix ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}MarkDuplicates.jar INPUT=${output}.sorted.bam OUTPUT=${output}.sorted.dp.bam METRICS_FILE=${output}.sorted.dp.metrix ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
LOG "#3 MarkDuplicates done."
  
# 4 Stat
echo "#4"
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T FlagStat -R ${ref} -L ${bait} -I ${output}.sorted.dp.bam -o ${stat}.stat"
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T FlagStat -R ${ref} -L ${bait} -I ${output}.sorted.dp.bam -o ${stat}.stat
LOG "#4 FlagStat done."

# 5
echo "#5"
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T RealignerTargetCreator -nt 8 -R ${ref} -I ${output}.sorted.dp.bam -o ${output}.intervals -known ${mills} -known ${KG}"
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T RealignerTargetCreator -nt 8 -R ${ref} -I ${output}.sorted.dp.bam -o ${output}.intervals -known ${mills} -known ${KG}
LOG "#5 RealignerTargetCreator done."

# 6
echo "#6"
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T IndelRealigner -R ${ref} -I ${output}.sorted.dp.bam -o ${output}.sorted.dp.ir.bam -targetIntervals ${output}.intervals  -known ${mills} -known ${KG}"
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T IndelRealigner -R ${ref} -I ${output}.sorted.dp.bam -o ${output}.sorted.dp.ir.bam -targetIntervals ${output}.intervals  -known ${mills} -known ${KG}
LOG "#6 IndelRealigner done."

# 7
echo "#7"
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T BaseRecalibrator -nct 8 -I ${output}.sorted.dp.ir.bam -R ${ref} -knownSites ${dbsnp} -knownSites ${mills} -knownSites ${KG} -o ${output}.sorted.dp.ir.grp"
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T BaseRecalibrator -nct 8 -I ${output}.sorted.dp.ir.bam -R ${ref} -knownSites ${dbsnp} -knownSites ${mills} -knownSites ${KG} -o ${output}.sorted.dp.ir.grp
LOG "#7 BaseRecalibrator done."

# 8
echo "#8"
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T PrintReads -nct 8 -R ${ref} -I ${output}.sorted.dp.ir.bam -BQSR ${output}.sorted.dp.ir.grp -o ${output}.recal.bam"
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T PrintReads -nct 8 -R ${ref} -I ${output}.sorted.dp.ir.bam -BQSR ${output}.sorted.dp.ir.grp -o ${output}.recal.bam
LOG "#8 PrintReads done."

# 9 Stat
echo "#9"
echo "java -Xmx8g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T DepthOfCoverage -R ${ref} -I ${output}.recal.bam -L ${bait} -omitBaseOutput -omitLocusTable -omitIntervals -o ${stat}.cov.out --minMappingQuality 20 --minBaseQuality 20 --logging_level ERROR --summaryCoverageThreshold 8 --summaryCoverageThreshold 15 --summaryCoverageThreshold 30 --summaryCoverageThreshold 50"
(java -Xmx8g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}CollectInsertSizeMetrics.jar VALIDATION_STRINGENCY=LENIENT INPUT=${output}.recal.bam OUTPUT=${stat}.stat.insertsizematrix HISTOGRAM_FILE=${stat}.stat.hist.pdf)&
(java -Xmx8g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T DepthOfCoverage -R ${ref} -I ${output}.recal.bam -L ${bait} -omitBaseOutput -omitLocusTable -omitIntervals -o ${stat}.cov.out --minMappingQuality 20 --minBaseQuality 20 --logging_level ERROR --summaryCoverageThreshold 8 --summaryCoverageThreshold 15 --summaryCoverageThreshold 30 --summaryCoverageThreshold 50)&
wait
LOG "#9 DepthOfCoverage done."

# 10 upload to g-Storage

LOG "Done sample ${sample_id}, total elasped time: $((${time_now} - ${time_start})) sec."

# cleanup
#rm $output.sam && rm $output.sorted.bam
#rm $output.sorted.dp.bam && rm $output.sorted.dp.ir.bam && rm $output.sorted.dp.ir.grp 
