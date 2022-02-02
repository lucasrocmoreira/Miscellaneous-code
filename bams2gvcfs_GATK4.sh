#!/bin/bash
#%V     SCRIPTNAME 0.1
#%V     Copyright (C) 2021 Lucas Moreira

source /broad/software/scripts/useuse
reuse Java-1.8
use UGER

GATK_EXECUTABLE=/seq/vgb/software/gatk/gatk-4.2.2.0/gatk
tmpdir=/home/unix/lmoreira/vgb/gorilla/tmp

Help()
{
   # Display Help
   echo
   echo "Script to call genotypes with GATK4."
   echo
   echo "Syntax: bams2gvcf_GATK4.sh REF BAM_FILE"
   echo "REF = reference genome file (.fasta)"
   echo "BAM_FILE = aligned bam file (.bam)"
   echo
}

# If argument is provided, print help message
# $@ here means all the arguments of the script ($1 and $2)
if [[ $@ ]]; then
    true
else
    Help
    exit 1
fi

file=`basename $2`
name=`echo $file | cut -d '.' -f1`

if [ ! -d GVCF/ ]; then
  mkdir GVCF/;
fi

export _JAVA_OPTIONS=-Djava.io.tmpdir=$tmpdir

#chroms=(CM017847.1 CM017848.1 CM017849.1 CM017850.1 CM017851.1 CM017852.1 CM017853.1 CM017854.1 CM017855.1 CM017856.1 CM017857.1 CM017858.1 CM017859.1 CM017860.1 CM017861.1 CM017862.1 CM017863.1 CM017864.1 CM017865.1 CM017866.1 CM017867.1 CM017868.1 CM017869.1 CM017870.1 chrU)

# -b tells the cluster to execute the job as a binary. This saves a lot of time because the cluster will not transfer your script but instead just the pointer to the path.
# -p defines the priority of the job
# -V declares that all environment variables in the qsub command's environment	are to be exported to the batch job.
# -j declares if the standard error stream of the job will be merged with the standard output stream of 
# -cwd tells that the job should be executed in the same directory that qsub was called.
qsub -l h_vmem=32g -l h_rt=48:00:00 -b y -p -10 -N $name -cwd -o GVCF/$name.haplotype_caller.ERC.out -j y -V $GATK_EXECUTABLE --java-options "-Xmx10g" HaplotypeCaller --tmp-dir $tmpdir -R $1 -I $2 -O GVCF/$name.g.vcf.gz -ERC GVCF -GQB 20 -GQB 100 --min-pruning 2 #-GQB specifies the exclusive upper bounds for reference confidence GQ bands
