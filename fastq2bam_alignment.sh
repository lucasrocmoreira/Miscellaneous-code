#!/bin/bash
. /seq/vgb/jason/scripts/bashFunctions.sh
source /broad/software/scripts/useuse

#Argument Handling Code; see /seq/vgb/jason/scripts/bashFunctions.sh if more information is required
description="Aligns the provided fastq files and performs preprocessing, namely trimming adapters, MarkDuplicates, IndelRealigner, and, optionally, BQSR.\\nWritten by Jason Turner-Maier (jturner), modified by Lucas R Moreira, but heavily based off of http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently"
argumentDef=()
argumentDef+=("-fastq" fastqFile "The fastq file (read 1) to be aligned / preprocessed" 1)
argumentDef+=("-ref" refFile "The genome reference to align the bam to. Should also be the reference for the knownSites variants, or problems will occur. A sequence dictionary for the reference should also be present in the same folder, as well as a fasta index." 1)
argumentDef+=("-knownVars" knownVariants "The variant file(s) to use for BQSR; if this argument is not passed, BQSR will be skipped." "unlimited" "skip")
argumentDef+=("-outDir" outDir "The directory in which to place output files; will attempt to create it if it does not exist" 1 "./")
argumentDef+=("-tempDir" tempDir "The directory in which to place temporary files; will attempt to create it if it does not exist" 1 "./")
argumentDef+=("-numCPUs" numCPUs "The number of cpus to use for running bwa mem. Note that, due to how the queueing system works, requesting fewer cores will also request less memory. It's probably a good idea to request at least 4 cpus to make sure your job doesn't get killed due to using too much memory." 1 "4")
#argumentDef+=("-queue" ugerQueue "The UGER queue to submit jobs to" 1 "vert") TODO: allow users to specify alternate queues
argumentDef+=("-gatk" gatkVersion "The directory containing the version of GATK (GenomeAnalysisTK.jar) to use" 1 "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-107-g35688f4/")
argumentDef+=("-trimmomatic" trimmomaticVersion "The directory containing the version of Trimmomatic (trimmomatic-0.39.jar) to use" 1 "/seq/vgb/software/Trimmomatic-0.39/")
argumentDef+=("-picard" picardVersion "The directory containing the version of picard (picard.jar) to use" 1 "/seq/software/picard/1.999/bin/")
argumentDef+=("-bwa" bwaVersion "The directory containing the version of bwa to use" 1 "/seq/vgb/software/bwa/bwa-0.7.17/bwa")
argumentDef+=("-fastqc" fastqcVersion "The directory containing the version of fastqc to use" 1 "/seq/vgb/lmoreira/software/FastQC/fastqc")

if [ "$#" == "0" ]; then
    argumentsArr=("-h")
else
    argumentsArr=("$@")
fi
handleArguments "$description" argumentDef[@] argumentsArr[@]

#Set up dotkits and variables
# reuse BamTools
reuse .bamtools-2.5.1
reuse GCC-5.2
reuse Java-1.8
reuse Samtools
reuse UGER
if [ ! -d $outDir ]; then
    mkdir $outDir
fi
if [ ! -d $tempDir ]; then
    mkdir $tempDir
fi

if [ "$knownVariants" != "skip" ]; then
    knownVariants=`echo $knownVariants | sed -e 's: : --knownSites :g' -e 's:^:--knownSites :'`
fi

currID=`basename $fastqFile | sed 's:.R1_001.fastq.gz::'`
read1=$fastqFile
read2=`echo $fastqFile | sed 's/R1/R2/'`

#Get read group information
header=`zcat $read1 | head -1`
IFS=':' read -a header <<< "$header"
INSTRUMENT=${header[0]}
RUN_ID=${header[1]}
FLOWCELL_BARCODE=${header[2]}
LANE=${header[3]}
ID=$FLOWCELL_BARCODE.$LANE
PU=$FLOWCELL_BARCODE.$LANE.$currID
SM=$currID
PL=ILLUMINA
LB=TrueSeq

#Process based on pipeline description at https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently#step1
#$currID.1.* files relate to running Trimmomatic, which removes Illumina adaptors
qsub -l h_vmem=32g -l h_rt=36:00:00 -b y -p -10 -N trim.$currID -cwd -o $outDir/$currID.1.trim.out -j y -V java -Xmx32G -jar $trimmomaticVersion/trimmomatic-0.39.jar PE $read1 $read2 $outDir/$currID.1.forward_paired.fq.gz $outDir/$currID.1.forward_unpaired.fq.gz $outDir/$currID.1.reverse_paired.fq.gz $outDir/$currID.1.reverse_unpaired.fq.gz ILLUMINACLIP:$trimmomaticVersion/adapters/TruSeq3-PE-2.fa:2:30:10:8:true

#$currID.2.* files relate to performing QC
qsub -l h_vmem=4g -l h_rt=6:00:00 -b y -p -10 -N fastqc.$currID -cwd -o $outDir/$currID.2.fastqc.out -j y -V -hold_jid trim.$currID $fastqcVersion $outDir/$currID.1.forward_paired.fq.gz $outDir/$currID.1.reverse_paired.fq.gz

#$currID.3.* files relate to performing the alignment
#We write the commands to a script so that we can bwa, and MergeBamAlignment piped into each other as below
#This enables us to not have to write out large temporary read files
echo "#!/bin/bash" > $outDir/$currID.3.align.sh
echo "$bwaVersion/bwa mem -M -t $numCPUs -p $refFile $outDir/$currID.1.forward_paired.fq.gz $outDir/$currID.1.reverse_paired.fq.gz | samtools sort -o $outDir/$currID.bam -@ $numCPUs" >> $outDir/$currID.3.align.sh
chmod u=rwx $outDir/$currID.3.align.sh
qsub -l h_rt=48:00:00 -p -10 -pe smp $numCPUs -binding linear:$numCPUs -S /bin/bash -N aln.$currID -cwd -o $outDir/$currID.3.align.out -j y -V -hold_jid fastqc.$currID -l h_vmem=6G $outDir/$currID.3.align.sh

#$currID.4.* files relate to adding read groups
qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N readgroup.$currID -cwd -o $outDir/$currID.4.readgroup.out -j y -V -hold_jid aln.$currID java -Xmx32G -jar $picardVersion/picard.jar AddOrReplaceReadGroups INPUT=$outDir/$currID.bam OUTPUT=$outDir/$currID.4.groups_added.bam RGID=$ID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM

#$currID.5.* files relate to running MarkDuplicates, which tries to find and mark redundant reads
qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N mdp.$currID -cwd -o $outDir/$currID.5.mark_dup.out -pe smp $numCPUs -binding linear:$numCPUs -j y -V -hold_jid readgroup.$currID java -Xmx12G -jar $picardVersion/picard.jar MarkDuplicates INPUT=$outDir/$currID.4.groups_added.bam OUTPUT=$outDir/$currID.5.mark_dup.bam METRICS_FILE=$outDir/$currID.5.mark_dup.metrics TMP_DIR=$tempDir
qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N ind.$currID -cwd -o $outDir/$currID.5.index.out -j y -V -hold_jid mdp.$currID samtools index $outDir/$currID.5.mark_dup.bam

#$currID.6.* files relate to running IndelRealigner, which fixes read alignments around indels
#qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N rtc.$currID -cwd -o $outDir/$currID.6.indel_targets.out -j y -V -hold_jid ind.$currID java -Xmx32g -Djava.io.tmpdir=$tempDir -jar $gatkVersion/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $refFile -I $outDir/$currID.5.mark_dup.bam -o $outDir/$currID.6.indel_targets.intervals
#qsub -l h_vmem=64g -l h_rt=20:00:00 -b y -p -10 -N idr.$currID -cwd -o $outDir/$currID.6.indel_realigned.out -j y -V -hold_jid rtc.$currID java -Xmx64g -Djava.io.tmpdir=$tempDir -jar $gatkVersion/GenomeAnalysisTK.jar -T IndelRealigner -R $refFile -I $outDir/$currID.5.mark_dup.bam -targetIntervals $outDir/$currID.6.indel_targets.intervals -o $outDir/$currID.6.indel_realigned.bam

#$currID.7.* files relate to running BaseQualityScoreRecalibration, which modifies base qualities to make them more accurate for performing variant calling
#Since our variant calling pipeline includes BQSR, we don't always want to run it; thus it is opt-in, argument-wise
if [ "$knownVariants" != "skip" ]; then
    qsub -l h_vmem=32g -l h_rt=48:00:00 -b y -p -10 -N brc.$currID -cwd -o $outDir/$currID.7.base_recal_targets.out -j y -V -hold_jid idr.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T BaseRecalibrator -R $refFile -I $outDir/$currID.6.indel_realigned.bam -o $outDir/$currID.7.bqsr_results.table $knownVariants
    qsub -l h_vmem=32g -l h_rt=48:00:00 -b y -p -10 -N brp.$currID -cwd -o $outDir/$currID.7.base_recal_execute.out -j y -V -hold_jid brc.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T PrintReads -R $refFile -I $outDir/$currID.6.indel_realigned.bam -BQSR $outDir/$currID.7.bqsr_results.table -o $outDir/$currID.7.base_recalibrated.bam
    qsub -l h_vmem=32g -l h_rt=48:00:00 -b y -p -10 -N brt.$currID -cwd -o $outDir/$currID.7.base_recal_targets_second_pass.out -j y -V -hold_jid brc.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T BaseRecalibrator -R $refFile -I $outDir/$currID.6.indel_realigned.bam -BQSR $outDir/$currID.7.bqsr_results.table -o $outDir/$currID.7.bqsr_results.second_pass.table $knownVariants
    reuse R-3.2
    qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N bqc.$currID -cwd -o $outDir/$currID.7.base_recal_qc.out -j y -V -hold_jid brt.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $refFile -before $outDir/$currID.7.bqsr_results.table -after $outDir/$currID.7.bqsr_results.second_pass.table -plots $outDir/$currID.7.bqsr_QC.pdf
fi
