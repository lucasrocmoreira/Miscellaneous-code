#!/bin/bash
. /seq/vgb/jason/scripts/bashFunctions.sh
source /broad/software/scripts/useuse

#Argument Handling Code; see /seq/vgb/jason/scripts/bashFunctions.sh if more information is required
description="Realigns the provided bam file and performs preprocessing, namely MarkDuplicates, IndelRealigner, and, optionally, BQSR.\\nWritten by Jason Turner-Maier (jturner), but heavily based off of http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently"
argumentDef=()
argumentDef+=("-bam" bamFile "The bam file to be realigned / preprocessed" 1)
argumentDef+=("-ref" refFile "The genome reference to align the bam to. Should also be the reference for the knownSites variants, or problems will occur. A sequence dictionary for the reference should also be present in the same folder, as well as a fasta index." 1)
argumentDef+=("-knownVars" knownVariants "The variant file(s) to use for BQSR; if this argument is not passed, BQSR will be skipped." "unlimited" "skip")
argumentDef+=("-outDir" outDir "The directory in which to place output files; will attempt to create it if it does not exist" 1 "./")
argumentDef+=("-tempDir" tempDir "The directory in which to place temporary files; will attempt to create it if it does not exist" 1 "./")
argumentDef+=("-numCPUs" numCPUs "The number of cpus to use for running bwa mem. Note that, due to how the queueing system works, requesting fewer cores will also request less memory. It's probably a good idea to request at least 4 cpus to make sure your job doesn't get killed due to using too much memory." 1 "4")
argumentDef+=("-retainBams" retainIBams "If passed, intermediate bams are retained instead of being deleted." 0 "no")
#argumentDef+=("-queue" ugerQueue "The UGER queue to submit jobs to" 1 "vert") TODO: allow users to specify alternate queues
argumentDef+=("-gatk" gatkVersion "The directory containing the version of GATK (GenomeAnalysisTK.jar) to use" 1 "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-107-g35688f4/")
argumentDef+=("-picard" picardVersion "The directory containing the version of picard (picard.jar) to use" 1 "/seq/software/picard/1.999/bin/")
argumentDef+=("-bwa" bwaVersion "The directory containing the version of bwa to use" 1 "/seq/vgb/software/bwa/bwa-0.7.17/bwa")

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

currID=`basename $bamFile | sed 's:.bam::'`

#Read group processing
#echo "The processing pipeline requires the bam file's read group to have, at least, an ID tag, a platform tag, and a sample tag. Here are the read group(s) from the input bam file:"
#bamtools header -in $bamFile | grep '@RG'
#echo "If the read group(s) have all three of the above tags (ID, PL, and SM), enter \"yes\". If there is exactly one read group and it is missing tags, enter \"no\". If tags are missing, but there are more than one read group, enter \"multi\"."
#read -p "> " rgResults
#while [ $rgResults != "yes" ] && [ $rgResults != "no" ] && [ $rgResults != "multi" ]; do
#    echo "Please enter \"yes\", \"no\", or \"multi\"."
#    read -p "> " rgResults
#done
#
#if [ $rgResults == "multi" ]; then
#    echo "This script is not currently configured to correctly update bams with missing information and multiple read groups. Please fix the read groups on your own and rerun. Sorry!"
#    exit 0
#elif [ $rgResults == "no" ]; then
#    read -p "Please enter the correct ID: " rgID
#    read -p "Please enter the correct platform: " rgPL
#    read -p "Please enter the correct sample: " rgSM
#    read -p "Please enter the correct library (or Unknown): " rgLB
#    read -p "Please enter the correct platform unit (or Unknown): " rgPU
#    qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N frg.$currID -cwd -o $outDir/$currID.0.fix_read_group.out -j y -V java -Xmx8G -jar $picardVersion/picard.jar AddOrReplaceReadGroups I=$bamFile O=$outDir/$currID.0.read_group_fixed.bam RGID=$rgID RGPL=$rgPL RGLB=$rgLB RGPU=$rgPU RGSM=$rgSM VALIDATION_STRINGENCY=LENIENT
#    qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N fin.$currID -cwd -o $outDir/$currID.0.fix_read_group_index.out -j y -V -hold_jid frg.$currID samtools index $outDir/$currID.0.read_group_fixed.bam
#    bamFile=$outDir/$currID.0.read_group_fixed.bam
#fi

#Process based on pipeline description at http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
#$currID.1.* files relate to running RevertSam, which clears various information out of the bam file
qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N rev.$currID -cwd -o $outDir/$currID.1.revert.out -j y -V -hold_jid fin.$currID java -Xmx8G -jar $picardVersion/picard.jar RevertSam I=$bamFile O=$outDir/$currID.1.revert.bam SANITIZE=true MAX_DISCARD_FRACTION=0.005 ATTRIBUTE_TO_CLEAR=XT ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=AS ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=true REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempDir

#$currID.2.* files relate to marking illumina adapters, which does what it says on the can
qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N mrk.$currID -cwd -o $outDir/$currID.2.mark.out -j y -V -hold_jid rev.$currID java -Xmx8G -jar $picardVersion/picard.jar MarkIlluminaAdapters I=$outDir/$currID.1.revert.bam O=$outDir/$currID.2.mark.bam M=$outDir/$currID.2.mark.metrics.txt TMP_DIR=$tempDir

#$currID.3.* files relate to performing the alignment
#We write the commands to a script so that we can run SamToFastq, bwa, and MergeBamAlignment piped into each other as below
#This enables us to not have to write out large temporary read files
echo "#!/bin/bash" > $outDir/$currID.3.align.sh
echo "java -Xmx32G -jar $picardVersion/picard.jar SamToFastq I=$outDir/$currID.2.mark.bam FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=$tempDir | $bwaVersion/bwa mem -M -t $numCPUs -p $refFile /dev/stdin | java -Xmx32G -jar $picardVersion/picard.jar MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=$outDir/$currID.1.revert.bam OUTPUT=$outDir/$currID.3.aligned.bam R=$refFile CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS TMP_DIR=$tempDir" >> $outDir/$currID.3.align.sh
if [ "$retainIBams" == "no" ]; then
    echo "if [ -f $outDir/$currID.0.read_group_fixed.bam ]; then" >> $outDir/$currID.3.align.sh
    echo "    rm $outDir/$currID.0.read_group_fixed.ba[im]" >> $outDir/$currID.3.align.sh
    echo "fi" >> $outDir/$currID.3.align.sh
    echo "rm $outDir/$currID.1.revert.ba[im]" >> $outDir/$currID.3.align.sh
    echo "rm $outDir/$currID.2.mark.ba[im]" >> $outDir/$currID.3.align.sh
fi
chmod u=rwx $outDir/$currID.3.align.sh
qsub -l h_rt=48:00:00 -p -10 -pe smp $numCPUs -binding linear:$numCPUs -S /bin/bash -N aln.$currID -cwd -o $outDir/$currID.3.align.out -j y -V -hold_jid mrk.$currID -l h_vmem=6G $outDir/$currID.3.align.sh

#$currID.4.* files relate to running MarkDuplicates, which tries to find and mark redundant reads
qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N mdp.$currID -cwd -o $outDir/$currID.4.mark_dup.out -pe smp $numCPUs -binding linear:$numCPUs -j y -V -hold_jid aln.$currID java -Xmx12G -jar $picardVersion/picard.jar MarkDuplicates INPUT=$outDir/$currID.3.aligned.bam OUTPUT=$outDir/$currID.4.mark_dup.bam METRICS_FILE=$outDir/$currID.4.mark_dup.metrics TMP_DIR=$tempDir
qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N ind.$currID -cwd -o $outDir/$currID.4.index.out -j y -V -hold_jid mdp.$currID samtools index $outDir/$currID.4.mark_dup.bam

#$currID.5.* files relate to running IndelRealigner, which fixes read alignments around indels
qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N rtc.$currID -cwd -o $outDir/$currID.5.indel_targets.out -j y -V -hold_jid ind.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $refFile -I $outDir/$currID.4.mark_dup.bam -o $outDir/$currID.5.indel_targets.intervals
qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N idr.$currID -cwd -o $outDir/$currID.5.indel_realigned.out -j y -V -hold_jid rtc.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T IndelRealigner -R $refFile -I $outDir/$currID.4.mark_dup.bam -targetIntervals $outDir/$currID.5.indel_targets.intervals -o $outDir/$currID.5.indel_realigned.bam

#$currID.6.* files relate to running BaseQualityScoreRecalibration, which modifies base qualities to make them more accurate for performing variant calling
#Since our variant calling pipeline includes BQSR, we don't always want to run it; thus it is opt-in, argument-wise
if [ "$knownVariants" != "skip" ]; then
    qsub -l h_vmem=32g -l h_rt=48:00:00 -b y -p -10 -N brc.$currID -cwd -o $outDir/$currID.6.base_recal_targets.out -j y -V -hold_jid idr.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T BaseRecalibrator -R $refFile -I $outDir/$currID.5.indel_realigned.bam -o $outDir/$currID.6.bqsr_results.table $knownVariants
    qsub -l h_vmem=32g -l h_rt=48:00:00 -b y -p -10 -N brp.$currID -cwd -o $outDir/$currID.6.base_recal_execute.out -j y -V -hold_jid brc.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T PrintReads -R $refFile -I $outDir/$currID.5.indel_realigned.bam -BQSR $outDir/$currID.6.bqsr_results.table -o $outDir/$currID.6.base_recalibrated.bam
    qsub -l h_vmem=32g -l h_rt=48:00:00 -b y -p -10 -N brt.$currID -cwd -o $outDir/$currID.6.base_recal_targets_second_pass.out -j y -V -hold_jid brc.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T BaseRecalibrator -R $refFile -I $outDir/$currID.5.indel_realigned.bam -BQSR $outDir/$currID.6.bqsr_results.table -o $outDir/$currID.6.bqsr_results.second_pass.table $knownVariants
    reuse R-3.2
    qsub -l h_vmem=32g -l h_rt=20:00:00 -b y -p -10 -N bqc.$currID -cwd -o $outDir/$currID.6.base_recal_qc.out -j y -V -hold_jid brt.$currID java -Xmx4g -jar $gatkVersion/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $refFile -before $outDir/$currID.6.bqsr_results.table -after $outDir/$currID.6.bqsr_results.second_pass.table -plots $outDir/$currID.6.bqsr_QC.pdf
fi
    
if [ "$retainIBams" == "no" ]; then
    echo "rm $outDir/$currID.3.aligned.ba[im]" > $outDir/$currID.omega.cleanup.sh
    echo "rm $outDir/$currID.4.mark_dup.ba[im]" >> $outDir/$currID.omega.cleanup.sh
    if [ "$knownVariants" != "skip" ]; then
	echo "rm $outDir/$currID.5.indel_realigned.bam" >> $outDir/$currID.omega.cleanup.sh
    fi
    chmod u=rwx $outDir/$currID.omega.cleanup.sh
    qsub -l h_vmem=32g -l h_rt=2:00:00 -p -10 -S /bin/bash -N cln.$currID -cwd -o $outDir/$currID.omega.cleanup.out -j y -V -hold_jid idr.$currID,bqc.$currID -l h_vmem=6G $outDir/$currID.omega.cleanup.sh
fi
