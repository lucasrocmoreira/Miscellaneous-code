#!/bin/bash
source /broad/software/scripts/useuse
reuse UGER
reuse Java-1.8

#for i in /home/unix/lmoreira/vgb/gorilla/fq/*R1.fastq.gz; do 
for i in /home/unix/lmoreira/vgb/gorilla/fq/gor1_R1.fastq.gz; do
    read1=$i;
    file1=`echo $read1 | cut -d '_' -f1`;
    name=`echo $file1 | cut -d '/' -f8`;
    read2ending="_R2.fastq.gz";
    read2=$file1$read2ending;
    header=`zcat $read1 | head -1`;
    echo $read1 $read2;
    IFS=':' read -a header <<< "$header"
    INSTRUMENT=${header[0]}
    RUN_ID=${header[1]}
    FLOWCELL_BARCODE=${header[2]}
    LANE=${header[3]}
    ID=$FLOWCELL_BARCODE.$LANE
    PU=$FLOWCELL_BARCODE.$LANE.$name
    SM=$name
    PL=ILLUMINA
    LB=HiSeq2000
    echo $name $read1 $read2 $ID $SM $LB $PU $PL

    java -Xmx8G -jar /home/unix/lmoreira/vgb/software/picard.jar FastqToSam FASTQ=$read1 FASTQ2=$read2 OUTPUT=$file1.unaligned.bam READ_GROUP_NAME=$ID SAMPLE_NAME=$SM LIBRARY_NAME=$LB PLATFORM_UNIT=$PU PLATFORM=$PL

done
