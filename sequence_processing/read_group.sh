#!/bin/bash
# Run with: bash read_group.sh [ref] [file_r1] [file_r2] [sample] [output]

HEADER=$(zcat $2 | head -n 1)
REF=$1
R1=$2
R2=$3
OUT=$5

FLWCL=$(echo $HEADER | head -n 1 | cut -d':' -f3)
LN=$(echo $HEADER | head -n 1 | cut -d':' -f4)
ADPT=$(echo $HEADER | head -n 1 | cut -d':' -f10)

ID=$(echo "$FLWCL.$LN")
PL="illumina"
PU=$(echo "$FLWCL.$LN.$ADPT")
LB=$(echo "$ADPT")
SM=$4

#printf "Read Group @RG\tID:${ID}\tPL:${PL}\tPU:${PU}\tLB:${LB}\tSM:${SM}"

bwa mem -t 8 -R $(echo "@RG\tID:${ID}\tPL:${PL}\tPU:${PU}\tLB:${LB}\tSM:${SM}") $REF $R1 $R2 | samtools sort -o $OUT
