#!/bin/bash
date
if [ ${#1} == 0 ]
then
	echo Required positional inputs: filtered enhancers, variant file list, motif file list, output name
	exit
fi

ENH_BED=$1
VAR=$2
MOT=$3
OUT_NAME=$4

VAR_FILES=($(cat $VAR))
MOT_FILES=($(cat $MOT))

module load bedtools

for M in ${MOT_FILES[@]}
do
	date
	echo $M
	awk 'BEGIN{FS=OFS="\t"};$2>0{print $1,$2-1,$3,$4,$5,$6,$7}' $M | intersectBed -a - -b $ENH_BED -wa -wb >> motif_enhancer_intersect.temp.bed
done

sort -k1,1 -k2,2n motif_enhancer_intersect.temp.bed > motif_enhancer_intersect_sorted.temp.bed
rm motif_enhancer_intersect.temp.bed

MOTIF_ENH_HEADER="chr\tmotif_start\tmotif_end\tmotif_name\tmotif_score\tmotif_strand\tmotif_file\tchr\tenh_start\tenh_end\tgene_id\tgene_name\tgene_strand\tcell_lines"

date
for V in ${VAR_FILES[@]}
do
	V_HEADER=$(head -1 $V | sed 's/\t/\\t/g')
	COMBINED_HEADER=$(echo "$MOTIF_ENH_HEADER\t$V_HEADER")
	case $V in
		*snp*) SAVE_NAME=$OUT_NAME\_snp.bed ;;
		*indel*) SAVE_NAME=$OUT_NAME\_indel.bed ;;
		*) SAVE_NAME=$OUT_NAME\_$V.bed ;;
	esac
	echo -e $COMBINED_HEADER > $SAVE_NAME
	intersectBed -a motif_enhancer_intersect_sorted.temp.bed -b $V -wa -wb >> $SAVE_NAME
done

rm motif_enhancer_intersect_sorted.temp.bed

date

