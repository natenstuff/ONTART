#!/bin/bash

fastqFile=$1
samples=("hs1" "mm39" "rheMac10")
minMapSize=800
path=/home/tnsmits/lustre/ONTART/minION/analysis

source /home/tnsmits/programs/miniconda3/etc/profile.d/conda.sh

rm -r finalTry
mkdir finalTry


conda activate seqkitEnv
seqkit fq2fa ${fastqFile} -o "${path}/finalTry/ontArt.reads.fa"
conda deactivate

conda activate samtoolsEnv
samtools faidx "${path}/finalTry/ontArt.reads.fa"
conda deactivate

##for every reference (==sample): prepare the readBed of the primary alignment
for sample in "${samples[@]}"; do

    echo "dealing with: ${sample}"
##map the fastq to each genome
    conda activate minimapEnv
    minimap2 -a -x map-ont -t 16 "/home/tnsmits/lustre/lib/${sample}.fa" "${fastqFile}" --secondary=no > "${path}/finalTry/${sample}.sam"
    conda deactivate

##samtools to keep only the primary alignment with mapQ==60
    conda activate samtoolsEnv
    samtools view -S -b -F 256  -F 0x800 -q 60 "${path}/finalTry/${sample}.sam" > "${path}/finalTry/${sample}.filtered.bam"
    conda deactivate

##bedtools to keep only the alignment mapping 800 bp or more
    conda activate bedtoolsEnv
    bedtools bamtobed -i "${path}/finalTry/${sample}.filtered.bam" > "${path}/finalTry/${sample}.filtered.bed"
    conda deactivate

    cat "${path}/finalTry/${sample}.filtered.bed" | awk -v size="${minMapSize}" '{if (($3 - $2) >= size ) { print $0}}' OFS='\t' > "${path}/finalTry/${sample}.filtered.${minMapSize}bp.bed"
    cat "${path}/finalTry/${sample}.filtered.${minMapSize}bp.bed" | awk '{ print $4}' OFS='\t' > "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.txt"

#Get the total length of every read
    cat "${path}/finalTry/ontArt.reads.fa.fai" | fgrep -f "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.txt" | cut -f1-2 > "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.Length.tsv"

#SAMTOOLS to create readSegments
    conda activate samtoolsEnv
    samtools view -F16 "${path}/finalTry/${sample}.filtered.bam" | awk '{ print $1, $6, "+"}' OFS='\t' |  sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\([0-9]*\)\([SH]\)\([0-9DMI]*[DMI]\)\([0-9]*\)\([SH]\)\([[:blank:]]*\)\(.\)/\1\2\3\2\6\8\9/' | sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\([0-9DMI]*[DMI]\)\([0-9]*\)\([SH]\)\([[:blank:]]*\)\(.\)/\1\20\2\4\6\7/' | sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\([0-9]*\)\([SH]\)\([0-9DMI]*[DMI]\)\([[:blank:]]*\)\(.\)/\1\2\3\20\6\7/' | sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\([0-9DMI]*[DMI]\)\([[:blank:]]*\)\(.\)/\1\20\20\4\5/' | fgrep -f "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.txt" > "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.readSegment.tsv"
    samtools view -f16 "${path}/finalTry/${sample}.filtered.bam" | awk '{ print $1, $6, "-"}' OFS='\t' |  sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\([0-9]*\)\([SH]\)\([0-9DMI]*[DMI]\)\([0-9]*\)\([SH]\)\([[:blank:]]*\)\(.\)/\1\2\6\2\3\8\9/' | sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\([0-9DMI]*[DMI]\)\([0-9]*\)\([SH]\)\([[:blank:]]*\)\(.\)/\1\2\4\20\2\6\7/' | sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\([0-9]*\)\([SH]\)\([0-9DMI]*[DMI]\)\([[:blank:]]*\)\(.\)/\1\20\2\3\2\6\7/' | sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\([0-9DMI]*[DMI]\)\([[:blank:]]*\)\(.\)/\1\20\20\4\5/'  |  fgrep -f "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.txt" >> "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.readSegment.tsv"

#Create readBeds
    join -t $'\t'  <(sort "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.readSegment.tsv" ) <(sort "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.Length.tsv") | awk '{ print $1, $2, ($5 - $3), $4  }' OFS='\t' > "${path}/finalTry/${sample}.readBed.bed"

#Get the alignmentScores for primary alignments
    samtools view  "${path}/finalTry/${sample}.filtered.bam" | awk '{ print $1, $14}' OFS='\t' | fgrep -f "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.txt" | sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\(AS:i:\)\([0-9]*\)/\1\2\4/' > "${path}/finalTry/${sample}.filtered.${minMapSize}bp.readNames.Scores.tsv"
    conda deactivate

#Create badReads based on alignments to reference gaps, mittochondrial and unknown DNA sections
    conda activate bedtoolsEnv
    bedtools intersect -wa -a "${path}/finalTry/${sample}.filtered.${minMapSize}bp.bed" -b "/home/tnsmits/lustre/lib/genomeGaps/${sample}Gaps.extended.bed" | awk '{print $4}' OFS='\t' >> "${path}/finalTry/badReads.txt"
    conda deactivate
    cat  "${path}/finalTry/${sample}.filtered.${minMapSize}bp.bed" | fgrep "chrM" | awk '{print $4}' OFS='\t' >> "${path}/finalTry/badReads.txt"
    cat  "${path}/finalTry/${sample}.filtered.${minMapSize}bp.bed" | fgrep "chrU" | awk '{print $4}' OFS='\t' >> "${path}/finalTry/badReads.txt"

done


#Identify those reads that map to at least two genomes
join -t $'\t'  <(sort "./finalTry/hs1.readBed.bed"  ) <(sort "./finalTry/mm39.readBed.bed" )     | fgrep -v -f "./finalTry/badReads.txt" > "./finalTry/temp.doubleMap.readBed.hs1.mm39.bed"
join -t $'\t'  <(sort "./finalTry/hs1.readBed.bed"  ) <(sort "./finalTry/rheMac10.readBed.bed" ) | fgrep -v -f "./finalTry/badReads.txt" > "./finalTry/temp.doubleMap.readBed.hs1.rheMac10.bed"
join -t $'\t'  <(sort "./finalTry/mm39.readBed.bed" ) <(sort "./finalTry/rheMac10.readBed.bed" ) | fgrep -v -f "./finalTry/badReads.txt" > "./finalTry/temp.doubleMap.readBed.mm39.rheMac10.bed"

#Identify those reads that map to three genomes
join -t $'\t'  <(sort "./finalTry/temp.doubleMap.readBed.hs1.mm39.bed") <(sort "./finalTry/rheMac10.readBed.bed" ) > "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.bed"
#cat "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.bed" | awk '{print $1}' OFS='\t' > "./finalTry/tripleMap.readNames.triple.hs1.mm39.rheMac10.txt"

#Identify those reads that map to just two genomes by removing triple mappers
cat "./finalTry/temp.doubleMap.readBed.hs1.mm39.bed"     | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.bed" ) > "./finalTry/doubleMap.readBed.hs1.mm39.bed"
cat "./finalTry/temp.doubleMap.readBed.hs1.rheMac10.bed" | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.bed" ) > "./finalTry/doubleMap.readBed.hs1.rheMac10.bed"
cat "./finalTry/temp.doubleMap.readBed.mm39.rheMac10.bed"| fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.bed" ) > "./finalTry/doubleMap.readBed.mm39.rheMac10.bed"

#Calculate the individual mapLengths and interval or overlapSizes for double mappers
cat "./finalTry/doubleMap.readBed.hs1.mm39.bed"      | awk '{ print $1, $2, $3, $4, ($3 - $2), $5, $6, $7, ($6 - $5)}' OFS='\t' | awk '{ if ($2 < $6) { if ($3 > ($7 - 1)) {print $0, ($7-$6)} else { print $0, ($3-$6)}} else { if ($7 > ($3 - 1)){ print $0, ($3-$2)} else { print $0, ($7-$2)}}}' OFS='\t' > "./finalTry/doubleMap.readBed.hs1.mm39.length.overlapSize.bed"
cat "./finalTry/doubleMap.readBed.hs1.rheMac10.bed"  | awk '{ print $1, $2, $3, $4, ($3 - $2), $5, $6, $7, ($6 - $5)}' OFS='\t' | awk '{ if ($2 < $6) { if ($3 > ($7 - 1)) {print $0, ($7-$6)} else { print $0, ($3-$6)}} else { if ($7 > ($3 - 1)){ print $0, ($3-$2)} else { print $0, ($7-$2)}}}' OFS='\t' > "./finalTry/doubleMap.readBed.hs1.rheMac10.length.overlapSize.bed"
cat "./finalTry/doubleMap.readBed.mm39.rheMac10.bed" | awk '{ print $1, $2, $3, $4, ($3 - $2), $5, $6, $7, ($6 - $5)}' OFS='\t' | awk '{ if ($2 < $6) { if ($3 > ($7 - 1)) {print $0, ($7-$6)} else { print $0, ($3-$6)}} else { if ($7 > ($3 - 1)){ print $0, ($3-$2)} else { print $0, ($7-$2)}}}' OFS='\t' > "./finalTry/doubleMap.readBed.mm39.rheMac10.length.overlapSize.bed"

#Reads where the overlap >= 85% of the mapLength of one of the alignments, are counted as a single mapper of unknownOrigin
cat "./finalTry/doubleMap.readBed.hs1.mm39.length.overlapSize.bed"      | awk '{ if ($5 <= $9) { if ( ($5 * 0.85) <=  $10 ) {print $1}} else { if (($9 * 0.85) <= $10) { print $1}}}'  OFS='\t' > "./finalTry/unknownOrigin.doubleMap.readName.hs1.mm39.txt"
cat "./finalTry/doubleMap.readBed.hs1.rheMac10.length.overlapSize.bed"  | awk '{ if ($5 <= $9) { if ( ($5 * 0.85) <=  $10 ) {print $1}} else { if (($9 * 0.85) <= $10) { print $1}}}'  OFS='\t' > "./finalTry/unknownOrigin.doubleMap.readName.hs1.rheMac10.txt"
cat "./finalTry/doubleMap.readBed.mm39.rheMac10.length.overlapSize.bed" | awk '{ if ($5 <= $9) { if ( ($5 * 0.85) <=  $10 ) {print $1}} else { if (($9 * 0.85) <= $10) { print $1}}}'  OFS='\t' > "./finalTry/unknownOrigin.doubleMap.readName.mm39.rheMac10.txt"

#Remove singleOverlappers from doubleMappers
cat "./finalTry/doubleMap.readBed.hs1.mm39.length.overlapSize.bed"      | fgrep -v -f "./finalTry/unknownOrigin.doubleMap.readName.hs1.mm39.txt"      > "./finalTry/ARTEFACT.doubleMap.readBed.hs1.mm39.length.overlapSize.bed"
cat "./finalTry/doubleMap.readBed.hs1.rheMac10.length.overlapSize.bed"  | fgrep -v -f "./finalTry/unknownOrigin.doubleMap.readName.hs1.rheMac10.txt"  > "./finalTry/ARTEFACT.doubleMap.readBed.hs1.rheMac10.length.overlapSize.bed"
cat "./finalTry/doubleMap.readBed.mm39.rheMac10.length.overlapSize.bed" | fgrep -v -f "./finalTry/unknownOrigin.doubleMap.readName.mm39.rheMac10.txt" > "./finalTry/ARTEFACT.doubleMap.readBed.mm39.rheMac10.length.overlapSize.bed"

#Add scores to singleOverlappers to determine their origin
join -t $'\t' <(sort "./finalTry/unknownOrigin.doubleMap.readName.hs1.mm39.txt")      <(sort  <( join -t $'\t' <(sort "./finalTry/hs1.filtered.${minMapSize}bp.readNames.Scores.tsv")  <(sort   "./finalTry/mm39.filtered.${minMapSize}bp.readNames.Scores.tsv")))     > "./finalTry/unknownOrigin.hs1.mm39.readName.Scores.tsv"
join -t $'\t' <(sort "./finalTry/unknownOrigin.doubleMap.readName.hs1.rheMac10.txt")  <(sort  <( join -t $'\t' <(sort "./finalTry/hs1.filtered.${minMapSize}bp.readNames.Scores.tsv")  <(sort   "./finalTry/rheMac10.filtered.${minMapSize}bp.readNames.Scores.tsv"))) > "./finalTry/unknownOrigin.hs1.rheMac10.readName.Scores.tsv"
join -t $'\t' <(sort "./finalTry/unknownOrigin.doubleMap.readName.mm39.rheMac10.txt") <(sort  <( join -t $'\t' <(sort "./finalTry/mm39.filtered.${minMapSize}bp.readNames.Scores.tsv") <(sort   "./finalTry/rheMac10.filtered.${minMapSize}bp.readNames.Scores.tsv"))) > "./finalTry/unknownOrigin.mm39.rheMac10.readName.Scores.tsv"

#Assign singleOverlappers to genome of origin
cat "./finalTry/unknownOrigin.hs1.mm39.readName.Scores.tsv"      | awk '{if ($2 > $3) {print $1}}' OFS='\t'  > "./finalTry/unknownOrigin.hs1.mm39.readName.HUMAN.txt"
cat "./finalTry/unknownOrigin.hs1.mm39.readName.Scores.tsv"      | awk '{if ($2 < $3) {print $1}}' OFS='\t'  > "./finalTry/unknownOrigin.hs1.mm39.readName.MOUSE.txt"
cat "./finalTry/unknownOrigin.hs1.mm39.readName.Scores.tsv"      | awk '{if ($2 == $3) {print $1}}' OFS='\t' > "./finalTry/unknownOrigin.hs1.mm39.readName.AMBIG.txt"

cat "./finalTry/unknownOrigin.hs1.rheMac10.readName.Scores.tsv"  | awk '{if ($2 > $3) {print $1}}' OFS='\t'  > "./finalTry/unknownOrigin.hs1.rheMac10.readName.HUMAN.txt"
cat "./finalTry/unknownOrigin.hs1.rheMac10.readName.Scores.tsv"  | awk '{if ($2 < $3) {print $1}}' OFS='\t'  > "./finalTry/unknownOrigin.hs1.rheMac10.readName.RHEMAC.txt"
cat "./finalTry/unknownOrigin.hs1.rheMac10.readName.Scores.tsv"  | awk '{if ($2 == $3) {print $1}}' OFS='\t' > "./finalTry/unknownOrigin.hs1.rheMac10.readName.AMBIG.txt"

cat "./finalTry/unknownOrigin.mm39.rheMac10.readName.Scores.tsv" | awk '{if ($2 > $3) {print $1}}' OFS='\t'  > "./finalTry/unknownOrigin.mm39.rheMac10.readName.MOUSE.txt"
cat "./finalTry/unknownOrigin.mm39.rheMac10.readName.Scores.tsv" | awk '{if ($2 < $3) {print $1}}' OFS='\t'  > "./finalTry/unknownOrigin.mm39.rheMac10.readName.RHEMAC.txt"
cat "./finalTry/unknownOrigin.mm39.rheMac10.readName.Scores.tsv" | awk '{if ($2 == $3) {print $1}}' OFS='\t' > "./finalTry/unknownOrigin.mm39.rheMac10.readName.AMBIG.txt"

#Calculate overlapSizes for tripleMappers
cat "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.bed" | awk '{ print $1, $2, $3, $4, ($3 - $2), $5, $6, $7, ($6 - $5), $8, $9, $10, ($9 - $8)}' OFS='\t' | awk '{ if ($2 < $6) { if ($3 > ($7 - 1)) {print $0, ($7-$6)} else { print $0, ($3-$6)}} else { if ($7 > ($3 - 1)){ print $0, ($3-$2)} else { print $0, ($7-$2)}}}' OFS='\t' | awk '{ if ($6 < $10) { if ($7 > ($11 - 1)) {print $0, ($11-$10)} else { print $0, ($7-$10)}} else { if ($11 > ($7 - 1)){ print $0, ($7-$6)} else { print $0, ($11-$6)}}}' OFS='\t' | awk '{ if ($2 < $10) { if ($3 > ($11 - 1)) {print $0, ($11-$10)} else { print $0, ($3-$10)}} else { if ($11 > ($3 - 1)){ print $0, ($3-$2)} else { print $0, ($11-$2)}}}' OFS='\t' | awk '{ if ($14 <= $15) { if ($14 <= $16) {print $0, $14} else { print $0, $16}} else { if ($15 <= $16) {print $0, $15} else {print $0, $16}}}' OFS='\t' | awk '{ if ($14 >= $15) { if ($14 >= $16) {print $0, $14} else { print $0, $16}} else { if ($15 >= $16) {print $0, $15} else {print $0, $16}}}' OFS='\t' > "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.length.overlapSizes.bed"

#find individual overlaps (where individual genomeMaps are completely overlapped) in tripleMappers
cat "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.length.overlapSizes.bed" |   awk '{ if ($14 >= $16) { if (($5 * 0.85) <= $14) { print $0 }} else {if (($5 * 0.85) <= $16) { print $0 }}}' OFS='\t'     >  "./finalTry/tripleMap.readBed.readBed.hs1.Overlap.bed"
cat "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.length.overlapSizes.bed" |   awk '{ if ($14 >= $15) { if (($9 * 0.85) <= $14) { print $0 }} else {if (($9 * 0.85) <= $15) { print $0 }}}' OFS='\t'     >  "./finalTry/tripleMap.readBed.readBed.mm39.Overlap.bed"
cat "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.length.overlapSizes.bed" |   awk '{ if ($15 >= $16) { if (($13 * 0.85) <= $15) { print $0 }} else {if (($13 * 0.85) <= $16) { print $0 }}}' OFS='\t'   >  "./finalTry/tripleMap.readBed.readBed.rheMac10.Overlap.bed"

#tripleMappers reads without an overlap are artefacts
cat "./finalTry/tripleMap.readBed.triple.hs1.mm39.rheMac10.length.overlapSizes.bed" |  fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.readBed.hs1.Overlap.bed" ) |  fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.readBed.mm39.Overlap.bed" ) |  fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.readBed.rheMac10.Overlap.bed" )  > "./finalTry/triple.ARTEFACT.tripleGenome.readbed.hs1.mm39.rheMac10.bed"

#find tripleMappers reads with all shared overlaps
cat "./finalTry/tripleMap.readBed.readBed.hs1.Overlap.bed"  | fgrep -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.readBed.mm39.Overlap.bed" )  | fgrep -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.readBed.rheMac10.Overlap.bed" )   > "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed"

##find tripleMappers reads with shared overlaps; these are artefacts with the third genome ##withLengthFilter
#cat "./finalTry/tripleMap.readBed.readBed.hs1.Overlap.bed"  | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") | fgrep -f <( awk '{print $1}' OFS='\t' "./finalTry/overlap.readBed.readBed.mm39.Overlap.bed" )     | awk -v size="${minMapSize}" '{if ($17 < (0.85*(size))) { print $0} else if ($15 == $16) {print $0}}' OFS='\t'  > "./finalTry/triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed"
#cat "./finalTry/tripleMap.readBed.readBed.hs1.Overlap.bed"  | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") | fgrep -f <( awk '{print $1}' OFS='\t' "./finalTry/overlap.readBed.readBed.rheMac10.Overlap.bed" ) | awk -v size="${minMapSize}" '{if ($17 < (0.85*(size))) { print $0} else if ($14 == $15) {print $0}}' OFS='\t'  > "./finalTry/triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed"
#cat "./finalTry/tripleMap.readBed.readBed.mm39.Overlap.bed" | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") | fgrep -f <( awk '{print $1}' OFS='\t' "./finalTry/overlap.readBed.readBed.rheMac10.Overlap.bed" ) | awk -v size="${minMapSize}" '{if ($17 <( 0.85*(size))) { print $0} else if ($14 == $16) {print $0}}' OFS='\t'  > "./finalTry/triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed"

#find tripleMappers reads with shared overlaps; these are artefacts with the third genome or bigMaps ##withOUTLengthFilter
cat "./finalTry/tripleMap.readBed.readBed.hs1.Overlap.bed"  | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") | fgrep -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.readBed.mm39.Overlap.bed" )     > "./finalTry/temp.triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed"
cat "./finalTry/tripleMap.readBed.readBed.hs1.Overlap.bed"  | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") | fgrep -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.readBed.rheMac10.Overlap.bed" ) > "./finalTry/temp.triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed"
cat "./finalTry/tripleMap.readBed.readBed.mm39.Overlap.bed" | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") | fgrep -f <( awk '{print $1}' OFS='\t' "./finalTry/tripleMap.readBed.readBed.rheMac10.Overlap.bed" ) > "./finalTry/temp.triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed"

#remove reads with double and triple overlaps from single overlap reads
cat "./finalTry/tripleMap.readBed.readBed.hs1.Overlap.bed"      | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/temp.triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed" ) | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/temp.triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed") > "./finalTry/triple.ARTEFACT.mouse.to.rheMac.bed"
cat "./finalTry/tripleMap.readBed.readBed.mm39.Overlap.bed"     | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/temp.triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed" ) | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/temp.triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed") > "./finalTry/triple.ARTEFACT.human.to.rheMac.bed"
cat "./finalTry/tripleMap.readBed.readBed.rheMac10.Overlap.bed" | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/temp.triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed")  | fgrep -v -f <( awk '{print $1}' OFS='\t' "./finalTry/temp.triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed") > "./finalTry/triple.ARTEFACT.human.to.mouse.bed"

#Find those reads where one bigMap also cover the alignments to other genomes, but is not considered an overlap
cat "./finalTry/temp.triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed" |   awk '{ if (($2 * 0.85) < $6) { if (($2 * 0.85) < $10) { if ($3 > ($7 * 0.85)) { if ($3 > ($11 * 0.85)) { print $0 }}}}}'   OFS='\t'  >  "./finalTry/tripleMap.bigMap.hs1.readBed.bed"
cat "./finalTry/temp.triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed" |   awk '{ if (($6 * 0.85) < $2) { if (($6 * 0.85) < $10) { if ($7 > ($3 * 0.85)) { if ($7 > ($11 * 0.85)) { print $0 }}}}}'   OFS='\t'  >  "./finalTry/tripleMap.bigMap.mm39.readBed.bed"
cat "./finalTry/temp.triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed" |   awk '{ if (($10 * 0.85) < $6) { if (($10 * 0.85) < $2) { if ($11 > ($7 * 0.85)) { if ($11 > ($3 * 0.85)) { print $0 }}}}}' OFS='\t'  >  "./finalTry/tripleMap.bigMap.rheMac10.readBed.bed"

#remove bigMaps from doubleOverlap tripleMappers
cat "./finalTry/temp.triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed" | fgrep -v -f  "./finalTry/tripleMap.bigMap.hs1.readBed.bed"      > "./finalTry/triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed"
cat "./finalTry/temp.triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed" | fgrep -v -f  "./finalTry/tripleMap.bigMap.mm39.readBed.bed"     > "./finalTry/triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed"
cat "./finalTry/temp.triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed" | fgrep -v -f  "./finalTry/tripleMap.bigMap.rheMac10.readBed.bed" > "./finalTry/triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed"

#figure out what genome the triple overlaps originate
join -t $'\t' <(sort <( awk '{print $1}' OFS='\t' "./finalTry/triple.overlap.sharedRegionBy.human.mouse.rheMac.Overlap.bed") ) <(sort    <(join -t $'\t' <(sort "./finalTry/hs1.filtered.${minMapSize}bp.readNames.Scores.tsv") <(sort  <(join -t $'\t' <(sort "./finalTry/mm39.filtered.${minMapSize}bp.readNames.Scores.tsv") <(sort    "./finalTry/rheMac10.filtered.${minMapSize}bp.readNames.Scores.tsv") )))) > "./finalTry/triple.unknownOrigin.tripleOverlap.Scores.tsv"

cat "./finalTry/triple.unknownOrigin.tripleOverlap.Scores.tsv" | awk '{if ($2 > $3) { if ($2 > $4) {print $1}}}' OFS='\t' > "./finalTry/triple.unknownOrigin.tripleOverlap.HUMAN.bed"
cat "./finalTry/triple.unknownOrigin.tripleOverlap.Scores.tsv" | awk '{if ($4 > $3) { if ($4 > $2) {print $1}}}' OFS='\t' > "./finalTry/triple.unknownOrigin.tripleOverlap.RHEMAC.bed"
cat "./finalTry/triple.unknownOrigin.tripleOverlap.Scores.tsv" | awk '{if ($3 > $2) { if ($3 > $4) {print $1}}}' OFS='\t' > "./finalTry/triple.unknownOrigin.tripleOverlap.MOUSE.bed"
cat "./finalTry/triple.unknownOrigin.tripleOverlap.Scores.tsv" | fgrep -v -f "./finalTry/triple.unknownOrigin.tripleOverlap.HUMAN.bed" | fgrep -v -f  "./finalTry/triple.unknownOrigin.tripleOverlap.MOUSE.bed" | fgrep -v -f "./finalTry/triple.unknownOrigin.tripleOverlap.RHEMAC.bed" > "./finalTry/triple.unknownOrigin.tripleOverlap.AMBIG.bed"

#deal with double overlaps: these are artefacts between the third genome and the best mapper over the shared region
join -t $'\t' <(sort <( awk '{print $1}' OFS='\t' "./finalTry/triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed" ) ) <(sort  <( join -t $'\t' <(sort  "./finalTry/hs1.filtered.${minMapSize}bp.readNames.Scores.tsv")   <(sort    "./finalTry/mm39.filtered.${minMapSize}bp.readNames.Scores.tsv")))     >  "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.bed"
join -t $'\t' <(sort <( awk '{print $1}' OFS='\t' "./finalTry/triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed" ) ) <(sort  <( join -t $'\t' <(sort  "./finalTry/mm39.filtered.${minMapSize}bp.readNames.Scores.tsv")  <(sort    "./finalTry/rheMac10.filtered.${minMapSize}bp.readNames.Scores.tsv"))) >  "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.bed"
join -t $'\t' <(sort <( awk '{print $1}' OFS='\t' "./finalTry/triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed" ) ) <(sort  <( join -t $'\t' <(sort  "./finalTry/hs1.filtered.${minMapSize}bp.readNames.Scores.tsv")   <(sort    "./finalTry/rheMac10.filtered.${minMapSize}bp.readNames.Scores.tsv"))) >  "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.bed"

cat "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.bed" | awk '{if ($2 > $3) {print $0}}' OFS='\t'  >  "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.HUMAN.tsv"
cat "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.bed" | awk '{if ($2 < $3) {print $0}}' OFS='\t'  >  "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.MOUSE.tsv"
cat "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.bed" | awk '{if ($2 == $3) {print $0}}' OFS='\t' >  "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.AMBIG.tsv"

cat "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.bed"  | awk '{if ($2 > $3) {print $0}}' OFS='\t'  >  "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.HUMAN.tsv"
cat "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.bed"  | awk '{if ($2 < $3) {print $0}}' OFS='\t'  >  "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.RHEMAC.tsv"
cat "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.bed"  | awk '{if ($2 == $3) {print $1}}' OFS='\t' >  "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.AMBIG.tsv"

cat "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.bed"  | awk '{if ($2 > $3) {print $0}}' OFS='\t'  >  "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.MOUSE.tsv"
cat "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.bed"  | awk '{if ($2 < $3) {print $0}}' OFS='\t'  >  "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.RHEMAC.tsv"
cat "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.bed"  | awk '{if ($2 == $3) {print $0}}' OFS='\t' >  "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.AMBIG.tsv"


cat "./finalTry/triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed" | fgrep -f <(awk '{print $1}' OFS='\t' "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.HUMAN.tsv") > "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.HUMAN.bed"
cat "./finalTry/triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed" | fgrep -f <(awk '{print $1}' OFS='\t' "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.MOUSE.tsv") > "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.MOUSE.bed"
cat "./finalTry/triple.rheMac.ARTEFACT.overlap.sharedRegionBy.human.mouse.Overlap.bed" | fgrep -f <(awk '{print $1}' OFS='\t' "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.AMBIG.tsv") > "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.AMBIG.bed"

cat "./finalTry/triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed" | fgrep -f <(awk '{print $1}' OFS='\t' "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.MOUSE.tsv") > "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.MOUSE.bed"
cat "./finalTry/triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed" | fgrep -f <(awk '{print $1}' OFS='\t' "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.RHEMAC.tsv") > "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.RHEMAC.bed"
cat "./finalTry/triple.human.ARTEFACT.overlap.sharedRegionBy.mouse.rheMac.Overlap.bed" | fgrep -f <(awk '{print $1}' OFS='\t' "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.AMBIG.tsv") > "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.AMBIG.bed"

cat "./finalTry/triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed" | fgrep -f <(awk '{print $1}' OFS='\t' "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.HUMAN.tsv") > "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.HUMAN.bed"
cat "./finalTry/triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed" | fgrep -f <(awk '{print $1}' OFS='\t' "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.RHEMAC.tsv") > "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.RHEMAC.bed"
cat "./finalTry/triple.mouse.ARTEFACT.overlap.sharedRegionBy.human.rheMac.Overlap.bed" | fgrep -f <(awk '{print $1}' OFS='\t' "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.AMBIG.tsv") > "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.AMBIG.bed"



#remapping readParts to reduce false Mapping results 
#Start by selecting the reads to be remapped
#starting with the base for the soloGenome Reads not mapped to multiple genomes
cat "./finalTry/hs1.readBed.bed"      | fgrep -v -f "./finalTry/badReads.txt" |  fgrep -v -f <(awk '{ print $1}' OFS='\t'  "./finalTry/temp.doubleMap.readBed.hs1.mm39.bed" )     | fgrep -v -f <(awk '{ print $1}' OFS='\t'   "./finalTry/temp.doubleMap.readBed.hs1.rheMac10.bed" )  | awk '{ print $1}' OFS='\t' > "./finalTry/hs1.Reads.names.txt"
cat "./finalTry/mm39.readBed.bed"     | fgrep -v -f "./finalTry/badReads.txt" |  fgrep -v -f <(awk '{ print $1}' OFS='\t'  "./finalTry/temp.doubleMap.readBed.hs1.mm39.bed" )     | fgrep -v -f <(awk '{ print $1}' OFS='\t'   "./finalTry/temp.doubleMap.readBed.mm39.rheMac10.bed" ) | awk '{ print $1}' OFS='\t' > "./finalTry/mm39.Reads.names.txt"
cat "./finalTry/rheMac10.readBed.bed" | fgrep -v -f "./finalTry/badReads.txt" |  fgrep -v -f <(awk '{ print $1}' OFS='\t'  "./finalTry/temp.doubleMap.readBed.hs1.rheMac10.bed" ) | fgrep -v -f <(awk '{ print $1}' OFS='\t'  "./finalTry/temp.doubleMap.readBed.mm39.rheMac10.bed" )  | awk '{ print $1}' OFS='\t' > "./finalTry/rheMac10.Reads.names.txt"

#soloGenome Reads mapped to two genomes
cat "./finalTry/unknownOrigin.hs1.mm39.readName.HUMAN.txt" >> "./finalTry/hs1.Reads.names.txt"
cat "./finalTry/unknownOrigin.hs1.mm39.readName.MOUSE.txt" >> "./finalTry/mm39.Reads.names.txt"

cat "./finalTry/unknownOrigin.hs1.rheMac10.readName.HUMAN.txt"  >> "./finalTry/hs1.Reads.names.txt"
cat "./finalTry/unknownOrigin.hs1.rheMac10.readName.RHEMAC.txt" >> "./finalTry/rheMac10.Reads.names.txt"

cat "./finalTry/unknownOrigin.mm39.rheMac10.readName.MOUSE.txt"  >> "./finalTry/mm39.Reads.names.txt"
cat "./finalTry/unknownOrigin.mm39.rheMac10.readName.RHEMAC.txt" >> "./finalTry/rheMac10.Reads.names.txt"

#soloGenome Reads from reads mapped three genomes with triple overlaps
cat "./finalTry/triple.unknownOrigin.tripleOverlap.HUMAN.bed" | awk '{print $1}' OFS='\t'  >> "./finalTry/hs1.Reads.names.txt"
cat "./finalTry/triple.unknownOrigin.tripleOverlap.MOUSE.bed" | awk '{print $1}' OFS='\t'  >> "./finalTry/mm39.Reads.names.txt"
cat "./finalTry/triple.unknownOrigin.tripleOverlap.RHEMAC.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/rheMac10.Reads.names.txt"

#sologenome Reads from reads mapped to three genomes, with double overlaps, but all is covered by a bigMap
cat "./finalTry/tripleMap.bigMap.hs1.readBed.bed"      >> "./finalTry/hs1.Reads.names.txt"
cat "./finalTry/tripleMap.bigMap.mm39.readBed.bed"     >> "./finalTry/mm39.Reads.names.txt"
cat "./finalTry/tripleMap.bigMap.rheMac10.readBed.bed" >> "./finalTry/rheMac10.Reads.names.txt"

#NOW ARTEFACTS
#doubleGenome ARTEFACTS from reads mapped to two genomes
cat "./finalTry/ARTEFACT.doubleMap.readBed.hs1.mm39.length.overlapSize.bed"      > "./finalTry/doubleMap.ARTEFACT.hs1.mm39.readBed.bed"
cat "./finalTry/ARTEFACT.doubleMap.readBed.hs1.rheMac10.length.overlapSize.bed"  > "./finalTry/doubleMap.ARTEFACT.hs1.rheMac10.readBed.bed"
cat "./finalTry/ARTEFACT.doubleMap.readBed.mm39.rheMac10.length.overlapSize.bed" > "./finalTry/doubleMap.ARTEFACT.mm39.rheMac10.readBed.bed"


cat "./finalTry/doubleMap.ARTEFACT.hs1.mm39.readBed.bed"      | awk '{print $1}' OFS='\t' > "./finalTry/hs1.ARTEFACT.readNames.txt"
cat "./finalTry/doubleMap.ARTEFACT.hs1.mm39.readBed.bed"      | awk '{print $1}' OFS='\t' > "./finalTry/mm39.ARTEFACT.readNames.txt"

cat "./finalTry/doubleMap.ARTEFACT.hs1.rheMac10.readBed.bed"  | awk '{print $1}' OFS='\t' >> "./finalTry/hs1.ARTEFACT.readNames.txt"
cat "./finalTry/doubleMap.ARTEFACT.hs1.rheMac10.readBed.bed"  | awk '{print $1}' OFS='\t' > "./finalTry/rheMac10.ARTEFACT.readNames.txt"

cat "./finalTry/doubleMap.ARTEFACT.mm39.rheMac10.readBed.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/mm39.ARTEFACT.readNames.txt"
cat "./finalTry/doubleMap.ARTEFACT.mm39.rheMac10.readBed.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/rheMac10.ARTEFACT.readNames.txt"


cat "./finalTry/doubleMap.ARTEFACT.hs1.mm39.readBed.bed"      | awk '{print $1}' OFS='\t' > "./finalTry/ARTEFACT.readNames.txt"
cat "./finalTry/doubleMap.ARTEFACT.hs1.rheMac10.readBed.bed"  | awk '{print $1}' OFS='\t' >> "./finalTry/ARTEFACT.readNames.txt"
cat "./finalTry/doubleMap.ARTEFACT.mm39.rheMac10.readBed.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/ARTEFACT.readNames.txt"

#doubleGenome ARTEFACTS from reads mapped to three genomes with 2 overlap
cat "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.HUMAN.bed" > "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.rheMac10.bed"
cat "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.RHEMAC.bed" >> "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.rheMac10.bed"

cat "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.HUMAN.bed" > "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.mm39.bed"
cat "./finalTry/triple.ARTEFACT.human.overlap.sharedRegionBy.FoundOverlap.MOUSE.bed" >> "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.mm39.bed"

cat "./finalTry/triple.ARTEFACT.mouse.overlap.sharedRegionBy.FoundOverlap.RHEMAC.bed" > "./finalTry/tripleMap.ARTEFACT.doubleOverlap.mm39.rheMac10.bed"
cat "./finalTry/triple.ARTEFACT.rheMac.overlap.sharedRegionBy.FoundOverlap.MOUSE.bed" >> "./finalTry/tripleMap.ARTEFACT.doubleOverlap.mm39.rheMac10.bed"


cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.rheMac10.bed"  | awk '{print $1}' OFS='\t' >> "./finalTry/hs1.ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.rheMac10.bed"  | awk '{print $1}' OFS='\t' >> "./finalTry/rheMac10.ARTEFACT.readNames.txt"

cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.mm39.bed"      | awk '{print $1}' OFS='\t' >> "./finalTry/hs1.ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.mm39.bed"      | awk '{print $1}' OFS='\t' >> "./finalTry/mm39.ARTEFACT.readNames.txt"

cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.mm39.rheMac10.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/mm39.ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.mm39.rheMac10.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/rheMac10.ARTEFACT.readNames.txt"


cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.rheMac10.bed"  | awk '{print $1}' OFS='\t' >> "./finalTry/ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.mm39.bed"      | awk '{print $1}' OFS='\t' >> "./finalTry/ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.mm39.rheMac10.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/ARTEFACT.readNames.txt"

#doubleGenome ARTEFACTS from reads mapped to three genomes with 1 overlap
cat "./finalTry/triple.ARTEFACT.human.to.mouse.bed"  > "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.mm39.bed"
cat "./finalTry/triple.ARTEFACT.human.to.rheMac.bed" > "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.rheMac10.bed"
cat "./finalTry/triple.ARTEFACT.mouse.to.rheMac.bed" > "./finalTry/tripleMap.ARTEFACT.singleOverlap.mm39.rheMac10.bed"


cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.mm39.bed"      | awk '{print $1}' OFS='\t' >> "./finalTry/hs1.ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.mm39.bed"      | awk '{print $1}' OFS='\t' >> "./finalTry/mm39.ARTEFACT.readNames.txt"

cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.rheMac10.bed"  | awk '{print $1}' OFS='\t' >> "./finalTry/hs1.ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.rheMac10.bed"  | awk '{print $1}' OFS='\t' >> "./finalTry/rheMac10.ARTEFACT.readNames.txt"

cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.mm39.rheMac10.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/mm39.ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.mm39.rheMac10.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/rheMac10.ARTEFACT.readNames.txt"


cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.mm39.bed"      | awk '{print $1}' OFS='\t' >> "./finalTry/ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.rheMac10.bed"  | awk '{print $1}' OFS='\t' >> "./finalTry/ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.mm39.rheMac10.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/ARTEFACT.readNames.txt"


##tripleGenome ARTEFACTS from reads mapped to three genomes with 0 overlaps
cat "./finalTry/triple.ARTEFACT.tripleGenome.readbed.hs1.mm39.rheMac10.bed" > "./finalTry/tripleMap.ARTEFACT.tripleGenome.hs1.mm39.rheMac10.readBed.bed"


cat "./finalTry/tripleMap.ARTEFACT.tripleGenome.hs1.mm39.rheMac10.readBed.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/hs1.ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.tripleGenome.hs1.mm39.rheMac10.readBed.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/mm39.ARTEFACT.readNames.txt"
cat "./finalTry/tripleMap.ARTEFACT.tripleGenome.hs1.mm39.rheMac10.readBed.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/rheMac10.ARTEFACT.readNames.txt"


cat "./finalTry/tripleMap.ARTEFACT.tripleGenome.hs1.mm39.rheMac10.readBed.bed" | awk '{print $1}' OFS='\t' >> "./finalTry/ARTEFACT.readNames.txt"

rm "./finalTry/badRegions.badReads.txt"
mkdir ./finalTry/readCounts

##for every reference (==sample): prepare the readBed of the primary alignment
for sample in "${samples[@]}"; do
    
    echo "remapping with: ${sample} reads"
    conda activate seqkitEnv
    #get sequences of indivudual mapping parts of putative artefact reads
    seqkit subseq --bed <(fgrep -f "./finalTry/${sample}.ARTEFACT.readNames.txt" "./finalTry/${sample}.readBed.bed")  ${fastqFile} > "./finalTry/putative.artefact.${sample}Part.fastq"
    #get parts of reads aligned to individual genomes
    seqkit subseq --bed <(fgrep -f "./finalTry/${sample}.Reads.names.txt" "./finalTry/${sample}.readBed.bed")  ${fastqFile} > "./finalTry/putative.${sample}Only.${sample}Part.fastq"
    conda deactivate

    conda activate minimapEnv
    #remapping readpartsto remove 
    for sampleSec in "${samples[@]}"; do
        #remap individual parts of putative artefact reads
        minimap2 -a -x map-ont -t 16  "/home/tnsmits/lustre/lib/${sampleSec}.fa" "./finalTry/putative.artefact.${sample}Part.fastq" --secondary=no > "${path}/finalTry/putative.artefact.${sample}Part.${sampleSec}Mapped.sam"
        #remap parts of reads aligned to individual genomes
        minimap2 -a -x map-ont -t 16  "/home/tnsmits/lustre/lib/${sampleSec}.fa" "./finalTry/putative.${sample}Only.${sample}Part.fastq" --secondary=no > "${path}/finalTry/putative.${sample}Only.${sampleSec}Mapped.sam"
    done
    conda deactivate
   
    #keep only good maps
    conda activate samtoolsEnv
    for sampleSec in "${samples[@]}"; do
        samtools view -S -b -F 256  -F 0x800 -q 60  "${path}/finalTry/putative.artefact.${sample}Part.${sampleSec}Mapped.sam" > "${path}/finalTry/putative.artefact.${sample}Part.${sampleSec}Mapped.bam"
        samtools view -S -b -F 256  -F 0x800 -q 60  "${path}/finalTry/putative.${sample}Only.${sampleSec}Mapped.sam"          > "${path}/finalTry/putative.${sample}Only.${sampleSec}Mapped.bam"
        
        samtools view  "${path}/finalTry/putative.artefact.${sample}Part.${sampleSec}Mapped.bam" | awk '{ print $1, $14}' OFS='\t' | sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\(AS:i:\)\([0-9]*\)/\1\2\4/' | sed 's/_[0-9]*\-[0-9]*:.//'   | sort > "${path}/finalTry/putative.artefact.${sample}Part.${sampleSec}Mapped.Scores.tsv"
        samtools view  "${path}/finalTry/putative.${sample}Only.${sampleSec}Mapped.bam"          | awk '{ print $1, $14}' OFS='\t' | sed 's/^\([^[:blank:]]*\)\([[:blank:]]*\)\(AS:i:\)\([0-9]*\)/\1\2\4/' | sed 's/_[0-9]*\-[0-9]*:.//'   | sort > "${path}/finalTry/putative.${sample}Only.${sampleSec}Mapped.Scores.tsv"

    done
    
    samtools sort  "${path}/finalTry/putative.${sample}Only.${sample}Mapped.bam" -o "${path}/finalTry/putative.${sample}Only.${sample}Mapped.sorted.bam"
    samtools index "${path}/finalTry/putative.${sample}Only.${sample}Mapped.sorted.bam"

    conda deactivate
    
    conda activate mosdepthEnv
    
    cd "${path}/finalTry"
    mosdepth -t 16 -b 1000 -n -m ${sample} "${path}/finalTry/putative.${sample}Only.${sample}Mapped.sorted.bam"
    cd ..
    meanCoverage=$(awk '{ if ($1 == "total") {if ($3 != 0.0){ if ( $3 <= 0.05) { print $2} else { print 1}} else { print 1}}}' "./finalTry/${sample}.mosdepth.global.dist.txt" OFS='\t' | sort -h | tail -n 1)
    echo "MeanCoverage ${sample}: ${meanCoverage}"
    
    zcat "./finalTry/${sample}.regions.bed.gz"  | awk -v COV="${meanCoverage}" '{ if ( $4 > COV) { print $0}}'  OFS='\t' > "./finalTry/badRegions.${sample}.bed"
    conda deactivate
    
    conda activate bedtoolsEnv
    bedtools bamtobed -i "${path}/finalTry/putative.${sample}Only.${sample}Mapped.bam" > "${path}/finalTry/putative.${sample}Only.${sample}Mapped.bed"
    bedtools intersect -wa -a "${path}/finalTry/putative.${sample}Only.${sample}Mapped.bed" -b  "./finalTry/badRegions.${sample}.bed" | awk '{print $4}' OFS='\t' >> "./finalTry/badRegions.badReads.txt"
    conda deactivate
    
    
    join -a 1 -a 2 -e'0' -o '0,1.2,2.2' "./finalTry/putative.artefact.${sample}Part.hs1Mapped.Scores.tsv"  "./finalTry/putative.artefact.${sample}Part.mm39Mapped.Scores.tsv" | join -a 1 -a 2 -e'0' -o '0,1.2,1.3,2.2' -  "./finalTry/putative.artefact.${sample}Part.rheMac10Mapped.Scores.tsv"  > "./finalTry/putative.artefact.${sample}Part.allScores.tsv"
    
    
    cat "./finalTry/putative.artefact.${sample}Part.allScores.tsv" | awk '{if ($2 > $3) { if ($2 > $4) {print $1}}}' OFS='\t' | fgrep -v -f "./finalTry/badRegions.badReads.txt" > "./finalTry/putative.artefact.${sample}Part.allScores.HUMAN.tsv"
    cat "./finalTry/putative.artefact.${sample}Part.allScores.tsv" | awk '{if ($3 > $4) { if ($3 > $2) {print $1}}}' OFS='\t' | fgrep -v -f "./finalTry/badRegions.badReads.txt" > "./finalTry/putative.artefact.${sample}Part.allScores.MOUSE.tsv"
    cat "./finalTry/putative.artefact.${sample}Part.allScores.tsv" | awk '{if ($4 > $2) { if ($4 > $3) {print $1}}}' OFS='\t' | fgrep -v -f "./finalTry/badRegions.badReads.txt" > "./finalTry/putative.artefact.${sample}Part.allScores.RHEMAC.tsv"
    cat "./finalTry/putative.artefact.${sample}Part.allScores.tsv" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f "./finalTry/putative.artefact.${sample}Part.allScores.HUMAN.tsv" | fgrep -v -f  "./finalTry/putative.artefact.${sample}Part.allScores.MOUSE.tsv" | fgrep -v -f "./finalTry/putative.artefact.${sample}Part.allScores.RHEMAC.tsv" > "./finalTry/putative.artefact.${sample}Part.allScores.AMBIG.tsv"

    join -a 1 -a 2 -e'0' -o '0,1.2,2.2' "./finalTry/putative.${sample}Only.hs1Mapped.Scores.tsv"  "./finalTry/putative.${sample}Only.mm39Mapped.Scores.tsv" | join -a 1 -a 2 -e'0' -o '0,1.2,1.3,2.2' -  "./finalTry/putative.${sample}Only.rheMac10Mapped.Scores.tsv"  > "./finalTry/putative.${sample}Only.allScores.tsv"

    cat "./finalTry/putative.${sample}Only.allScores.tsv" | awk '{if ($2 > $3) { if ($2 > $4) {print $1}}}' OFS='\t' | fgrep -v -f "./finalTry/badRegions.badReads.txt" > "./finalTry/putative.${sample}Only.allScores.HUMAN.tsv"
    cat "./finalTry/putative.${sample}Only.allScores.tsv" | awk '{if ($3 > $4) { if ($3 > $2) {print $1}}}' OFS='\t' | fgrep -v -f "./finalTry/badRegions.badReads.txt" > "./finalTry/putative.${sample}Only.allScores.MOUSE.tsv"
    cat "./finalTry/putative.${sample}Only.allScores.tsv" | awk '{if ($4 > $2) { if ($4 > $3) {print $1}}}' OFS='\t' | fgrep -v -f "./finalTry/badRegions.badReads.txt" > "./finalTry/putative.${sample}Only.allScores.RHEMAC.tsv"
    cat "./finalTry/putative.${sample}Only.allScores.tsv" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f "./finalTry/putative.${sample}Only.allScores.HUMAN.tsv" | fgrep -v -f  "./finalTry/putative.${sample}Only.allScores.MOUSE.tsv" | fgrep -v -f "./finalTry/putative.${sample}Only.allScores.RHEMAC.tsv" > "./finalTry/putative.${sample}Only.allScores.AMBIG.tsv"

    
done

for sample in "${samples[@]}"; do
    for s2 in "${samples[@]}"; do
    
        if [[ "${sample}" == "hs1" ]]; then
            sampleOne="HUMAN"
        elif [[ "${sample}" == "mm39" ]]; then
            sampleOne="MOUSE"
        else
            sampleOne="RHEMAC"
        fi

    
        echo "Archiving final ${sample} reads"
        #archiving nonArtefact readCounts
        cat "./finalTry/putative.${sample}Only.allScores.${sampleOne}.tsv" | awk '{ print $1}' OFS='\t'  | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt"  > "./finalTry/readCounts/${sample}.Reads.names.txt"
        if [[ "${sample}" != "${s2}" ]]; then
            #Adding nonArtefact reads that were originally wrongly assigned
            cat "./finalTry/putative.${s2}Only.allScores.${sampleOne}.tsv" | awk '{ print $1}' OFS='\t'  | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt"  >>  "./finalTry/readCounts/${sample}.Reads.names.txt"

        fi

    
    done
done
#NOW ARTEFACTS
cat "./finalTry/hs1.ARTEFACT.readNames.txt"        | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -f "./finalTry/putative.artefact.hs1Part.allScores.HUMAN.tsv"       > "./finalTry/readCounts/hs1.ARTEFACT.readNames.txt"
cat "./finalTry/mm39.ARTEFACT.readNames.txt"       | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -f "./finalTry/putative.artefact.mm39Part.allScores.MOUSE.tsv"      > "./finalTry/readCounts/mm39.ARTEFACT.readNames.txt"
cat "./finalTry/rheMac10.ARTEFACT.readNames.txt"   | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -f "./finalTry/putative.artefact.rheMac10Part.allScores.RHEMAC.tsv" > "./finalTry/readCounts/rheMac10.ARTEFACT.readNames.txt"
cat "./finalTry/ARTEFACT.readNames.txt"            | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" > "./finalTry/readCounts/ARTEFACT.readNames.txt" 
cat "./finalTry/readCounts/ARTEFACT.readNames.txt" | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f "./finalTry/putative.artefact.hs1Part.allScores.HUMAN.tsv" | fgrep -v -f "./finalTry/putative.artefact.mm39Part.allScores.MOUSE.tsv" | fgrep -v -f "./finalTry/putative.artefact.rheMac10Part.allScores.RHEMAC.tsv" > "./finalTry/readCounts/other.ARTEFACT.readNames.txt"

#doubleGenome ARTEFACTS from reads mapped to two genomes
cat "./finalTry/doubleMap.ARTEFACT.hs1.mm39.readBed.bed"      | fgrep -f "./finalTry/putative.artefact.hs1Part.allScores.HUMAN.tsv" | fgrep -f "./finalTry/putative.artefact.mm39Part.allScores.MOUSE.tsv"                 > "./finalTry/readCounts/doubleMap.ARTEFACT.hs1.mm39.readBed.bed"
cat "./finalTry/doubleMap.ARTEFACT.hs1.mm39.readBed.bed"      | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t' "./finalTry/readCounts/doubleMap.ARTEFACT.hs1.mm39.readBed.bed")       > "./finalTry/readCounts/doubleMap.ARTEFACT.other.readBed.bed"


cat "./finalTry/doubleMap.ARTEFACT.hs1.rheMac10.readBed.bed"  | fgrep -f "./finalTry/putative.artefact.hs1Part.allScores.HUMAN.tsv" | fgrep -f "./finalTry/putative.artefact.rheMac10Part.allScores.RHEMAC.tsv"            > "./finalTry/readCounts/doubleMap.ARTEFACT.hs1.rheMac10.readBed.bed"
cat "./finalTry/doubleMap.ARTEFACT.hs1.rheMac10.readBed.bed"  | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t' "./finalTry/readCounts/doubleMap.ARTEFACT.hs1.rheMac10.readBed.bed")   >> "./finalTry/readCounts/doubleMap.ARTEFACT.other.readBed.bed"


cat "./finalTry/doubleMap.ARTEFACT.mm39.rheMac10.readBed.bed" | fgrep -f "./finalTry/putative.artefact.mm39Part.allScores.MOUSE.tsv" | fgrep -f "./finalTry/putative.artefact.rheMac10Part.allScores.RHEMAC.tsv"             > "./finalTry/readCounts/doubleMap.ARTEFACT.mm39.rheMac10.readBed.bed"
cat "./finalTry/doubleMap.ARTEFACT.mm39.rheMac10.readBed.bed" | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t' "./finalTry/readCounts/doubleMap.ARTEFACT.mm39.rheMac10.readBed.bed")   >> "./finalTry/readCounts/doubleMap.ARTEFACT.other.readBed.bed"


#doubleGenome ARTEFACTS from reads mapped to three genomes with 2 overlap
cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.rheMac10.bed" | fgrep -f "./finalTry/putative.artefact.hs1Part.allScores.HUMAN.tsv" | fgrep -f "./finalTry/putative.artefact.rheMac10Part.allScores.RHEMAC.tsv"  > "./finalTry/readCounts/tripleMap.ARTEFACT.doubleOverlap.hs1.rheMac10.bed" 
cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.rheMac10.bed" | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t'  "./finalTry/readCounts/tripleMap.ARTEFACT.doubleOverlap.hs1.rheMac10.bed") >  "./finalTry/readCounts/tripleMap.ARTEFACT.doubleOverlap.other.bed"

cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.mm39.bed" | fgrep -f "./finalTry/putative.artefact.hs1Part.allScores.HUMAN.tsv" | fgrep -f "./finalTry/putative.artefact.mm39Part.allScores.MOUSE.tsv"  > "./finalTry/readCounts/tripleMap.ARTEFACT.doubleOverlap.hs1.mm39.bed" 
cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.hs1.mm39.bed" | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t'  "./finalTry/readCounts/tripleMap.ARTEFACT.doubleOverlap.hs1.mm39.bed") >  "./finalTry/readCounts/tripleMap.ARTEFACT.doubleOverlap.other.bed"

cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.mm39.rheMac10.bed" | fgrep -f "./finalTry/putative.artefact.mm39Part.allScores.MOUSE.tsv" | fgrep -f "./finalTry/putative.artefact.rheMac10Part.allScores.RHEMAC.tsv"  > "./finalTry/readCounts/tripleMap.ARTEFACT.doubleOverlap.mm39.rheMac10.bed" 
cat "./finalTry/tripleMap.ARTEFACT.doubleOverlap.mm39.rheMac10.bed" | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t'  "./finalTry/readCounts/tripleMap.ARTEFACT.doubleOverlap.mm39.rheMac10.bed") >  "./finalTry/readCounts/tripleMap.ARTEFACT.doubleOverlap.other.bed"

#doubleGenome ARTEFACTS from reads mapped to three genomes with 1 overlap
cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.mm39.bed"  | fgrep -f "./finalTry/putative.artefact.hs1Part.allScores.HUMAN.tsv" | fgrep -f "./finalTry/putative.artefact.mm39Part.allScores.MOUSE.tsv" > "./finalTry/readCounts/tripleMap.ARTEFACT.singleOverlap.hs1.mm39.bed"
cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.mm39.bed"  |  fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t' "./finalTry/readCounts/tripleMap.ARTEFACT.singleOverlap.hs1.mm39.bed" ) >"./finalTry/readCounts/tripleMap.ARTEFACT.singleOverlap.other.bed"

cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.rheMac10.bed"  | fgrep -f "./finalTry/putative.artefact.hs1Part.allScores.HUMAN.tsv" | fgrep -f "./finalTry/putative.artefact.rheMac10Part.allScores.RHEMAC.tsv" > "./finalTry/readCounts/tripleMap.ARTEFACT.singleOverlap.hs1.rheMac10.bed"
cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.hs1.rheMac10.bed"  |  fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t' "./finalTry/readCounts/tripleMap.ARTEFACT.singleOverlap.hs1.rheMac10.bed" ) >"./finalTry/readCounts/tripleMap.ARTEFACT.singleOverlap.other.bed"

cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.mm39.rheMac10.bed"  | fgrep -f "./finalTry/putative.artefact.mm39Part.allScores.MOUSE.tsv" | fgrep -f "./finalTry/putative.artefact.rheMac10Part.allScores.RHEMAC.tsv" > "./finalTry/readCounts/tripleMap.ARTEFACT.singleOverlap.mm39.rheMac10.bed"
cat "./finalTry/tripleMap.ARTEFACT.singleOverlap.mm39.rheMac10.bed"  |  fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t' "./finalTry/readCounts/tripleMap.ARTEFACT.singleOverlap.mm39.rheMac10.bed" ) >"./finalTry/readCounts/tripleMap.ARTEFACT.singleOverlap.other.bed"

##tripleGenome ARTEFACTS from reads mapped to three genomes with 0 overlaps
cat "./finalTry/tripleMap.ARTEFACT.tripleGenome.hs1.mm39.rheMac10.readBed.bed" | fgrep -f "./finalTry/putative.artefact.hs1Part.allScores.HUMAN.tsv" | fgrep -f "./finalTry/putative.artefact.mm39Part.allScores.MOUSE.tsv" | fgrep -f "./finalTry/putative.artefact.rheMac10Part.allScores.RHEMAC.tsv" > "./finalTry/readCounts/tripleMap.ARTEFACT.tripleGenome.hs1.mm39.rheMac10.readBed.bed"
cat "./finalTry/tripleMap.ARTEFACT.tripleGenome.hs1.mm39.rheMac10.readBed.bed" | fgrep -v -f "${path}/finalTry/badReads.txt" | fgrep -v -f "./finalTry/badRegions.badReads.txt" | fgrep -v -f <(awk {'print $1'} OFS='\t' "./finalTry/readCounts/tripleMap.ARTEFACT.tripleGenome.hs1.mm39.rheMac10.readBed.bed") > "./finalTry/readCounts/tripleMap.ARTEFACT.tripleGenome.other.readBed.bed"

