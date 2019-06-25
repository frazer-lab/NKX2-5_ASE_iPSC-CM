#!/bin/bash

# column 9 is peak
# column 10 is peak length

# add peak length to column 11

module load cardips
source activate cardips

awk -v OFS="\t" '{$11=$3-$2; print $0}' meta_macs2_callPeak_peaks.narrowPeak > meta_macs2_callPeak_peaks.narrowPeak.temp
mv meta_macs2_callPeak_peaks.narrowPeak.temp meta_macs2_callPeak_peaks.narrowPeak


# filter to keep q<0.01 peaks (non necessary because the algorithm was set to -q 0.01)
awk -v OFS="\t" '$9>2 {print $0}' meta_macs2_callPeak_peaks.narrowPeak > meta_macs2_callPeak_peaks.q001.narrowPeak


# collapse regions with multiple summits

awk -v OFS="\t" '{print $1,$2,$3,$6,$5,$7,$8,$9,$10,$11}' \
meta_macs2_callPeak_peaks.q001.narrowPeak \
| bedtools merge -delim ";" -c 4,5,6,7,8,9,10 \
-o count,mean,mean,mean,mean,collapse,mean -i stdin \
> meta_macs2_callPeak_peaks.q001.narrowPeak.collapse


# make bed file and sorted
awk -v OFS="\t" '{print $1,$2,$3}' \
meta_macs2_callPeak_peaks.q001.narrowPeak.collapse \
> meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.bed

sort -k1,1 -k2,2n meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.bed \
> meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.sorted.bed


# get number of fragments from each sample for each peak

python /frazer01/home/cdeboever/repos/cdeboever3/cdpipelines/cdpipelines/convert_bed_to_saf.py \
meta_macs2_callPeak_peaks.q001.narrowPeak.collapse \
meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.saf

bam_dir=/projects/CARDIPS/pipeline/ChIPseq/mark/NKX
bam_concat=""
while read name alias a
do
bam_concat=`echo "$bam_concat $bam_dir/$name/alignment/$name.filtered.querySorted.bam"`
done < sample.use.list

cmd="featureCounts -p -T 10 -F SAF --donotsort \
-a meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.saf \
-o meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.counts \
$bam_concat"

eval $cmd

sed -i '1d' meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.counts


# change header 
cat sample.use.list | while read uu alias name
do
	sed -i "1s@$bam_dir/$uu/alignment/$uu.filtered.querySorted.bam@$name@" meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.counts
done

# filter sex chromosomes
awk -v OFS="\t" '$2!="chrX" && $2!="chrY" {print $0}' meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.counts > meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.noXY.counts
