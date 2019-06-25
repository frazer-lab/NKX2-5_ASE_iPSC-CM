#!/bin/bash

#$ -N job.MACS2_NKX
#$ -o /home/paola/Family1070/private_output/PeakCalling/NKX25/MACS2.out
#$ -e /home/paola/Family1070/private_output/PeakCalling/NKX25/MACS2.err
#$ -l week
#$ -pe smp 12
#$ -V


source activate cardips

control=32bae3e4-8cb1-11e5-9b81-00259029e99d
readLength=100
plusMinus=100

cmd=""
while read name a b
do
	cmd=`echo "$cmd /frazer01/projects/CARDIPS/pipeline/ChIPseq/mark/NKX/$name/alignment/$name.filtered.cordSorted.bam"`
done < /home/paola/Family1070/private_output/PeakCalling/NKX25/sample.use.list

cmd2=`echo "macs2 callpeak -t $cmd -c \
/frazer01/projects/CARDIPS/pipeline/ChIPseq/sample/$control/alignment/$control.filtered.cordSorted.bam \
-f BAMPE -g hs --outdir /frazer01/home/paola/Family1070/private_output/PeakCalling/NKX25 -n meta_macs2_callPeak -B --SPMR --verbose 3 -s $readLength --call-summits -q 0.01"`

eval $cmd2

err=$(echo $?)
if [ $err != 0 ]
then
	echo "Error occurred when call narrow peaks." >> Error.log
	exit
fi


