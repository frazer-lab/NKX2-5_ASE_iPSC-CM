#!/bin/bash

#$ -N job.MACS2_H3K27AC_CM
#$ -o /home/paola/Family1070/private_output/PeakCalling/H3K27AC_CM/MACS2.out
#$ -e /home/paola/Family1070/private_output/PeakCalling/H3K27AC_CM/MACS2.err
#$ -l week
#$ -pe smp 12
#$ -V

source activate cardips

rootDir=/projects/CARDIPS/pipeline/ChIPseq
control1=5baf99b1-8cb1-11e5-9b81-00259029e99d
control2=1234c342-de08-43b8-aa4a-67b897f4b53d-merged	
readLength=100
mark=H3K27AC_CM

cmd=""
while read name a
do
	cmd=`echo "$cmd $rootDir/sample/$name/alignment/$name.filtered.cordSorted.bam"`
done < /home/paola/Family1070/private_output/PeakCalling/H3K27AC_CM/sample.use.list

cmd2=`echo "macs2 callpeak -t $cmd -c \
$rootDir/sample/$control1/alignment/$control1.filtered.cordSorted.bam $rootDir/sample/$control2/alignment/$control2.filtered.cordSorted.bam \
-f BAMPE -g hs --outdir /home/paola/Family1070/private_output/PeakCalling/H3K27AC_CM/ -n meta_macs2_callPeak -B --SPMR --verbose 3 -s $readLength --call-summits -q 0.01"`

eval $cmd2

err=$(echo $?)
if [ $err != 0 ]
then
	echo "Error occurred when call narrow peaks." >> Error.log
	exit
fi
