#!/bin/bash

#$ -N job.homerMOT
#$ -o /home/paola/Family1070/private_output/Homer_motif_enrichment
#$ -e /home/paola/Family1070/private_output/Homer_motif_enrichment
#$ -l week
#$ -pe smp 12

module load cardips

cd /home/paola/Family1070/private_output/Homer_motif_enrichment

findMotifsGenome.pl /home/paola/Family1070/private_output/PeakCalling/NKX25/meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.bed hg19r NKX25
findMotifsGenome.pl /home/paola/Family1070/private_output/PeakCalling/ATAC_CM/meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.bed hg19r ATAC_CM 
findMotifsGenome.pl /home/paola/Family1070/private_output/PeakCalling/ATAC_IPSC/meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.bed hg19r ATAC_IPSC 
findMotifsGenome.pl /home/paola/Family1070/private_output/PeakCalling/H3K27AC_CM/meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.bed hg19r H3K27AC_CM
findMotifsGenome.pl /home/paola/Family1070/private_output/PeakCalling/H3K27AC_IPSC/meta_macs2_callPeak_peaks.q001.narrowPeak.collapse.bed hg19r H3K27AC_IPSC


