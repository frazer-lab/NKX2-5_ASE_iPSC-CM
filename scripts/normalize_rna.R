home = "/frazer01/home/paola/Family1070/private_output/Validation/data_tables"
setwd(home)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(limma))

main = read.table ("/frazer01/projects/CARDIPS/data/CARDiPS_production_metadata/assay_mastertable.txt", sep= "\t", header=T)
main= main[!duplicated(main),]
rna_cm = droplevels(subset(main, dataset=="production" &assay =="rna" & cell=="CM"))
raw<-read.table("/projects/CARDIPS/pipeline/RNAseq/combined_files/rsem_expected_counts.tsv", header=T, check.names=F, 
                sep="\t", row.names=1)
tpm<-read.table("/projects/CARDIPS/pipeline/RNAseq/combined_files/rsem_tpm.tsv", header=T, check.names=F, 
                sep="\t", row.names=1)

gene_info<-read.table("/publicdata/gencode_v19_20151104/gene_info.tsv", header=T, sep="\t")

chromo<-str_split_fixed(gene_info$chrom, "chr", 2)[,2]
unique(chromo)
autoso<- as.character(1:22)

raw_sel<-raw[chromo %in% autoso,]
tpm_sel<-tpm[chromo %in% autoso,]

raw_sel<-subset(raw_sel, select=c(as.character(rna_cm$assay_uuid)))
tpm_sel<-subset(tpm_sel, select=c(as.character(rna_cm$assay_uuid)))
raw_sel<-raw_sel[(rowMeans(tpm_sel)>2),]
tpm_sel<-tpm_sel[(rowMeans(tpm_sel)>2),]

rna_cm<-rna_cm[,c(6:length(rna_cm), 1:5)]
rna_cm<-rna_cm[as.character(rna_cm$assay_uuid) %in% colnames(raw),]
ctnt<- tpm_sel[as.character(gene_info[gene_info$gene_name=="TNNT2","gene_id"]),]
ctnt<- subset(ctnt, select=as.character(rna_cm$assay_uuid))
rna_cm$cTNT_TPM = t(ctnt)[,1]

write.csv(rna_cm, "Production_cm_samples.csv")


counts<-round(raw_sel,0)
workingData<-DESeqDataSetFromMatrix(counts, rna_cm, design= ~ subject_uuid)
workingData<-estimateSizeFactors(workingData)
workingData<-estimateDispersions(workingData, fitType="parametric")
vst_workingData<-varianceStabilizingTransformation(workingData, blind=TRUE)

mat<-as.data.frame(assay(vst_workingData))
colnames(mat)<-as.character(rna_cm$assay_uuid)
write.table(mat, file="Production_cms_vst_counts.txt", quote=F, row.names=T, col.names=T, sep="\t")

fit <- lmFit( mat, model.matrix(~ cTNT_TPM, rna_cm))
res <- residuals(fit, mat)
mat_corrected<- res + rowMeans(as.matrix(mat_cm))
write.table(mat, file="Production_cms_vst_counts_corrected_ctnt.txt", quote=F, row.names=T, col.names=T, sep="\t")


