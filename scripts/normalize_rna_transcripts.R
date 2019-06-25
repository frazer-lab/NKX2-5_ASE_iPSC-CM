home = "/frazer01/home/paola/Family1070/private_output/Validation/data_tables"
setwd(home)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(limma))

rna_cm = read.csv("Production_cm_samples.csv", row.names=1)
raw<-read.table("/projects/CARDIPS/pipeline/RNAseq/combined_files/rsem_expected_counts_isoforms.tsv", header=T, check.names=F, 
                sep="\t", row.names=1)
tpm<-read.table("/projects/CARDIPS/pipeline/RNAseq/combined_files/rsem_tpm_isoforms.tsv", header=T, check.names=F, 
                sep="\t", row.names=1)
gene_info<-read.table("/publicdata/gencode_v19_20151104/gene_info.tsv", header=T, sep="\t")


t_to_g = read.table("/publicdata/gencode_v19_20151104/transcript_to_gene.tsv", header=F, sep="\t")
t_to_g =merge(t_to_g, gene_info, by.x="V2", by.y="gene_id")
t_to_g = t_to_g[match(rownames(raw), t_to_g$V1),]

chromo<-str_split_fixed(t_to_g$chrom, "chr", 2)[,2]
autoso<- as.character(1:22)

raw_sel<-raw[chromo %in% autoso,]
tpm_sel<-tpm[chromo %in% autoso,]

raw_sel<-subset(raw_sel, select=c(as.character(rna_cm$assay_uuid)))
tpm_sel<-subset(tpm_sel, select=c(as.character(rna_cm$assay_uuid)))
raw_sel<-raw_sel[(rowMeans(tpm_sel)>2),]
tpm_sel<-tpm_sel[(rowMeans(tpm_sel)>2),]


counts<-round(raw_sel,0)
workingData<-DESeqDataSetFromMatrix(counts, rna_cm, design= ~ subject_uuid)
workingData<-estimateSizeFactors(workingData)
workingData<-estimateDispersions(workingData, fitType="parametric")
vst_workingData<-varianceStabilizingTransformation(workingData, blind=TRUE)

mat<-as.data.frame(assay(vst_workingData))
colnames(mat)<-as.character(rna_cm$assay_uuid)
write.table(mat, file="Production_cms_vst_counts_transcripts.txt", quote=F, row.names=T, col.names=T, sep="\t")

fit <- lmFit( mat, model.matrix(~ cTNT_TPM, rna_cm))
res <- residuals(fit, mat)
mat_corrected<- res + rowMeans(as.matrix(mat))
write.table(mat_corrected, file="Production_cms_vst_counts_corrected_ctnt_transcripts.txt", quote=F, row.names=T, col.names=T, sep="\t")


