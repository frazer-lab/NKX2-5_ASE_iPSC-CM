#module load R/3.1.1 (R version 3.1.1 (2014-07-10) -- "Sock it to Me")

library("TFMPvalue")
library("seqinr")
library("stringr")
library("Biostrings")
library("seqLogo")


pwms_files <- list.files("PWMs/")
homerScore<- read.table("HomerScores.txt", row.names=1)
seq <- read.table("../TestMAtr.txt", header = T, sep = "\t")

#seq=seq[1:100,] #test first few lines

seq$REfSeq <- as.character(seq$REfSeq)
seq$ALtSeq <- as.character(seq$ALtSeq)

TFMaxscores <- data.frame()


i <- 1
for (f in 1:length(pwms_files)) {


mot <- read.table(paste("PWMs/", pwms_files[f],sep=""), header = T)
mat <- t(mot)

#NucFreq <- c(0.295, 0.205, 0.205, 0.295) #erin's original parameters
NucFreq <- c(0.25, 0.25, 0.25, 0.25) #homer parameters
pwm <- log(mat/NucFreq)

#################################### 
# Method 1: Homer Scores 
#################################### 
tfm_score<-homerScore[pwms_files[f],1]
#################################### 
# Method 2: Haploreg 
#################################### 
#tfm_score <- TFMLazyScore(pwm,  0.00000004, bg = c(A=0.25, C=0.25, G=0.25,T=0.25), type = "PWM")


for (s in 1:length(seq$REfSeq)) {
	maxRefScore <- 0
	maxRefScorePos <- 0
	maxRefScoreAllele = NA
	maxRefScoreStrand = NA

	maxAltScore <- 0
	maxAltScorePos <- 0
	maxAltScoreAllele = NA
	maxAltScoreStrand = NA


	for (n in 1:(nchar(seq$REfSeq[s])-length(pwm[1,]))) {
		sub <- substr(seq$REfSeq[s], n, n-1+length(pwm[1,]))
		score <- 0
		for (nuc in 1:length(pwm[1,])) {
			score <- score + pwm[substr(sub,nuc,nuc),nuc]
		}

		if (score > maxRefScore) {
			maxRefScore = score
			maxRefScorePos = n
			maxRefScoreAllele = "Ref"
			maxRefScoreStrand = "+"

		}
	}
	revComp <- paste(toupper(comp(unlist(strsplit(seq$REfSeq[s], "") ))), sep = "", collapse = "")
	for (n in 1:(nchar(revComp)-length(pwm[1,]))) {
		sub <- substr(revComp, n, n-1+length(pwm[1,]))
		score <- 0
		for (nuc in 1:length(pwm[1,])) {
			score <- score + pwm[substr(sub,nuc,nuc),nuc]
		}

		if (score > maxRefScore) {
			maxRefScore = score
			maxRefScorePos = n
			maxRefScoreAllele = "Ref"
			maxRefScoreStrand = "-"
		}


	}
	for (n in 1:(nchar(seq$ALtSeq[s])-length(pwm[1,]))) {
		sub <- substr(seq$ALtSeq[s], n,n-1+length(pwm[1,]))
		score <- 0
		for (nuc in 1:length(pwm[1,])) {
			score <- score + pwm[substr(sub,nuc,nuc),nuc]
		}

		if (score > maxAltScore) {
			maxAltScore = score
			maxAltScorePos = n
			maxAltScoreAllele = "Alt"
			maxAltScoreStrand = "+"
			
		}
	}

	revComp <- paste(toupper(comp(unlist(strsplit(seq$ALtSeq[s], "") ))), sep = "", collapse = "")
	for (n in 1:(nchar(revComp)-length(pwm[1,]))) {
		sub <- substr(revComp, n, n-1+length(pwm[1,]))
		score <- 0
		for (nuc in 1:length(pwm[1,])) {
			score <- score + pwm[substr(sub,nuc,nuc),nuc]
		}

		if (score > maxAltScore) {
			maxAltScore = score
			maxAltScorePos = n
			maxAltScoreAllele = "Alt"
			maxAltScoreStrand = "-"
		}

		
	}
		TFMaxscores[i,"Motif"] <- pwms_files[f]
		TFMaxscores[i,"TFM_score"] <- tfm_score
		TFMaxscores[i,"Seq"] <- s
		TFMaxscores[i,"SNP"] <- seq$Variant[s]
		TFMaxscores[i,"maxRefScore"] <- maxRefScore 
		TFMaxscores[i,"maxRefScorePos"] <- maxRefScorePos
		TFMaxscores[i,"maxRefScoreAllele"] <- maxRefScoreAllele
		TFMaxscores[i,"maxRefScoreStrand"] <- maxRefScoreStrand

		TFMaxscores[i,"maxAltScore"] <- maxAltScore 
		TFMaxscores[i,"maxAltScorePos"] <- maxAltScorePos
		TFMaxscores[i,"maxAltScoreAllele"] <- maxAltScoreAllele
		TFMaxscores[i,"maxAltScoreStrand"] <- maxAltScoreStrand
		
		i <- i+1

}
}

write.table(TFMaxscores, "All_siteScores_1.txt", quote = F, sep = "\t", row.names = F)

sigSites <- subset(TFMaxscores, (TFMaxscores$maxRefScore >= TFMaxscores$TFM_score |  TFMaxscores$maxAltScore >= TFMaxscores$TFM_score) & abs(TFMaxscores$maxRefScore - TFMaxscores$maxAltScore) > 0 )

write.table(sigSites, "AlleleAlteringSites_1.txt", quote = F, sep = "\t", row.names = F)






################## ~~~~~~~ ########################## ~~~~~~~~~~~~ ####################
##for each site, calculate the Ref and Alt score at the same positions
################## ~~~~~~~ ########################## ~~~~~~~~~~~~ ####################

for (s in 1:length(sigSites$Seq)) {


f <- match(sigSites$Motif[s], pwms_files)

mot <- read.table(paste("PWMs/", pwms_files[f],sep=""), header = T)
mat <- t(mot)
#NucFreq <- c(0.295, 0.205, 0.205, 0.295)
NucFreq <- c(0.25, 0.25, 0.25, 0.25) 
pwm <- log(mat/NucFreq)


#################################### 
# Method 1: Homer Scores 
#################################### 
tfm_score<-homerScore[pwms_files[f],1]
#################################### 
# Method 2: Haploreg 
#################################### 
#tfm_score <- TFMLazyScore(pwm,  0.0000004, bg = c(A=0.25, C=0.25, G=0.25,T=0.25), type = "PWM")

##ref allele at maxAlt Site
if (sigSites$maxRefScore[s] > 0) {
	query_seq <- seq$ALtSeq[sigSites$Seq[s]]

	if (sigSites$maxRefScoreStrand[s] == "-") {
		query_seq <-  paste(toupper(comp(unlist(strsplit(seq$ALtSeq[sigSites$Seq[s]], "") ))), sep = "", collapse = "")
	}

	sub <- substr(query_seq, sigSites$maxRefScorePos[s],sigSites$maxRefScorePos[s]-1+length(pwm[1,]))
	sigSites[s,"Seq_atMaxRefSite"] <- sub

	score <- 0

	for (nuc in 1:length(pwm[1,])) {
		score <- score + pwm[substr(sub,nuc,nuc),nuc]
	}

	sigSites[s,"AltScore_atMaxRefSite"] <- score
}

##alt

if (sigSites$maxAltScore[s] > 0) {
	query_seq <- seq$REfSeq[sigSites$Seq[s]]

	if (sigSites$maxAltScoreStrand[s] == "-") {
		query_seq <-  paste(toupper(comp(unlist(strsplit(seq$REfSeq[sigSites$Seq[s]], "") ))), sep = "", collapse = "")
	}

	sub <- substr(query_seq, sigSites$maxAltScorePos[s],sigSites$maxAltScorePos[s]-1+length(pwm[1,]))
	sigSites[s,"Seq_atMaxAltSite"] <- sub
	score <- 0

	for (nuc in 1:length(pwm[1,])) {
		score <- score + pwm[substr(sub,nuc,nuc),nuc]
	}
	sigSites[s,"RefScore_atMaxAltSite"] <- score

}
}
write.table(sigSites, "AlleleAlteringSites_scores_1.txt", quote = F, sep = "\t", row.names = F)

################## ~~~~~~~ ########################## ~~~~~~~~~~~~ ####################
## PLOTS SUMMARY AT ALL SIGNIFICANT SITES
################## ~~~~~~~ ########################## ~~~~~~~~~~~~ ####################

# sigSites$Motif <- factor(sigSites$Motif)
# 
# up=NA
# down=NA
# for (f in 1:length(levels(sigSites$Motif))) {
#   mo=sigSites[sigSites$Motif==levels(sigSites$Motif)[f],]
#   up[f]<-nrow(mo[mo$maxRefScore > mo$maxAltScore,])
#   down[f]<-nrow(mo[mo$maxRefScore < mo$maxAltScore,])
# }
# 
# tab=rbind(up,down)
# colnames(tab)=levels(sigSites$Motif)
# #tab=(tab/nrow(usSites))*100
# 
# tab[2,]=-tab[2,]
# tab=tab[,order(tab[1,], decreasing=T)]
# 
# pdf("SummaryNumberofMotifs_1.pdf")
# par(pin=c(3,3))
# barplot(as.numeric(tab[1,]), beside=TRUE, names.arg=colnames(tab), 
#         las=2, ylim=c(min(tab), max(tab)), col="orangered")
# barplot(as.numeric(tab[2,]), beside=TRUE, add=T, axes=FALSE, col="dodgerblue3")
# dev.off()
# 
# ################## ~~~~~~~ ########################## ~~~~~~~~~~~~ ####################
# ## UNIQUE SIGNIFICANT SITES = keep SNPs where difference between the Ref and Alt Score is higher
# ################## ~~~~~~~ ########################## ~~~~~~~~~~~~ ####################
# 
# usSites= sigSites[order((abs(sigSites$maxRefScore-sigSites$maxAltScore)), decreasing=TRUE),]
# usSites=subset(usSites, !duplicated(usSites$SNP))
# usSites[is.na(usSites)]<-0
# write.table(usSites, "UniqAlteringSites_1.txt", quote = F, sep = "\t", row.names = F)
# 
# ###
# up=NA
# down=NA
# for (f in 1:length(levels(sigSites$Motif))) {
#   mo=usSites[usSites$Motif==levels(sigSites$Motif)[f],]
#   up[f]<-nrow(mo[mo$maxRefScore > mo$maxAltScore,])
#   down[f]<-nrow(mo[mo$maxRefScore < mo$maxAltScore,])
# }
# 
# tab=rbind(up,down)
# colnames(tab)=levels(sigSites$Motif)
# #tab=(tab/nrow(usSites))*100
# 
# tab[2,]=-tab[2,]
# tab=tab[,order(tab[1,], decreasing=T)]
# 
# pdf("SummaryNumberodMotifs_unique_1.pdf")
# par(pin=c(3,3))
# barplot(as.numeric(tab[1,]), beside=TRUE, names.arg=colnames(tab), 
#         las=2, ylim=c(min(tab), max(tab)), col="orangered")
# barplot(as.numeric(tab[2,]), beside=TRUE, add=T, axes=FALSE, col="dodgerblue3")
# dev.off()
# 

################## ~~~~~~~ ########################## ~~~~~~~~~~~~ ####################
## LOGOS
################## ~~~~~~~ ########################## ~~~~~~~~~~~~ ####################
# ##### Seq_atMaxRefSite is the sequence with the alterante (minor) allele and vice-versa####
# 
# pdf("Logos_hisograms_1.pdf")
# par(mar=c(10,10,10,10), pin=c(2,2))
# 
# for (f in 1:length(levels(sigSites$Motif))) {
#   mo <- subset(sigSites, Motif == levels(sigSites$Motif)[f])
#   mo=mo[(mo$maxRefScore > mo$maxAltScore), ]
#   sequence_maj= mo$Seq_atMaxAltSite
#   position=16-mo$maxRefScorePos
#   consensus_maj=consensusMatrix(sequence_maj, as.prob=T)
#   mo=mo[mo$maxRefScorePos==mo$maxAltScorePos,]
#   sequence_min=mo$Seq_atMaxRefSite
#   consensus_min=consensusMatrix(sequence_min, as.prob=T)
#   seqLogo(consensus_maj, yaxis=F)
#   seqLogo(consensus_min, yaxis=F)
#   hist(position, breaks=seq(0,ncol(consensus_maj),by=1), xaxt='n', col="gray", 
#        freq = TRUE, main=levels(sigSites$Motif)[f])
#   axis(side=1, at=seq(0.5,ncol(consensus_maj)-0.5), labels=seq(1,ncol(consensus_maj)))
# }
# 
# 
# dev.off()
# 
# 
# pdf("Logos_hisograms_unique_1.pdf")
# par(mar=c(10,10,10,10), pin=c(2,2))
# 
# for (f in 1:length(levels(sigSites$Motif))) {
#   mo <- subset(usSites, Motif == levels(sigSites$Motif)[f])
#   mo=mo[(mo$maxRefScore > mo$maxAltScore), ]
#   sequence_maj= mo$Seq_atMaxAltSite
#   position=16-mo$maxRefScorePos
#   consensus_maj=consensusMatrix(sequence_maj, as.prob=T)
#   mo=mo[mo$maxRefScorePos==mo$maxAltScorePos,]
#   sequence_min=mo$Seq_atMaxRefSite
#   consensus_min=consensusMatrix(sequence_min, as.prob=T)
#   seqLogo(consensus_maj, yaxis=F)
#   seqLogo(consensus_min, yaxis=F)
#   hist(position, breaks=seq(0,ncol(consensus_maj),by=1), xaxt='n', col="gray", 
#        freq = TRUE, main=levels(sigSites$Motif)[f])
#   axis(side=1, at=seq(0.5,ncol(consensus_maj)-0.5), labels=seq(1,ncol(consensus_maj)))
# }
# 
# 
# dev.off()
