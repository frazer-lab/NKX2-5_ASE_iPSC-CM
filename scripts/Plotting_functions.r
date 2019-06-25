########### ChromHMM plot ##############
########################################

getChromHMM = function(eids, chrom, from, to)
{
    roadmap = list()

    for (eid in eids)
    {
        indata           = read.table(paste("/publicdata/roadmap_15_state_20151104/", eid, "_15_coreMarks_mnemonics_sorted.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE) 
        colnames(indata) = c("chrom", "from", "to", "name")
        this_roadmap     = indata[indata$chrom == chrom & indata$to > from & indata$from <= to,]

		if (length(this_roadmap[this_roadmap$from < from, "from"]) > 0) {this_roadmap[this_roadmap$from < from, "from"] = from}
		if (length(this_roadmap[this_roadmap$to   > to  , "to"  ]) > 0) {this_roadmap[this_roadmap$to   > to  , "to"  ] = to  }
		
		roadmap[[eid]]   = this_roadmap
    }
    return(roadmap)
}

plotChromHMM = function(chrom, from, to, eid_list, eids_file, annots_file)
{
    eids     = read.table(eids_file  , header = TRUE , sep = "\t", stringsAsFactors = FALSE, comment.char = "") 
    annots   = read.table(annots_file, header = TRUE , sep = "\t", stringsAsFactors = FALSE, comment.char = "") 

    annots$name2 = paste("mark", annots$name, sep = ".")
    annots$name2 = gsub("/", ".", annots$name2)

    #annots$chromatin = c("active", "active", "active", "active", "active", "active", "active", "?", "repressed", "?", "?", "?", "repressed", "repressed", "repressed")

    eids2 = eids[eids$eid %in% eid_list, ]

    roadmap = getChromHMM(eids = eids2$eid, chrom = chrom, from = from, to = to) # move into loop

    plot(1, 1, type = "n",
         xlim = c(from, to), ylim = c(0, length(roadmap)),
         xlab = ""         , ylab = ""                   ,
         axes = FALSE
        )
    #axis(1)
    axis(2, at = c(1:length(eids2$eid)) - 0.5, labels = eids2$name, las = 2, tick = FALSE)

    for (ii in 1:length(eids2$eid))
    {
        eid   = eids2[ii, "eid"        ]
        color = eids2[ii, "color"      ]
        desc  = eids2[ii, "description"]
        this  = roadmap[[eid]]

        this  = merge(this, annots[, c("name", "color")])
		
		if (length(this[this$from < from, "from"]) > 0) {this[this$from < from, "from"] = from}
		if (length(this[this$to   > to  , "to"  ]) > 0) {this[this$to   > to  , "to"  ] = to  }
		
        if (length(this$name) > 0)
		{
			rect(xleft = this$from, xright = this$to, ybottom = ii - 0.95, ytop = ii - 0.05, col = this$color, border = this$color)
		}
    }
}


########### plot genes ##################
########################################


plot_genes_all = function(x, from, to, gene_id, exp)
{
    plot(1,1, type = "n", xlim = c(from, to), ylim = c(1, 0), xlab = "", ylab = "", axes = FALSE)

    x = x[x$end >= from & x$start <= to,]
    if(nrow(x)>0) {
    x$color = "#aaaaaa"
    
    if (nrow(x[x$gene_type == "protein_coding", ]) > 0) {x[x$gene_type == "protein_coding", "color"] = "#25285a"}
    if (nrow(x[x$gene_type == "lincRNA"       , ]) > 0) {x[x$gene_type == "lincRNA"       , "color"] = "#8EE5EE"}
    if (nrow(x[x$gene_id   == gene_id         , ]) > 0) {x[x$gene_id   == gene_id         , "color"] = "#ff0000"}
    #if (nrow(x[x$gene_type == "antisense"     , ]) > 0) {x[x$gene_type == "antisense"     , "color"] = "#CDC673"}
    
    if (nrow(x[x$color != "#aaaaaa",]) > 0)
    {
        x      = x[x$color != "#aaaaaa",]
        x$from = x$start
        x$to   = x$end

        if (nrow(x[x$strand == "-",]) > 0)
        {
            x[x$strand == "-", "from"] = x[x$strand == "-", "end"  ]
            x[x$strand == "-", "to"  ] = x[x$strand == "-", "start"]
        }

        x = x[order(x$start),]
        x$y = exp
      offset = (to - from)/20

        for (ii in 1: nrow(x))
        {
            this_y    = exp
            this_from = x[ii, "start"] - offset
            this_to   = x[ii, "end"  ] + offset

            if (ii == 1)
            {
                occupied = data.frame(from = this_from, to = this_to, y = this_y)
            }else
            {
                while(1)
                {
                    if (nrow(occupied[occupied$y == this_y & this_from <= occupied$to,]) > 0)
                    {
                        this_y = this_y + 1
                    }else
                    {
                        x[ii, "y"] = this_y
                        occupied   = rbind(occupied, data.frame(from = this_from, to = this_to, y = this_y))
                        break
                    }
                }
            }
        }
        arrows(x0 = x$from, x1 = x$to, y0 = 1 - (x$y/ max(x$y)), col = x$color, length = 0.05)
        text(x = rowMeans(x[,c("from", "to")]), y = 1 - (x$y/ max(x$y)), labels = x$gene_name, font = 4, pos = 1, col = x$color)
        
    }
        }
}



####### X axis ##########
#########################

plot_xaxis = function(chrom, from, to)
{
    myticks = pretty(seq(from, to, by = (to - from)/ 6))

    plot(1,1, type = "n", xlim = c(from, to), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE)
    abline(h = 1)
    segments(x0 = myticks              , y0 = 1, x1 = myticks, y1 = 0.75)
    text    (x  = myticks              , y  = 0.75, labels = paste(myticks / 1e6, "Mb"), pos = 1)
    text    (x  = (to - from)/ 2 + from, y  = 0.5 , labels = chrom                     , pos = 1)
}



#### Protgtex data from query snps #####
########################################

plot_gtexdata = function(gtex_dir, tissue, snps, from, to, margins, vline) {
if ("gtex" %in% list.files() ) { 
    setwd("gtex") } else {
    dir.create("gtex")
    setwd("gtex")
    }
    
file1 = paste(gtex_dir, tissue, '.v7.signif_variant_gene_pairs.txt', sep="" )
writeLines(snps, "query.snps")
file2 = paste(tissue , 'gtex.snps.queried' , sep=".")
system(paste('grep -f query.snps', file1, '>', file2))

if (file.size(file2)>30) {
fi             = read.table(file2)
writeLines(as.character(unique(fi$V2)), "egenes")
file3          = paste(tissue , 'gtex.snps.all' , sep=".")
system(paste('grep -f egenes', file1, '>', file3))
gtex           = read.table(file3)
colnames(gtex) = gtex_colnames
gtex$logP      = -log(gtex$pval_nominal, 10)
gtex$pos       = as.numeric(str_split_fixed(gtex$variant_id, "_", 5)[,2])
gtex$fm        = gtex$variant_id %in% snps

    
par(mar=margins)    
plot(NA, xlim= c(from, to), ylim= c(0, max(gtex$logP)*1.05), xlab = "", ylab="GTEx -LogP", axes = FALSE)
axis(2)
abline(h = 0, lwd=2)
    
palette = alpha(brewer.pal("Dark2", n=8),0.5)
for (gene in unique(gtex$gene_id)){
    this_gene = subset(gtex, gene_id == gene)
    this_col  = palette[match(gene ,unique(gtex$gene_id))]
    points(logP~pos, this_gene, col=this_col, pch=c(4,19)[this_gene$fm+1])
}
abline(v = vline, lty=2)
center = from + ((to - from) /2)
mtext(side=3,at=from, text=tissue, cex=0.8, adj=0, line=-0.5)   
gene_names = gene_info$gene_name[match(unique(gtex$gene_id), gene_info$gene_id )]
par(mar=c(0,0,0,0))
plot.new()
legend("topleft", pch=19,col =palette[1:length(unique(gtex$gene_id))] , legend = gene_names, bty="n")
#legend("bottom", pch=4,col ="black" , legend = "credible", bty="n")

    }else {
    par(mar=margins)    
    plot(NA, xlim= c(from, to), ylim= c(0, 1), xlab = "", ylab="GTEx -LogP", axes = FALSE)
    axis(2)
    abline(h = 0, lwd=2)
    abline(v = vline, lty=2)
    mtext(side=3,at=from, text=tissue, cex=0.8, adj=0, line=-0.5)
    plot.new()
}
    setwd("../")
    } 


########## Plot Islet eQTL from metaanalysis ##########
#######################################################


plot_islet_eqtls = function(dir, snps, credi_snps, chr, pos,  from, to, margins, vline) {

    if ("islet_eqtls" %in% list.files() ) { 
    setwd("islet_eqtls") } else {
    dir.create("islet_eqtls")
    setwd("islet_eqtls")
    }
    
file1 = paste(dir, 'eQTL.Greenwald_2018_biorxiv.txt', sep="/" )
writeLines(as.character(snps), "query.snps")

file2 = 'eqtl.snps.queried'
system(paste('grep -F -f query.snps', file1, '>', file2))

# if (file.size(file2)>30) {
fi = read.table(file2)
fi = data.frame(str_split_fixed(fi$V1,  "\\-", 2), fi$V6)    
colnames(fi) = c('snps', 'gene_id', 'pvalue')    
ag = aggregate(pvalue < 1e-4 ~ gene_id, fi, sum )
fi = droplevels(subset(fi, gene_id %in%  as.character(ag[ag[,2]>0, 1] )))    

coord = data.frame(snps=snps, chr=chr, pos=pos)
gtex           = merge (fi, coord, by="snps")

    if( nrow(gtex)>0) {
gtex$pvalue[gtex$pvalue==0]= min(gtex$pvalue[gtex$pvalue>0])
gtex$logP      = -log(gtex$pvalue, 10)
gtex$fm        = gtex$snps %in% credi_snps

    
par(mar=margins)    
plot(NA, xlim= c(from, to), ylim= c(0, max(gtex$logP)*1.05), xlab = "", ylab="-LogP", axes = FALSE)
axis(2)
abline(h = 0, lwd=2)
    
palette = alpha(brewer.pal("Dark2", n=8),0.5)
for (gene in unique(gtex$gene_id)){
    this_gene = subset(gtex, gene_id == gene)
    this_col  = palette[match(gene ,unique(gtex$gene_id))]
    points(logP~pos, this_gene, col=this_col, pch=c(4,19)[this_gene$fm+1])
}
abline(v = vline, lty=2)
center = from + ((to - from) /2)
mtext(side=3,at=from, text="Islets", cex=0.8, adj=0, line=-0.5)   
gene_names = gene_info$gene_name[match(unique(gtex$gene_id), gene_info$gene_id_noverson) ]
par(mar=c(0,0,0,0))
plot.new()
legend("topleft", pch=19,col =palette[1:length(unique(gtex$gene_id))] , legend = gene_names, bty="n")
#legend("bottom", pch=4,col ="black" , legend = "credible", bty="n")

    }else {
    par(mar=margins)    
    plot(NA, xlim= c(from, to), ylim= c(0, 1), xlab = "", ylab="-LogP", axes = FALSE)
    axis(2)
    abline(h = 0, lwd=2)
    abline(v = vline, lty=2)
    mtext(side=3,at=from, text="Islets", cex=0.8, adj=0, line=-0.5)
    plot.new()
}
    setwd("../")
    } 





########### GWAS plot ##################
########################################

plot_gwas = function(x, y, from, to, vline, plot_margins, lab="GWAS"){
par(mar=plot_margins)
plot(x, y, ylim=c(0, max(y)*1.05), xlim=c(from,to),
     pch=19, col=alpha("gray20",0.5), cex=1.2,
     axes = FALSE, ylab = "-LogP")
mtext(side=3, at=from, text=lab, cex=0.8, adj=0, line=-0.5)
axis(2)
abline(h = -0.2, lwd=2)
abline(v = vline, lty=2)
}


########### Selex values plot ##################
########################################

get_selex = function (merf, test, from=NA, to=NA, type="only_credible" ) {
if (type == "only_credible"){    
selex = merge(merf, test[,c("varID", 'REF', 'ALT')], by.x = c("varID", 'ref', 'alt'), by.y=1:3)
    }
if (type == "all_interval"){  
    selex = subset(merf, chr == unique(test$CHR) & pos < to & pos > from)
    }
selex$labels = selex$protein
selex = selex[order(selex$variant_p),]
selex$labels[duplicated(selex$varID)] <- NA
selex = selex[order(abs(selex$deltab), decreasing=T),]
lab2 = selex$protein  
lab2[duplicated(selex$varID)] <- NA
selex$labels[is.na(selex$labels)] <- lab2[is.na(selex$labels)]    ## best pvalue and best p per snp
 if( sum(!is.na(selex$labels))> 20){                              ## limit to pbSNP Bonf when many snps in a window
 selex$labels[selex$bonf_p >=0.05 ] <-NA
 selex = selex[order(selex$variant_p),]
 selex$labels[21:nrow(selex) ] <-NA
    }
   
return(selex)
}


plot_selex = function(x, y, from, to, color_values,  colour ="Red", 
                      labels=NA , plot_margins=c(0.5,5,0.5,0), vline=NA, title="Selex") { 
if (length(x) >0) {
par(mar = plot_margins)
my_palette   =  colorRampPalette(brewer.pal(colour, n=9), bias=2) (10) [2:10]
this_palette = alpha(my_palette[as.numeric(cut(color_values, breaks = 9))], 0.8)

plot(x, y, xlim=c(from,to), ylim = c(0, max(y)*1.05),
     pch=19, col= this_palette,  cex=1.2,
     axes = FALSE, ylab = "Delta")

text(x, y,label = labels, pos=4 , cex = 0.8, srt=0)
axis(2)
mtext(side=3,at=from, text=title, cex=0.8, adj=0, line=-0.5)
abline(h = 0, lwd=2)
abline(v = vline, lty=2)
#######
par(mar=c(0,0,0,0))
colkey(my_palette, clim=range(color_values, na.rm=T), clab="-log10(P)\n imbalance", add=F,  length=0.6, width =6, side=1, line.clab=0.5)
    } else {
   par(mar=plot_margins)    
    plot(NA, xlim= c(from, to), ylim= c(0, 1), xlab = "", ylab="Delta", axes = FALSE)
    axis(2)
    abline(h = 0, lwd=2)
    abline(v = vline, lty=2)
    mtext(side=3,at=from, text=title, cex=0.8, adj=0, line=-0.5)
    plot.new()
    }
}

################ layouts ################
#########################################

layout_one_legend_per_row = function(start = 1, n_rows=8, legend_fraction=5){
y = start:(n_rows*2)
i = y%%2==1
m = t (sapply( y[i] , function(x)  c(rep(x,legend_fraction), x+1) ) )
    return(m)
    }
              
layout_one_legend_multiple_rows = function(start = 1, n_rows=3, legend_fraction=5){
y = start:(start+n_rows-1)
m = t(sapply(y , function(x) c(rep(x,legend_fraction), start+n_rows) ) )
    return(m)
     }
             
################ Plot selex individual imbalance ################
#################################################################
            
             
plot_imb = function(snptest, well , allele=NA, type, label )  { 
     
    if(type == "snv"){
      dir ='/projects/T2D/analysis/selex_t2d_t1d/results/preferential_binding'
      colrs= c("green4", "red", "darkgoldenrod1", "blue3")
      names(colrs) = c('A', 'T', 'G', 'C')  
        
}
    
    if(type == "indel"){
        dir ='/projects/T2D/analysis/selex_t2d_t1d/results/preferential_binding_indels/'
        colrs = c("cyan2", "darkmagenta")
        names(colrs) = c("REF", "ALT")  
       
}
    
prot = read.table("/projects/T2D/analysis/selex_t2d_t1d/info/TF_families.txt", header=TRUE, stringsAsFactor=F, sep="\t")
test  = read.table(paste(dir, well, 'Sum_of_alleles_over_cycles.txt', sep="/"), header=TRUE)
pval  = read.table(paste(dir, well, 'pvalues_allele_imbalance.txt', sep="/"), header=TRUE, row.names=1)
   
test = subset(test, snp==as.character(snptest))


    
    
plot(1, xlim=c(0,4), ylim=c(0, max(test[,2:5])), type="n", las=1, xlab="", ylab="",  cex.main=0.8)
        
#    for( n in 2:(length(colrs)+1)){    
#            lines(test$cycle,(test[,n]/test$read_count),  type="b", pch=20, # col=colrs[n-1], lwd=1.5)
#            }
     
     for( n in 2:(length(colrs)+1)) {    
            lines(test$cycle, test[,n],  type="b", pch=20, col=colrs[n-1], lwd=1.5)
            }

    if(is.na(allele)){
    text(1,0.95, paste("p=", signif(pval[as.character(snptest),as.numeric(which.max(test[,n]))],2)), cex=0.8,
                col=colrs[as.numeric(which.max(test[,n]))] )
    } else {
     text(1,0.95, paste("p=", signif(pval[as.character(snptest),as.character(allele)],2)), cex=0.8,
                col=colrs[as.character(allele)] )
 }
    
    mtext(side=3,at=2, text=prot$protein[prot$well==well], cex=0.8, line=0)
    mtext(side=3,at=2, text=label, cex=0.8, line=1)
           }


########### #################################
############################################
writeHeader = function(id,
                      folder,
                      ppn    = 1,
                      que    = "short"
                       )
{
    
   header = paste("#!/bin/bash",
               "",
               paste("#$ -N ", id, sep = ""),
               paste("#$ -o ", folder, "/", id, ".out", sep = ""),
               paste("#$ -e ", folder, "/", id, ".err", sep = ""),
               paste("#$ -l", que ),
               paste("#$ -pe smp", ppn),
                   "#$ -V",
               sep = "\n"
              )
   return(header);
}

             
########### #################################
############################################
             
plot_pbs = function(snptest, well, type , rep="" )  { 
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(DescTools))
    
calc_area = function(xx, yy){
  N   = length(xx)
  dx  = xx[2:N] - xx[1:(N-1)]
  my  = 0.5*(yy[2:N] + yy[1:(N-1)])
  A   = sum(dx * my)
return(A)
    }

   
if(type == "snv"){
      dir ='/projects/T2D/analysis/selex_t2d_t1d/results/obs_pbs_snvs/'
} 
if(type == "indel"){ 
        dir ='/projects/T2D/analysis/selex_t2d_t1d/results/obs_pbs_indels/'  
}
    
prot = read.table("/projects/T2D/analysis/selex_t2d_t1d/info/TF_families.txt", header=TRUE, stringsAsFactor=F, sep="\t")

filename = paste("Sum_of_alleles_over_cycles", rep  , ".txt", sep="")
    test  = read.table(paste(dir, well, filename, sep="/"), header=TRUE)
if(type == "snv") { test = subset(test, snp_var==as.character(snptest))}
if(type == "indel") { test = subset(test, snp==as.character(snptest))}
tab = data.frame( "rest" = ( test$total_reads-test$read_count ),  "reads" = test$read_count )+ 0.1
p_val     = CochranArmitageTest(tab, "increasing")$p.value
tab$logOR = sapply(1:5, function(p) log((tab[1,1] * tab[p,2])/(tab[1,2] * tab[p,1]),10))  ## caclc the log ODDS
obs       = calc_area(1:5, tab$logOR)

tab_ref = data.frame( "rest" = ( test$read_count-test$ref_count ),  "reads" = test$ref_count )+ 0.1
tab_alt = data.frame( "rest" = ( test$read_count-test$alt_count ),  "reads" = test$alt_count )+ 0.1
rownames(tab_ref)= rownames(tab_alt)=test$cycle 
var_p = prop.trend.test(tab_ref$reads, tab_ref$reads + tab_alt$reads , score = c(0,1,1,1,1))$p.value ## this is a chi square test   
logOR_ref = sapply(1:5, function(p) log((tab_ref[1,1] * tab_ref[p,2])/(tab_ref[1,2] * tab_ref[p,1]),10))
logOR_alt = sapply(1:5, function(p) log((tab_alt[1,1] * tab_alt[p,2])/(tab_alt[1,2] * tab_alt[p,1]),10))    
pbs = calc_area(0:4, logOR_ref) - calc_area(0:4, logOR_alt)
    

plot(0:4, tab$logOR,  xlab="", ylab= "log10OR", las=1, cex.axis=0.8)
polygon(c(0:4, 4), c(tab$logOR, 0), col=alpha("gray", 0.5), border = alpha("black", 0.5))
mtext(side=3,at=2, text=prot$protein[prot$well==well], cex=0.8, line=0)
mtext(side=3,at=2, text=snptest, cex=0.8, line=1)
mtext(side=1,at=2, text=paste("p=", signif(p_val, 2), " OBS=", round(obs,2), sep=""), cex=0.7, line=2)

    if(type == "indel") {colors = c("cyan2", "darkmagenta")  }
    
if(type == "snv") { 
    colors = c("green4", "red", "darkgoldenrod1", "blue3") [match(  as.character(t(test[1,c('ref','alt')]) ) , c('A','T','G','C'))]}

plot(c(0:4,0:4), c(logOR_ref, logOR_alt ), pch=rep(c(19, 15), each=5), col=rep(c(colors[1], colors[2]), each=5),
     xlab="", ylab= "log10OR", las=1, cex.axis=0.8)
mtext(side=3,at=2, text=prot$protein[prot$well==well], cex=0.8, line=0)
mtext(side=3,at=2, text=snptest, cex=0.8, line=1)
mtext(side=1,at=2, text=paste("p=", signif(var_p, 2), " PBS=", round(pbs,2), sep=""), cex=0.7, line=2)
   
lines(0:4, logOR_ref, , col=colors[1] )
lines(0:4,  logOR_alt,  col=colors[2])
polygon(c(0:4, 4:0), c(logOR_ref, rev(logOR_alt)), col=alpha("gray", 0.2), border = NA)


    }      
                   
                   
                   
plot_pbs_with_LM = function(snptest, well, type  )  { 
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(DescTools))
    
calc_area = function(xx, yy){
  N   = length(xx)
  dx  = xx[2:N] - xx[1:(N-1)]
  my  = 0.5*(yy[2:N] + yy[1:(N-1)])
  A   = sum(dx * my)
return(A)
    }

   
if(type == "snv"){
      dir ='/projects/T2D/analysis/selex_t2d_t1d/results/obs_pbs_snvs/'
} 
if(type == "indel"){ 
        dir ='/projects/T2D/analysis/selex_t2d_t1d/results/obs_pbs_indels/'  
}
    
prot = read.table("/projects/T2D/analysis/selex_t2d_t1d/info/TF_families.txt", header=TRUE, stringsAsFactor=F, sep="\t")
test  = read.table(paste(dir, well, 'Sum_of_alleles_over_cycles.txt', sep="/"), header=TRUE)
if(type == "snv") { test = subset(test, snp_var==as.character(snptest))}
if(type == "indel") { test = subset(test, snp==as.character(snptest))}
tab = data.frame( "rest" = ( test$total_reads-test$read_count ),  "reads" = test$read_count )+ 0.1
p_val     = CochranArmitageTest(tab, "increasing")$p.value
tab$logOR = sapply(1:5, function(p) log((tab[1,1] * tab[p,2])/(tab[1,2] * tab[p,1]),10))  ## caclc the log ODDS
obs       = calc_area(1:5, tab$logOR)

tab_ref = data.frame( "rest" = ( test$read_count-test$ref_count ),  "reads" = test$ref_count )+ 0.1
tab_alt = data.frame( "rest" = ( test$read_count-test$alt_count ),  "reads" = test$alt_count )+ 0.1
rownames(tab_ref)= rownames(tab_alt)=test$cycle 
#var_p = prop.trend.test(tab_ref$reads, tab_ref$reads + tab_alt$reads , score = c(0,1,1,1,1))$p.value ## this is a chi square test   
logOR_ref = sapply(1:5, function(p) log((tab_ref[1,1] * tab_ref[p,2])/(tab_ref[1,2] * tab_ref[p,1]),10))
logOR_alt = sapply(1:5, function(p) log((tab_alt[1,1] * tab_alt[p,2])/(tab_alt[1,2] * tab_alt[p,1]),10))    
pbs = calc_area(0:4, logOR_ref) - calc_area(0:4, logOR_alt)
    

plot(0:4, tab$logOR,  xlab="", ylab= "log10OR", las=1, cex.axis=0.8)
polygon(c(0:4, 4), c(tab$logOR, 0), col=alpha("gray", 0.5), border = alpha("black", 0.5))
#mtext(side=3,at=2, text=prot$protein[prot$well==well], cex=0.8, line=0)
#mtext(side=3,at=2, text=snptest, cex=0.8, line=1)
#mtext(side=1,at=2, text=paste("p=", signif(p_val, 2), " OBS=", round(obs,2), sep=""), cex=0.7, line=2)

    if(type == "indel") {colors = c("cyan2", "darkmagenta")  }
    
if(type == "snv") { 
    colors = c("green4", "red", "darkgoldenrod1", "blue3") [match(  as.character(t(test[1,c('ref','alt')]) ) , c('A','T','G','C'))]}

plot(c(0:4,0:4), c(logOR_ref, logOR_alt ), pch=rep(c(19, 15), each=5), col=rep(c(colors[1], colors[2]), each=5),
     xlab="", ylab= "log10OR", las=1, cex.axis=0.8)
#mtext(side=3,at=2, text=prot$protein[prot$well==well], cex=0.8, line=0)
#mtext(side=3,at=2, text=snptest, cex=0.8, line=1)
mtext(side=3,at=2, text=paste(" PBS=", round(pbs,2), sep=""), cex=0.7, line=2)

#var_p =  padjtab [snptest,1]
#co_p4 =  padjtab [snptest,2]
#co_p1 =  padjtab [snptest,3]
    
#mtext(side=3,at=2, text=paste("FDR_trend 0-4=", signif(co_p4, 2)),  cex=0.7, line=0) 
#mtext(side=3,at=2, text=paste("FDR_trend 0-1=", signif(co_p1, 2)),  cex=0.7, line=1) 
#mtext(side=1,at=2, text=paste("FDR_permu=", signif(var_p, 2), " PBS=", round(pbs,2), sep=""), cex=0.7, line=2)
    
    
lines(0:4, logOR_ref, , col=colors[1] )
lines(0:4,  logOR_alt,  col=colors[2])
polygon(c(0:4, 4:0), c(logOR_ref, rev(logOR_alt)), col=alpha("gray", 0.2), border = NA)

    
##### LM model########
se_ref = sapply(1:5, function(p) sqrt((1/tab_ref[1,1]) + (1/tab_ref[p,2]) + (1/tab_ref[1,2]) + (1/tab_ref[p,1]) ))
se_alt = sapply(1:5, function(p) sqrt((1/tab_alt[1,1]) + (1/tab_alt[p,2]) + (1/tab_alt[1,2]) + (1/tab_alt[p,1]) ))
    
z_ref  =   logOR_ref/se_ref
z_alt  =   logOR_alt/se_alt


dz   = z_ref-z_alt

    cycles = 0:4
    l = lm(dz ~ cycles -1)

plot(cycles, dz, ylab = "d_zscores", las=1, pch=21, bg="gray80")    
    abline(l, col="red")
mtext(side=1, at = 2, line=2, text = paste("b=",round(summary(l)[[4]][1],3)), cex=0.7)
mtext(side=3, at = 2, line=0, text = paste("r=",round(summary(l)$r.squared,2 ),"p=",signif(summary(l)[[4]][4],2)) , cex=0.7) 

    }      