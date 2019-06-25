### Define tools
samtools   = "/software/samtools-1.2/samtools"
bcftools   = "/software/bcftools-1.2/bcftools"
python     = "/frazer01/home/matteo/anaconda2/bin/python"
java       = "/usr/lib/jvm/java-1.7.0/bin/java"
rscript    = "/frazer01/home/matteo/software/R-3.2.2/bin/Rscript"

### Get table from the database
getTableFromDb = function(query, add_dash = TRUE)
{    
    con <- dbConnect(MySQL(), user="cardips", password="G00dC@kes", dbname="cardips", host="localhost")

    rs <- dbSendQuery(con, query)
    db_table <- fetch(rs, n = -1)

    dbClearResult(rs)
    dbDisconnect(con)
    
    if (add_dash == TRUE)
    {
        db_table$assay_id = unlist(lapply(db_table$assay_id, addDashUUID))
        db_table$id       = unlist(lapply(db_table$id      , addDashUUID))
        db_table$wgs_id   = unlist(lapply(db_table$wgs_id  , addDashUUID))
    }

    return(db_table)
}

### Calculate counts
### bed: bed file with peaks
### bam_files: list of bam files to calculate counts on
getCountsFromPeaks = function(id, bed, bam_files,
                              python        = "/frazer01/home/matteo/anaconda2/bin/python",
                              featureCounts = "/software/subread-1.5.0-Linux-x86_64/bin/featureCounts"
                             )
{
    saf         = sub(".bed", ".saf"   , bed)
    output_file = sub(".bed", ".counts", bed)
    ppn         = 8
    
    head= paste("#!/bin/bash",
                "",
                paste("#\\$ -N", id),
                paste("#\\$ -o ", getwd(), "/", id, ".out", sep = ""),
                paste("#\\$ -e ", getwd(), "/", id, ".err", sep = ""),
                "#\\$ -l week",
                paste("#\\$ -pe smp", ppn),
                sep = "\n"
               )
    
    change_folder = paste("cd", getwd())
    
    python1 = paste(python, "/frazer01/home/cdeboever/repos/cdeboever3/cdpipelines/cdpipelines/convert_bed_to_saf.py",
                    bed,
                    saf,
                    sep = " \\\n"
                   )
    
    bam_concat = paste(bam_files, collapse = " ")
    
    feature_counts1 = paste(featureCounts, "-p -T 10 -F SAF --donotsort",
                            paste("-a", saf),
                            paste("-o", output_file),
                            bam_concat,
                            sep = " \\\n"
                           )

    out = paste(head, 
                change_folder,
                python1,
                feature_counts1,
                sep = "\n\n"
               )
    
    return(out);
}

### Substitute the bam file name with uuid in the counts file header
subUuidsInHeader = function(out_table, counts_file)
{
    out_table_bams = out_table[out_table$bam_file != "", ]

    #counts_file = paste("chipseq_peaks/", chipseq_folder, ".counts", sep = "") # move out of function
    out_counts  = gsub (".counts", ".uuid.counts", counts_file)

    counts = read.table(counts_file, header = TRUE , sep = "\t", stringsAsFactors = FALSE) 

    for (ii in 1: length(out_table_bams$name))
    {
        bam_file = out_table_bams$bam_file[[ii]]
        assay_id = out_table_bams$assay_id[[ii]]

        bam_file2 = paste("X", gsub("/", ".", bam_file), sep = "")
        bam_file3 = gsub("-", ".", bam_file2)

        colnames(counts) = gsub(bam_file3, assay_id, colnames(counts))
    }

    write.table(counts, file = out_counts, append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
}

### Write header of SH file
writeHeader = function(id, folder, log_folder, que, ppn)
{
    header = paste("#!/bin/bash",
                "",
                paste("#$ -N job.", id, sep = ""),
                paste("#$ -o ", log_folder, "/", id, ".out", sep = ""),
                paste("#$ -e ", log_folder, "/", id, ".err", sep = ""),
                paste("#$ -l", que),
                paste("#$ -pe smp", ppn),
                 "export PATH=/frazer01/home/paola/anaconda2/bin:$PATH",
				 "export PATH=/frazer01/home/paola/anaconda2/bin/python:$PATH",
                #"module load bedtools/2.25.0",
                #"module load tabix",
                #"module load STAR",
                #"module load picard",
                #"module load GATK",
                #"module load samtools",
                #"module load ucsc.linux",
                #"source activate py27",
				"module load cardips",
				"source activate cardips",
                sep = "\n"
               )
			   
	outdir  = paste(folder, "wasp", sep = "/")
	
    create_outdir = paste(paste("mkdir -p", outdir),
                          paste("cd"      , outdir),
                          sep = "\n"
                         )

	out = paste(header, 
				create_outdir,
				sep = "\n\n"
			   )
			   
    return(out);
}

### Merge BAM files
mergeBams = function(bam_files, subject_dir, out_bam, ppn,
                     samtools = "/software/samtools-1.2/samtools",
					 java     = "/usr/lib/jvm/java-1.7.0/bin/java"
                    )
{
	start_bam = sub("bam", "no_rg.bam", out_bam)
	
	java_mem = paste("-Xmx", round(ppn * 4 * 0.75, digits = 0), "g", sep = "")
	
    #rsync_command = ""
    #for (bam_file in bam_files)
    #{
    #    rsync_command = paste(rsync_command, 
    #                          paste("rsync -avz",
    #                                bam_file,
    #                                paste(subject_dir, ".", sep = "/"),
    #                                "&"
    #                               ),
    #                          sep = "\n"
    #                         )
    #}
    #rsync_command = paste(rsync_command, "wait;", sep = "\n")
	
    # 2. merge bam files
	if (length(bam_files) == 1)
	{
		samtools_merge_command = paste("rsync -avz",
									   bam_files[[1]],
									   start_bam
									  )
	}else
	{
		samtools_merge_command = paste(paste(java, java_mem, "-jar"), 
									   "$picard MergeSamFiles",
									   paste("I=", bam_files, sep = "", collapse = " "),
									   paste("O=", start_bam, sep = ""),
									   sep = " \\\n"
									  )

		#samtools_merge_command = paste(samtools,
		#							   "merge -r",
		#							   start_bam,
		#							   paste(subject_dir, "*.bam", sep = "/")
		#							  )
    }
	
	add_rg_command = paste(paste(java, java_mem, "-jar"),
						   "$picard AddOrReplaceReadGroups",
						   paste("I=", start_bam, sep = ""),
						   paste("O=", out_bam  , sep = ""),
						   "RGID=4",
						   "RGLB=lib1",
						   "RGPL=illumina",
						   "RGPU=unit1",
						   "RGSM=20",
						   sep = " \\\n"
						  )
	
	out = paste("# Merge BAMs", 
				#rsync_command,
				samtools_merge_command,
				add_rg_command,
				sep = "\n\n"
			   )
			   
    return(out)
}
### from Chris: WASP allele swap
waspAlleleSwap = function(id, wgs_id, folder, bam, vcf, bed, is.gz = FALSE,
                          samtools   = "/software/samtools-1.2/samtools", 
                          genome_fai = "/publicdata/gatk_bundle_2.8/hg19/ucsc.hg19.fasta.fai", 
                          bcftools   = "/software/bcftools-1.2/bcftools", 
                          chrom_conv = "/repos/cardips-pipelines/RNA/chrom_conv.tsv",
                          python     = "/frazer01/home/matteo/anaconda2/bin/python"
                         )
{
    new_vcf = paste(id, "input.vcf", sep = ".")
	outdir  = paste(folder, "wasp", sep = "/")
    
    if (is.gz == FALSE)
	{
		gzip1 = paste(paste("rsync -az", vcf, new_vcf),
					  paste("bgzip", new_vcf),
					  paste("tabix -p vcf ", new_vcf, ".gz", sep = ""),
					  sep = "\n")
	}else
	{
		gzip1 = paste(paste("rsync -az", vcf                         , paste(new_vcf, "gz"    , sep = ".")),
					  paste("rsync -az", paste(vcf, "tbi", sep = "."), paste(new_vcf, "gz.tbi", sep = ".")),
					  sep = "\n")
	}
    
    python1 = paste(#paste(python, "/frazer01/home/cdeboever/repos/cdeboever3/cdpipelines/cdpipelines/make_wasp_input.py"),
                    "make_wasp_input",
                    paste(outdir, "/", id, "_hets.vcf", sep = ""),
                    wgs_id,
                    paste(outdir, "snps", sep = "/"),
                    bed,
                    paste("-v", paste(new_vcf, "gz", sep = ".")),
                    paste("-g", genome_fai),
                    paste("-b", bcftools  ),
                    paste("-t", outdir    ),
                    paste("-c", chrom_conv),
                    sep = " \\\n"
                   )
    
    samtools1 = paste(paste(samtools, "view -b -F 1024"), # eliminated "-q 255" because ChIP-seq gives zero otherwise
                      bam,
                      paste("> ", outdir, "/", id, "_uniq.bam", sep = ""),
                      sep = " \\\n"
                     )
    
    python2 = paste(paste(python, "/repos/cardips-pipelines/submodules/WASP/mapping/find_intersecting_snps.py -s -p"), 
                    paste(outdir, "/", id, "_uniq.bam", sep = ""),
                    paste(outdir, "snps", sep = "/"),
                    sep = " \\\n"
                   )
    
    rm_uselsess = paste("rm -r", 
                        paste(outdir, "/", id, "_uniq_keep.bam", sep = ""),
                        paste(outdir, "snps", sep = "/"),
                        sep = " \\\n"
                       )

    out = paste("# waspAlleleSwap",
                gzip1,
                python1,
                samtools1,
                python2,
                #rm_uselsess,
                sep = "\n\n"
               )
    
    return(out);
}

### From Chris: WASP remap
waspRemap = function(id, folder, ppn,
                     genomeDir = "/publicdata/star_index_hg19_sorted_20151123"
                    )
{
    outdir = paste(folder, "wasp", sep = "/")
    
    star1  = paste("STAR",
                   paste("--genomeDir", genomeDir),
                   "--genomeLoad NoSharedMemory",
                   #"--genomeLoad LoadAndKeep",
                   "--readFilesCommand zcat",
                   paste("--readFilesIn", paste(outdir, "/", id, "_uniq.remap.fq1.gz", sep = ""), paste(outdir, "/", id, "_uniq.remap.fq2.gz", sep = "")),
                   "--outSAMattributes All",
                   "--outSAMunmapped Within",
                   paste("--outSAMattrRGline ID:1 PL:ILLUMINA PU:150522_D00611_0122_AC6G7GANXX LB:", id, " SM:", id, sep = ""),
                   "--outFilterMultimapNmax 20",
                   "--outFilterMismatchNmax 999",
                   "--alignIntronMin 20",
                   "--alignIntronMax 1000000",
                   "--alignMatesGapMax 1000000",
                   "--outSAMtype BAM Unsorted",
                   sep = " \\\n"
                  )
    
    if_cond = "if [ -d _STARtmp ] ; then rm -r _STARtmp ; fi"
    
    move_files = paste(paste("mv Aligned.out.bam" , paste(outdir, "/", id, ".bam"             , sep = "")),
                       paste("mv Log.out"         , paste(outdir, "/", id, "_Log.out"         , sep = "")),
                       paste("mv Log.final.out"   , paste(outdir, "/", id, "_Log.final.out"   , sep = "")),
                       paste("mv Log.progress.out", paste(outdir, "/", id, "_Log.progress.out", sep = "")),
                       paste("mv SJ.out.tab"      , paste(outdir, "/", id, "_SJ.out.tab"      , sep = "")),
                       sep = "\n"
                      )
    
    rm_uselsess = paste("rm -r", 
                        paste(outdir, "/", id, "_uniq.remap.fq1.gz", sep = ""),
                        paste(outdir, "/", id, "_uniq.remap.fq2.gz", sep = ""),
                        paste(outdir, "/", id, "_SJ.out.tab"       , sep = ""),
                        sep = " \\\n"
                       )
    
    out = paste("# waspRemap",
                star1,
                if_cond,
                move_files,
                #rm_uselsess,
                sep = "\n\n"
               )
    
    return(out);
}  

### From Chris: WASP compare alignments
waspAlignmentCompare = function(id, folder, ppn,
                          samtools   = "/software/samtools-1.2/samtools", 
                          genome     = "/publicdata/gatk_bundle_2.8/hg19/ucsc.hg19.fasta", 
                          java       = "/usr/lib/jvm/java-1.7.0/bin/java",
                          python     = "/frazer01/home/matteo/anaconda2/bin/python"
                         )
{
    outdir = paste(folder, "wasp", sep = "/")
    
    java_mem = paste("-Xmx", round(ppn * 4 * 0.75, digits = 0), "g", sep = "")
    
    python1 = paste(paste(python, "/repos/cardips-pipelines/submodules/WASP/mapping/filter_remapped_reads.py -p"),
                    paste(outdir, "/", id, "_uniq.to.remap.bam"   , sep = ""),
                    paste(outdir, "/", id, ".bam"                 , sep = ""),
                    paste(outdir, "/", id, "_filtered.bam"        , sep = ""),
                    paste(outdir, "/", id, "_uniq.to.remap.num.gz", sep = ""),
                    sep = " \\\n"
                   )
    
    java1 = paste(paste(java, java_mem, "-jar"),
                  paste("-XX:ParallelGCThreads=", ppn, sep = ""),
                  paste("-Djava.io.tmpdir", outdir, sep = ""),
                  "-jar $picard SortSam",
                  "VALIDATION_STRINGENCY=SILENT",
                  paste("I=", outdir, "/", id, "_filtered.bam", sep = ""),
                  paste("O=", outdir, "/", id, "_sorted.bam"  , sep = ""),
                  "SO=coordinate",
                  "CREATE_INDEX=TRUE",
                  sep = " \\\n"
                 )
    
    move1 = paste("mv",
                  paste(outdir, "/", id, "_sorted.bai"    , sep = ""),
                  paste(outdir, "/", id, "_sorted.bam.bai", sep = "")
                 )
    
    java2 = paste(paste(java, java_mem, "-jar"),
                  paste("-XX:ParallelGCThreads=", ppn, sep = ""),
                  paste("-Djava.io.tmpdir", outdir, sep = ""),
                  "-jar $picard ReorderSam",
                  "VALIDATION_STRINGENCY=SILENT",
                  paste("I=", outdir, "/", id, "_sorted.bam"   , sep = ""),
                  paste("O=", outdir, "/", id, "_reordered.bam", sep = ""),
                  paste("R=", genome, sep = ""),
                  sep = " \\\n"
                 )
    
    samtools_index = paste(samtools, "index", paste(outdir, "/", id, "_reordered.bam", sep = ""))
    
    java3 = paste(paste(java, java_mem, "-jar"),
                  paste("-XX:ParallelGCThreads=", ppn, sep = ""),
                  paste("-Djava.io.tmpdir", outdir, sep = ""),
                  "-jar $GATK",
                  paste("-R", genome),
                  "-T ASEReadCounter",
                  paste("-o "    , outdir, "/", id, "_allele_counts.tsv", sep = ""),
                  paste("-I "    , outdir, "/", id, "_reordered.bam"    , sep = ""),
                  paste("-sites ", outdir, "/", id, "_hets.vcf"         , sep = ""),
                  "-overlap COUNT_FRAGMENTS_REQUIRE_SAME_BASE",
                  "-U ALLOW_N_CIGAR_READS",
                  sep = " \\\n"
                 )
    
    rm_uselsess = paste("rm -r", 
                        paste(outdir, "/", id, "_uniq.to.remap.bam"   , sep = ""),
                        paste(outdir, "/", id, "_uniq.to.remap.num.gz", sep = ""),
                        paste(outdir, "/", id, ".bam"                 , sep = ""),
                        paste(outdir, "/", id, "_filtered.bam"        , sep = ""),
                        paste(outdir, "/", id, "_reordered.bam"       , sep = ""),
                        paste(outdir, "/", id, "_reordered.bam.bai"   , sep = ""),
                        sep = " \\\n"
                       )
    
    out = paste("# waspAlignmentCompare",
                python1,
                java1,
                move1,
                java2,
                samtools_index,
                java3,
                #rm_uselsess,
                sep = "\n\n"
               )
    
    return(out);
}

### From Chris: MBASED
mbasedAnalysis = function(id, wgs_id, folder, bed, ppn,
                          mapability = "/publicdata/mapability_20151104/wgEncodeCrgMapabilityAlign100mer.bigWig",
                          chrom_conv = "/repos/cardips-pipelines/RNA/chrom_conv.tsv",
                          python     = "/frazer01/home/matteo/anaconda2/bin/python",
                          rscript    = "/frazer01/home/matteo/software/R-3.2.2/bin/Rscript"
                         )
{
    mbaseddir = paste(folder, "mbased", sep = "/")
    waspdir   = paste(folder, "wasp", sep = "/")
    
    create_outdir = paste(paste("mkdir -p", mbaseddir),
                          paste("cd"      , mbaseddir),
                          sep = "\n"
                         ) 
    
    python1 = paste(#paste(python   , "/frazer01/home/cdeboever/repos/cdeboever3/cdpipelines/cdpipelines/make_mbased_input.py"),
                    "make_mbased_input",
					paste(waspdir  , "/", id, "_allele_counts.tsv", sep = ""),
                    paste(mbaseddir, "/", id, "_mbased_input.tsv" , sep = ""),
                    bed,
                    paste("-v ", waspdir  , "/", id, ".input.vcf.gz", sep = ""),
                    paste("-s" , wgs_id    ),
                    paste("-m" , mapability),
                    "-p bigWigAverageOverBed",
                    paste("-c" , chrom_conv),
                    sep = " \\\n"
                   )
    
    rscript1 = paste(paste(rscript, "/frazer01/home/cdeboever/repos/cdeboever3/cdpipelines/cdpipelines/scripts/mbased.R"),
                     paste(mbaseddir, "/", id, "_mbased_input.tsv", sep = ""),
                     paste(mbaseddir, "/", id, "_locus.tsv"       , sep = ""),
                     paste(mbaseddir, "/", id, "_snv.tsv"         , sep = ""),
                     wgs_id,
                     "FALSE",
                     "100000",
                     ppn,
                     sep = " \\\n"
                    )
    
    out = paste("# MBASED", 
                create_outdir,
                python1,
                rscript1,
                sep = "\n\n"
               )
    
    return(out);
} 

### MBASED only
mbasedAnalysisOnly = function(id, wgs_id, folder, bed, ppn,
                              mapability = "/publicdata/mapability_20151104/wgEncodeCrgMapabilityAlign100mer.bigWig",
                              chrom_conv = "/repos/cardips-pipelines/RNA/chrom_conv.tsv",
                              python     = "/frazer01/home/matteo/anaconda2/bin/python",
                              rscript    = "/frazer01/home/matteo/software/R-3.2.2/bin/Rscript"
                             )
{
    mbaseddir = paste(folder, "mbased_genes", sep = "/")
    waspdir   = paste(folder, "wasp"        , sep = "/")
    
    create_outdir = paste(paste("mkdir -p", mbaseddir),
                          paste("cd"      , mbaseddir),
                          sep = "\n"
                         ) 
    
    python1 = paste(#paste(python   , "/frazer01/home/cdeboever/repos/cdeboever3/cdpipelines/cdpipelines/make_mbased_input.py"),
                    "make_mbased_input",
					paste(waspdir  , "/", id, "_allele_counts.tsv", sep = ""),
                    paste(mbaseddir, "/", id, "_mbased_input.tsv" , sep = ""),
                    bed,
                    paste("-v ", waspdir  , "/", id, ".input.vcf.gz", sep = ""),
                    paste("-s" , wgs_id    ),
                    paste("-m" , mapability),
                    "-p bigWigAverageOverBed",
                    paste("-c" , chrom_conv),
                    sep = " \\\n"
                   )
    
    rscript1 = paste(paste(rscript, "/frazer01/home/cdeboever/repos/cdeboever3/cdpipelines/cdpipelines/scripts/mbased.R"),
                     paste(mbaseddir, "/", id, "_mbased_input.tsv", sep = ""),
                     paste(mbaseddir, "/", id, "_locus.tsv"       , sep = ""),
                     paste(mbaseddir, "/", id, "_snv.tsv"         , sep = ""),
                     wgs_id,
                     "FALSE",
                     "100000",
                     ppn,
                     sep = " \\\n"
                    )
    
    out = paste("# MBASED", 
                create_outdir,
                python1,
                rscript1,
                sep = "\n\n"
               )
    
    return(out);
} 

### MBASED only
runMbased =  function (full_name, sh_folder, id, log_folder, wgs_id, folder, bam_files, bam, vcf, bed, ppn,
				       is.gz      = FALSE, # if false, gzip
				       run        = FALSE, # if true, run qsub
				       samtools   = "/software/samtools-1.2/samtools", 
                       genome_fai = "/publicdata/gatk_bundle_2.8/hg19/ucsc.hg19.fasta.fai", 
                       genome     = "/publicdata/gatk_bundle_2.8/hg19/ucsc.hg19.fasta", 
                       bcftools   = "/software/bcftools-1.2/bcftools", 
                       chrom_conv = "/repos/cardips-pipelines/RNA/chrom_conv.tsv",
                       python     = "/frazer01/home/matteo/anaconda2/bin/python",
				       genomeDir  = "/publicdata/star_index_hg19_sorted_20151123",
                       java       = "/usr/lib/jvm/java-1.7.0/bin/java",
				       mapability = "/publicdata/mapability_20151104/wgEncodeCrgMapabilityAlign100mer.bigWig",
                       rscript    = "/frazer01/home/matteo/software/R-3.2.2/bin/Rscript"
				      )
{
	sh_file = paste(writeHeader         (full_name, folder, log_folder, ppn),
                    mbasedAnalysisOnly  (id, wgs_id, folder, bed, ppn),
                    sep = "\n\n"
                   )

	output_file = paste(sh_folder, "/job.", full_name, ".sh", sep = "")
	write (sh_file, file = output_file)
	
	if (run == TRUE)
	{
		system(paste("qsub", output_file), intern = FALSE)
	}
	
	return(output_file)
}  

### Run ASE
### Calls: writeHeader, waspAlleleSwap, waspRemap, waspAlignmentCompare, mbasedAnalysis
runAse = function (full_name, sh_folder, id, log_folder, wgs_id, folder, bam_files, bam, vcf, bed,que, ppn,
				   is.gz      = FALSE, # if false, gzip
				   run        = FALSE, # if true, run qsub
				   samtools   = "/software/samtools-1.2/samtools", 
                   genome_fai = "/publicdata/gatk_bundle_2.8/hg19/ucsc.hg19.fasta.fai", 
                   genome     = "/publicdata/gatk_bundle_2.8/hg19/ucsc.hg19.fasta", 
                   bcftools   = "/software/bcftools-1.2/bcftools", 
                   chrom_conv = "/repos/cardips-pipelines/RNA/chrom_conv.tsv",
                   python     = "/frazer01/home/matteo/anaconda2/bin/python",
				   genomeDir  = "/publicdata/star_index_hg19_sorted_20151123",
                   java       = "/usr/lib/jvm/java-1.7.0/bin/java",
				   mapability = "/publicdata/mapability_20151104/wgEncodeCrgMapabilityAlign100mer.bigWig",
                   rscript    = "/frazer01/home/matteo/software/R-3.2.2/bin/Rscript"
				  )
{
	sh_file = paste(writeHeader         (full_name, folder, log_folder, que, ppn),
					mergeBams           (bam_files, folder, bam, ppn),
                    waspAlleleSwap      (id, wgs_id, folder, bam, vcf, bed, is.gz = is.gz),
                    waspRemap           (id, folder, ppn),
                    waspAlignmentCompare(id, folder, ppn),
                    mbasedAnalysis      (id, wgs_id, folder, bed, ppn),
                    sep = "\n\n"
                   )

	output_file = paste(sh_folder, "/job.", full_name, ".sh", sep = "")
	write (sh_file, file = output_file)
	
	if (run == TRUE)
	{
		system(paste("qsub", output_file), intern = FALSE)
	}
	
	return(output_file)
}
