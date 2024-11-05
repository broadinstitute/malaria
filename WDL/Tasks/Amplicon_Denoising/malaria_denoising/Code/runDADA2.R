#!/bin/R env

########################################
# LOAD LIBRARIES AND PREPARE VARIABLES #
########################################

print("Entering DADA2 script")

if (!require("dada2")) {
  install.packages("dada2", repos="http://cran.rstudio.com/")
  library("dada2")
}

if (!require("limma")) {
  install.packages("limma", repos="http://cran.rstudio.com/")
  library("limma")
}

if (!require("data.table")) {
  install.packages("data.table", repos="http://cran.rstudio.com/")
  library("data.table")
}

if (!require("argparse")) {
  install.packages("argparse", repos="http://cran.rstudio.com/")
  library("argparse")
}

if (!require("viridisLite")) {
  install.packages("viridis", repos="http://cran.rstudio.com/")
  library("viridis")
}

#qprofile <- function(fastq, work_dir) {
#  # Plot Quality profiles before filtering
#  fastq_name = sub("\\.fq\\.gz$", "", basename(fastq))
#  pdf(paste0(work_dir,"/QProfile/", fastq_name, ".pdf"))
#  try(print(plotQualityProfile(fastq)), silent = TRUE)
#  dev.off()
#}

# Custom filtering, denoising parameters (if not default) can be provided as a separate config file?
parser <- ArgumentParser()
parser$add_argument("-b", "--sample_ids", nargs="+", help="Samples names")
parser$add_argument("-f1", "--orig_f", nargs="+", help="Forward fastq files before any processing (required)")
parser$add_argument("-f2", "--orig_r", nargs="+", help="Reverse fastq files before any processing (required)")
parser$add_argument("-p1", "--fnFs", nargs="+", help="Forward fastq files after primer removal (required)")
parser$add_argument("-p2", "--fnRs", nargs="+", help="Reverse fastq files after primer removal (required)")
parser$add_argument("-a1", "--adap_f", nargs="+", help="Forward fastq files after adaptor removal (required)")
parser$add_argument("-a2", "--adap_r", nargs="+", help="Reverse fastq files after adaptor removal (required)")
parser$add_argument("-c", "--class", help="Class specifying 'parasite' or 'vector' (required if '--default' is specified)")
parser$add_argument("-d", "--work_dir", help="Working directory path for writing all dada2 output files")
parser$add_argument("-o", "--output_filename", help="output tab-separated filename (required)")
parser$add_argument("-s", "--save_run", help="save Run as R workspace image")
parser$add_argument("-ee", "--maxEE", help="Maximum expected errors for filtering forward and reverse read")
parser$add_argument("-tR", "--trimRight", help="Length for trimming from right for both forward and reverse reads")
parser$add_argument("-mL", "--minLen", type="integer", help="Minimum length required for reads on both end. Shorter reads are discarded")
parser$add_argument("-tQ", "--truncQ", help="Truncate reads to first occurence of truncQ. All filtered reads have quality >= truncQ")
parser$add_argument("-id", "--matchIDs", type="integer", help="Match ids on fastqs to make sure reads on forward and reverse end are in same order")
parser$add_argument("-mC", "--max_consist", type="integer", help="Maximum cycles for error model until consistency. If no convergence, error values at max_consist cycle are used")
parser$add_argument("-wA", "--omega_a", type="double", help="P-value threshold in sample inference for forming a new partition")
parser$add_argument("-jC", "--justConcatenate", type="integer", help="Specify whether ASVs need to be concatenated with Ns instead of merging")
parser$add_argument("-mM", "--maxMismatch", type="integer", help="Specify the maximum number of mismatches allowed during merging")
parser$add_argument("--bimera", action='store_true', help="Optionally output list of sequences identified as bimeras")
args <- parser$parse_args()

########################################
# 	PARSE PARAMETERS               #
########################################

# Universal parameters
sample_ids = unlist(strsplit(args$sample_ids, " ")) 
work_dir = args$work_dir
output_filename = args$output_filename
save_run = args$save_run
randomize = TRUE
selfConsist = TRUE
filter = TRUE
matchIDs <- args$matchIDs
class = args$class

# Parameters for merging
justConcatenate <- args$justConcatenate
bimera = args$bimera
maxMismatch = args$maxMismatch

# DADA2 and Filtering parameters
maxEE <- args$maxEE
trimRight <- args$trimRight
minLen <- args$minLen
truncQ <- args$truncQ
max_consist <- args$max_consist
omega_a <- args$omega_a

fnFs <- unlist(strsplit(args$fnFs, " "))
fnRs <- unlist(strsplit(args$fnRs, " "))

if (is.null(matchIDs)||matchIDs == '') {
  matchIDs = TRUE
} else {
  matchIDs = as.logical(as.numeric(matchIDs))
}
if (is.null(trimRight)||trimRight == '') {
  trimRight = c(0,0)
} else {
  trimRight <- as.numeric(strsplit(trimRight,',')[[1]])
}
if (is.null(truncQ)||truncQ == '') {
  truncQ = c(5,5)
} else {
  truncQ <- as.numeric(strsplit(truncQ,',')[[1]])
}

if (class == "parasite") {
  # Parameters for filtering
  if (is.null(maxEE)||maxEE == '') {
    maxEE = c(5,5)
  } else {
    maxEE <- as.numeric(strsplit(maxEE,',')[[1]])
  }
  if (is.null(minLen)||minLen == '') {
    minLen=30
  } else {
    minLen = as.numeric(minLen)
  }

  # Parameters for Denoising
  if (is.null(max_consist)||max_consist == '') {
    max_consist=10
  } else {
    max_consist = as.numeric(max_consist)
  }
  if (is.null(omega_a)||omega_a == '') {
    omega_a=1e-120
  } else {
    omega_a = as.numeric(omega_a)
  }

  #Parameters for merging
  if (is.null(justConcatenate)||justConcatenate == '') {
    justConcatenate = FALSE
  } else {
    justConcatenate = as.logical(as.numeric(justConcatenate))
  }
} else if (class == "vector") {
  # Parameters for filtering
  if (is.null(maxEE)||maxEE == '') {
    maxEE = c(2,2)
  } else {
    maxEE <- as.numeric(strsplit(maxEE,',')[[1]])
  }
  if (is.null(minLen)||minLen == '') {
    minLen=75
  } else {
    minLen = as.numeric(minLen)
  }

  # Parameters for Denoising
  if (is.null(max_consist)||max_consist == '') {
    max_consist=20
  } else {
    max_consist = as.numeric(max_consist)
  }
  if (is.null(omega_a)||omega_a == '') {
    omega_a=1e-40
  } else {
    omega_a = as.numeric(omega_a)
  }

  #Parameters for merging
  if (is.null(justConcatenate)||justConcatenate == '') {
    justConcatenate = TRUE
  } else {
    justConcatenate = as.logical(as.numeric(justConcatenate))
  }
  
  if (is.null(maxMismatch)||maxMismatch == '') {
    maxMismatch=30
  } else {
    maxMismatch = as.numeric(maxMismatch)
  }
  
} else {
  stop("Please provide valid option for the '--class' argument")
}

if (dirname(output_filename) != ".") {
  output_filename <- output_filename
  } else {
  output_filename <- paste0(work_dir, "/", output_filename)
}

#Datatable to summarize parmeters
parameter_df <- data.frame(maxEE=maxEE,
	trimRight=trimRight,
	minLen=minLen,
	truncQ=truncQ,
	matchIDs=matchIDs,
	max_consist=max_consist,
	randomize=randomize,
	selfConsist=selfConsist,
	OMEGA_A=omega_a,
	justConcatenate=justConcatenate)

print(parameter_df)

########################################
#                DADA2                 #
########################################

# List files and sample names
if (length(fnFs) == 0 || length(fnFs) != length(fnRs)) {
	stop("fastq files incomplete or not found")
}

#TO DELETE? THESE PLOTS ARE CURRENTLY NOT BEING USED. DISCUSS WITH UCSF TEAM
# Plot Quality profiles before filtering
#lapply(fnFs, qprofile, work_dir)
#lapply(fnRs, qprofile, work_dir)

# Create paths for filtered fastq
filtFs <- file.path(work_dir, "filtered", paste0(sample_ids, "_filt_R1.fastq.gz"))
filtRs <- file.path(work_dir, "filtered", paste0(sample_ids, "_filt_R2.fastq.gz"))

names(filtFs) <- sample_ids
names(filtRs) <- sample_ids

## Filter read
if (filter == TRUE) {
	print("filtering samples...")
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
	maxN=0, maxEE=maxEE, trimRight=trimRight, truncQ=truncQ, minLen=minLen,
	rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose=TRUE,
	matchIDs=matchIDs)
	print("filtering done!")
} else {
	print("skipping filter except mandatory removal of N's... ")
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=c(0,0), maxN=0, rm.phix=TRUE,
	compress=TRUE, multithread=TRUE, verbose=TRUE, matchIDs=matchIDs)
}

# Report and Correct for samples with zero reads after filter
zeros <- row.names(out)[out[,2] == 0]
zeros = sub("^([^_]+_[^_]+)_.*", "\\1", zeros)
write.table(zeros, paste0(work_dir, "/zeroReadSamples.txt"), sep = "\t", quote = FALSE)
filtFs <- filtFs[out[,2] != 0]
filtRs <- filtRs[out[,2] != 0]
sample_ids <- sample_ids[out[,2] != 0]

# Update Out table
out <- out[(out[,2] != 0),]

#Compute the error model
print("starting error model learning for forward reads...")
errF <- learnErrors(filtFs, multithread=TRUE, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)
print("starting error model learning for reverse reads...")
errR <- learnErrors(filtRs, multithread=TRUE, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)

#TO DELETE? THESE PLOTS ARE CURRENTLY NOT BEING USED. DISCUSS WITH UCSF TEAM
#Plot the Errors
#pdf(paste0(work_dir,"/errF.pdf"))
#try(print(plotErrors(errF, nominalQ=TRUE)), silent = TRUE)
#dev.off()
#pdf(paste0(work_dir,"/errR.pdf"))
#try(print(plotErrors(errR, nominalQ=TRUE)), silent = TRUE)
#dev.off()

#DeReplicate the reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample_ids
names(derepRs) <- sample_ids

#Run core DADA2 algorithm
print("starting dada2 for forward reads...")
dadaFs <- dada(derepFs, err=errF, selfConsist=selfConsist, multithread=TRUE, verbose=TRUE, OMEGA_A=omega_a)
print("starting dada2 for reverse reads...")
dadaRs <- dada(derepRs, err=errR, selfConsist=selfConsist, multithread=TRUE, verbose=TRUE, OMEGA_A=omega_a)

# Merge reads
print("merging paired ends...")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=justConcatenate, trimOverhang = TRUE, maxMismatch = maxMismatch)

#Generate sequence table
print("generating sequence table...")
seqtab <- makeSequenceTable(mergers)
print("Number of sequences in table")
print(dim(seqtab))
# Inspect distribution of sequence lengths
print(table(nchar(getSequences(seqtab))))

#Remove Chimeras
if(bimera) {
  print("identifying bimeric sequences...")
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  print("Number of non-bimeric sequences:")
  print(dim(seqtab.nochim)[2])
  print("Percentage of reads which are non-bimeric:")
  print(sum(seqtab.nochim)/sum(seqtab))
  bimeras <- !(colnames(seqtab) %in% colnames(seqtab.nochim))
  write.table(data.frame(sequence = colnames(seqtab), bimera = bimeras), file=paste0(work_dir,"/ASVBimeras.txt"),
    quote=FALSE, sep="\t", row.names=FALSE)
} else {
  print("skipping Bimera identification..")
}

# Track reads through the pipeline
print("tracking reads through the pipeline...")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_ids

#Get reads from original files
reads_report = data.frame()
count_reads <- function(file_path) {

  gz_conn <- gzfile(file_path, "r")  
  read_count <- 0
  
  # Iterate over the lines in the file
  while (length(line <- readLines(gz_conn, n = 4)) > 0) {
    if (length(line) == 4 && substr(line[1], 1, 1) == "@") {
      read_count <- read_count + 1
    }
  }
  
  close(gz_conn)  
  return(read_count)
}

reads_report <- data.frame()
for(sample.name in names(track[,1])) {
  orig_f <- unlist(strsplit(args$orig_f, " "))[grep(sample.name, unlist(strsplit(args$orig_f, " ")))]
  orig_r <- unlist(strsplit(args$orig_r, " "))[grep(sample.name, unlist(strsplit(args$orig_r, " ")))]
  orig_f_reads = count_reads(orig_f)
  orig_r_reads = count_reads(orig_r)

  adap_f <- unlist(strsplit(args$adap_f, " "))[grep(sample.name, unlist(strsplit(args$adap_f, " ")))]
  adap_r <- unlist(strsplit(args$adap_r, " "))[grep(sample.name, unlist(strsplit(args$adap_r, " ")))]  
  adap_f_reads = count_reads(adap_f)
  adap_r_reads = count_reads(adap_r)

  prim_f <- unlist(strsplit(args$fnFs, " "))[grep(sample.name, unlist(strsplit(args$fnFs, " ")))]
  prim_r <- unlist(strsplit(args$fnRs, " "))[grep(sample.name, unlist(strsplit(args$fnRs, " ")))]  
  prim_f_reads = count_reads(prim_f)
  prim_r_reads = count_reads(prim_r)
  
  reads_report = rbind(reads_report,
                       data.frame(sample.name,
                                  orig_f_reads,
                                  orig_r_reads,
                                  adap_f_reads,
                                  adap_r_reads,
                                  prim_f_reads,
                                  prim_r_reads))
}

stopifnot(names(track[,1]) == reads_report$sample_ids)
rownames(reads_report) = reads_report$sample_ids
reads_report = data.matrix(reads_report)

track = cbind(track, data.matrix(reads_report))
original_discarded = track[,"orig_f_reads"] - track[,"adap_f_reads"] 
adaptor_discarded = track[,"adap_f_reads"] - track[,"prim_f_reads"] 
primer_discarded = track[,"prim_f_reads"] - track[,"filtered"]
filtered_discarded = track[,"filtered"] - track[,"denoisedF"]
merged_discarded = track[,"denoisedF"] - track[,"merged"]
track = cbind(track, original_discarded, adaptor_discarded, primer_discarded, filtered_discarded, merged_discarded)

# sink summary from stdout to a file
sink(paste0(work_dir,"/reads_summary.txt"))
print(track)
sink()

#TO DELETE? THESE PLOTS ARE CURRENTLY NOT BEING USED. DISCUSS WITH UCSF TEAM
#Subset the table to the desired experiments
#samples_order = read.csv(path_to_flist, sep = ",", header = FALSE)$V1
#track_plot = as.data.frame(track[, c("merged", "merged_discarded", "filtered_discarded", "primer_discarded", "adaptor_discarded", "original_discarded")])
#track_plot = track_plot[row.names(track_plot) %in% samples_order,] 
#
#track_plot <- track_plot[order(-track_plot[,1], 
#                               -track_plot[,2], 
#                               -track_plot[,3],
#                               -track_plot[,4],
#                               -track_plot[,5],
#                               -track_plot[,6]
#                               ),]
#
#track_plot_per = sweep(track_plot, 1, rowSums(track_plot), "/")
#
#track_plot_per <- track_plot_per[order(-track_plot_per[,1], 
#                               -track_plot_per[,2], 
#                               -track_plot_per[,3],
#                               -track_plot_per[,4],
#                               -track_plot_per[,5],
#                               -track_plot_per[,6]
#),]
#
#color_vector = viridis(nrow(t(as.matrix(track_plot))), option = "D")
#
#pdf(paste0(work_dir,"/stacked_barplot.pdf"), width = 12)
#par(mar = c(14, 4, 6, 2), xpd = TRUE) # increase the bottom margin
#barplot(
#  t(as.matrix(track_plot)),
#  col = c(color_vector, "red"),
#  border = NA,
#  ylab = "Reads Count",
#  xlab = "",
#  main = "DADA2 performance - Absolute",
#  las = 2,
#  cex.names = 0.5,
#  cex.axis = 0.7,
#  cex.main = 2
#)
#
#mtext("Experiment ID", side = 1, line = 12) # lower x-axis label
#
#legend(
#  "top",
#  inset=c(0,-0.1),
#  horiz = TRUE,
#  legend = c("Merged", "Merged discarded", "Filtered discarded", "Primer discarded", "Adaptor discarded", "Original discarded"),
#  fill = c(color_vector, "red"),
#  bty = "n",
#  cex = 1
#)
#dev.off()
#
#pdf(paste0(work_dir,"/stacked_barplot_per.pdf"), width = 12)
#par(mar = c(14, 4, 6, 2), xpd = TRUE) # increase the bottom margin
#barplot(
#  t(as.matrix(track_plot_per)),
#  col = c(color_vector, "red"),
#  border = NA,
#  ylab = "Reads Percentage",
#  xlab = "",
#  main = "DADA2 performance - Percentage",
#  las = 2,
#  cex.names = 0.5,
#  cex.axis = 0.7,
#  cex.main = 2
#)
#
#mtext("Experiment ID", side = 1, line = 12) # lower x-axis label
#
#legend(
#  "top",
#  inset=c(0,-0.1),
#  horiz = TRUE,
#  legend = c("Merged", "Merged discarded", "Filtered discarded", "Primer discarded", "Adaptor discarded", "Original discarded"),
#  fill = c(color_vector, "red"),
#  bty = "n",
#  cex = 1
#)
#
#dev.off()

#Show the barplot of length distribution
pdf(paste0(work_dir,"/sequences_barplot.pdf"))
print(barplot(table(nchar(getSequences(seqtab)))))
dev.off()

#Add Zero Read Samples
zeros.df = data.frame(matrix(ncol = ncol(seqtab), nrow = length(zeros)))
colnames(zeros.df) = colnames(seqtab)
rownames(zeros.df) = zeros
seqtab = rbind(seqtab, zeros.df)

# Create a new column for Sample_ID using row names and move the last column to the first position
seqtab$Sample_ID <- row.names(seqtab)
row.names(seqtab) = NULL
col_names <- names(seqtab)
seqtab <- seqtab[, c(length(col_names), 1:(length(col_names)-1))]
seqtab[is.na(seqtab)] = 0

#Generate output: sequence table to a tsv
write.table(seqtab, file=output_filename, quote = FALSE, sep = "\t", row.names = FALSE)

# Save Run as R workspace image (Optional)
if (is.null(save_run)||save_run == '') {
    print("--save_run not found or empty. skip saving Rdata image to a file")
  } else {
    save.image(paste0(work_dir,"/",save_run))
}

print("Leaving DADA2 script")
