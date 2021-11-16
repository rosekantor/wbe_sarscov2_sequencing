#!/usr/bin/env Rscript

# This script follows the pipeline laid out in the dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
## running this script:
## Rscript dada2_pipeline.R amplicon indir outdir
### amplicon is a string (e.g. RBD)
### indir and outdir are full paths
### NOTE: you must provide reads that have already been trimmed (e.g. by cutadapt) to remove adapters and primers
### trimmed reads must be in the indir path provided

## outputs
### fasta
#### dada2_out.fasta # contains the unique Amplicon Sequence Variants (ASVs) produced by the pipeline

### tables
#### dada2_out.csv # columns are ASVs, rows are samples, fill is read count for each ASV in each sample, ASVs are ordered as in the fasta
#### read_tracking.csv # rows are samples, columns are steps in the pipeline, fill is read count for each sample x step, use to see where reads were lost

### Plots
#### dada2_read_quality.pdf # shows the read quality for the first pair of read files
#### dada2_error_rates.pdf # shows the error models made by dada2 using your data
# The error rates for each possible transition (A→C, A→G, …) are shown.
# Points are the observed error rates for each consensus quality score.
# The black line shows the estimated error rates after convergence of the machine-learning algorithm.
# The red line shows the error rates expected under the nominal definition of the Q-score.
# Check that the estimated error rates (black line) are a good fit to the observed rates (points),
# and that the error rates drop with increased quality as expected.


library('ggplot2')
library('dada2')

# take options for command line arguments
args = commandArgs(trailingOnly=TRUE)
amplicon <- args[1]
path <- args[2]
outdir <- args[3]

# get lists of forward and reverse reads
fnFs <- sort(list.files(path, pattern=".1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)

print("processing samples:")
print(sample.names)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(275,275),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=FALSE,
              compress=TRUE, multithread=TRUE)
colnames(out) <- c('input', 'filtered')
row.names(out) <- sample.names

# some samples were completely removed by filtering, so
# we need to remake the named lists filtFs and filtRs to include only the files that exist
filtered_path <- paste(path, "filtered", sep='/')

filtFs <- sort(list.files(filtered_path, pattern="F_filt.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(filtered_path, pattern="R_filt.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(filtFs), "_F_filt.fastq.gz"), `[`, 1)

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# merge forward and reverse reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

# remove sequences that are too short or too long
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 450:570]

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab_dim <- dim(seqtab.nochim)
print(paste('final sample count is',seqtab_dim[1]))
print(paste('final ASV count is',seqtab_dim[2]))

# show what fraction of reads remain after removing bimeras (should be most of them!)
frac_nonchimeric = sum(seqtab.nochim)/sum(seqtab)
print(paste('fraction non-chimeric reads was', frac_nonchimeric))

# save final table
outfile <- paste(toupper(amplicon), "dada2_out.csv", sep="_")
output <- paste(outdir, outfile, sep="/")
write.table(seqtab.nochim, output, sep=',')

# save fasta of sequences
outfile <- paste(toupper(amplicon), "dada2_out.fasta", sep="_")
output <- paste(outdir, outfile, sep="/")
uniquesToFasta(getUniques(seqtab.nochim), fout=output, ids=paste0("Seq", seq(length(getUniques(seqtab.nochim)))))

# save file tracking reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
## merge with the output from filtering (counts before and after filtering reads)
tracking <- merge(out, track, by.x=0, by.y=0, all=TRUE)
outfile <- paste(toupper(amplicon), "read_tracking.csv", sep="_")
output <- paste(outdir, outfile, sep="/")
write.table(tracking, output, sep=',')

# save plots

plt_errors = plotErrors(errF, nominalQ=TRUE)
outfile <- paste(toupper(amplicon), "dada2_error_rates.pdf", sep="_")
output <- paste(outdir, outfile, sep="/")
ggsave(filename=output, plot=plt_errors)

# this plot sometimes breaks, depending on whether there was a sample with no reads
# therefore putting it last
plt_readQual = plotQualityProfile(fnFs[1:2])
outfile <- paste(toupper(amplicon), "dada2_read_quality.pdf", sep="_")
output <- paste(outdir, outfile, sep="/")
ggsave(filename=output, plot=plt_readQual)
