#!/usr/bin/env Rscript
library('dada2')

# hardcoded script to merge all sequence tables for each amplicon across sequencing batches

amplicons <- c('NTD', 'RBD', 'S1S2')
# for each amplicon, load all files using glob to collect files from all dates
for (amplicon in amplicons) {
    files <- Sys.glob(paste0('/Users/rosekantor/data/wbe_scv/qb3_sgene_*/results/dada2_out/', amplicon, '_dada2_out.csv'), dirmark = FALSE)

    # make an empty vector of the correct length to hold all the seqtabs as they are read in
    seqtabs <- vector("list", length(files))

    # read in the seqtabs
    for (i in seq(1,length(files))) {
        seqtab <- as.matrix(read.table(files[i], row.names=1, header = 1, sep=','))
        seqtabs[[i]] <- seqtab
    }

    # merge the seqtabs into one table
    merged <- mergeSequenceTables(tables=seqtabs)
    # collapse the table so that counts from sequences that completely contain other identical sequences get combined
    collapsed <- collapseNoMismatch(merged)

    # save the final seqtab and the fasta
    write.table(collapsed, paste0('/Users/rosekantor/data/wbe_scv/results_sgene/', amplicon, '_combined_dada2_out.csv'), sep=',')
    uniquesToFasta(getUniques(collapsed), fout=paste0('/Users/rosekantor/data/wbe_scv/results_sgene/', amplicon, '_combined_dada2_out.fasta'),
               ids=paste0("Seq", seq(length(getUniques(collapsed)))))
}
