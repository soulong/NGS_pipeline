
library(profileplyr)

setwd('Z:\\Data1\\NGS\\2025-05-20_H358_ATACseq\\results')

proplyrObject <- import_deepToolsMat(example) 


signalfiles <- list.files('bowtie2\\merged_replicate\\bigwig', '.bigWig$', full.names = T)
signalfiles <- list.files('bowtie2\\merged_replicate', '.bam$', full.names = T)
consensus <- rtracklayer::import('bowtie2\\merged_replicate\\macs2\\narrow_peak\\consensus\\consensus_peaks.mRp.clN.bed')

chipProfile <- BamBigwig_to_chipProfile(signalfiles, 
                                        consensus, 
                                        format = "bam",
                                        style="percentOfRegion",
                                        nOfWindows=100,
                                        distanceAround=100)
proplyrObject <- as_profileplyr(chipProfile)

