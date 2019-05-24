source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
library(dada2)
path <- "~/scratch/Biogeography_Paper/"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,190),
                     trimLeft=c(20,20),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, pool=TRUE, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, pool=TRUE, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
taxa <- assignTaxonomy(seqtab.nochim, "~/data/kjduncan/rdp_train_set_16.fa.gz", multithread=TRUE)
taxa.plus <- addSpecies(taxa, "~/data/kjduncan/rdp_species_assignment_16.fa.gz",verbose=TRUE)
colnames(taxa.plus)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
unname(head(taxa.plus))
save.image(file="Biogeography_Paper_DADA2.RData")
