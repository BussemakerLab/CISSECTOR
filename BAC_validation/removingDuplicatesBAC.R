library(GenomicRanges)
datafile <- "Documents/SuRE/Vince/SuRE49_predictions.rds"
fragments <- readRDS(datafile)
red <- reduce(fragments,ignore.strand=T)


c1 <- read.csv("Documents/SuRE/Vince/BAC/chr1_6.csv",header = T)
c1 <- c1[c1$end > 109468616,]
c1 <- c1[c1$start < 109646847,]

c1r <- GRanges(seqnames = "chr1", ranges = IRanges(c1$start,c1$end), strand = c1$strand)

rm(c1)

c1 <- read.csv("Documents/SuRE/Vince/BAC/chr1_7.csv",header = T)
c1 <- c1[c1$end > 155087487,]
c1 <- c1[c1$start < 155235691,]

c1r <- c(c1r,GRanges(seqnames = "chr1", ranges = IRanges(c1$start,c1$end), strand = c1$strand))


c1 <- read.csv("Documents/SuRE/Vince/BAC/chr5_7.csv",header = T)
c1 <- c1[c1$end > 139961656,]
c1 <- c1[c1$start < 140166144,]

c1r <- c(c1r,GRanges(seqnames = "chr5", ranges = IRanges(c1$start,c1$end), strand = c1$strand))


c1 <- read.csv("Documents/SuRE/Vince/BAC/chr6_2.csv",header = T)
c1 <- c1[c1$end > 26117719,]
c1 <- c1[c1$start < 26260178,]

c1r <- c(c1r,GRanges(seqnames = "chr6", ranges = IRanges(c1$start,c1$end), strand = c1$strand))

rm(c1)

oa <- overlapsAny(fragments, c1r, type = "equal")
fragments <- fragments[!oa]
saveRDS(fragments, "Documents/SuRE/Vince/SuRE49_predictions_dup_removed.rds")
length(fragments)
#5802 removed