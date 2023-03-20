library(rtracklayer)

homedir <- "/rigel/hblab/users/vdf2104/"

oldTrackFiles <- paste0(homedir, "git/SuRE/Vince/results/", c("plus","minus"), "Strand_poisson_GLM_K562_sure23_hg19_VDF_281117.bw")
newTrackFiles <- paste0(homedir, "processed/SuRE_42-45/K562_",c("fwd","rev"),"_coef_track.bw")
sure49File <- "/rigel/hblab/vega/hblab/projects/SuRE/rawdata/SuRE49_50_BACs_for_wholeGenomeGLM_validation_LP170420/LP170328_SuRE49/SuRE-pipelineOutput/SuRE-counts.txt.gz"
outFile <- paste0(homedir, "processed/SuRE_42-45/SuRE49_predictions.rds")
#read in sure49 data

#make ranges
gr <- makeGRangesFromDataFrame(read.table(sure49File,header=T), keep.extra.columns = T)
red <- reduce(gr,ignore.strand =T)
red <- red[order(width(red),decreasing = T)[1:4]]

gr <- subsetByOverlaps(gr, red)

gr$old <- 0
gr$new <- 0

isPlus <- as.logical(strand(gr) == "+")

gr$old[isPlus] <- sum(import.bw(oldTrackFiles[1], selection = BigWigSelection(gr[isPlus]), as = "NumericList"))
gr$old[!isPlus] <- sum(import.bw(oldTrackFiles[2], selection = BigWigSelection(gr[!isPlus]), as = "NumericList"))
gr$new[isPlus] <- sum(import.bw(newTrackFiles[1], selection = BigWigSelection(gr[isPlus]), as = "NumericList"))
gr$new[!isPlus] <- sum(import.bw(newTrackFiles[2], selection = BigWigSelection(gr[!isPlus]), as = "NumericList"))

saveRDS(gr, file = outFile)
