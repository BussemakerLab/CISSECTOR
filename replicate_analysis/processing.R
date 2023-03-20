chrN <- 17
chr <- readLines("/Users/vince/Documents/SuRE/Vince/newData/chrList.txt")[chrN]

suppressMessages(library(GenomicRanges))

rawdata_path <- "/Users/vince/Documents/SuRE2/rawdata/"
output_path <- "/Users/vince/Documents/SuRE2/processed/SuRE_42-45/"


files <- list.files(paste0(rawdata_path,"SuRE_42-45"), full.names = T, pattern = paste0(chr,"[.]"),recursive=T)
libraries <- gsub(".*_(SuRE[0-9_]{,4})/.*","\\1",files)
libraries[nchar(libraries) == 6] <- paste(libraries[nchar(libraries) == 6],"1",sep="_")
nlib <- length(files)

df <- do.call("rbind", lapply(1:length(files), function(i){
  
  x <- read.delim(files[i], header = T, stringsAsFactors = F, nrows = 1000000)[,c(2:4,10:15)]
  
  count_cols <- grep("SuRE",colnames(x))
  hep_cols <- grep("HEPG2",colnames(x)) 
  k_cols <- count_cols[!(count_cols %in% hep_cols)]   
  
  x <- x[,c(1:4,k_cols,hep_cols)]
  colnames(x) <- c("start","end","strand","iPCR","K562_1","K562_2","K562_3","HEPG2_1","HEPG2_2")
  
  x$lib <- rep(i,nrow(x))
  return(x)
}))

df <- df[df$end <= 5000000,]

df$strand <- factor(df$strand, levels = c("+","-"))
df <- df[order(df$strand,df$start),]

grp <- 1

grp_name <- paste(chr,grp,sep = "_")
grp_output_path <- paste0(output_path, grp_name)
dir.create(grp_output_path)

bin_annotation_file <- paste0(grp_output_path,"/bin_annotation.rds")
element_annotation_file <- paste0(grp_output_path,"/element_annotation.rds")



gr <- GRanges(chr, IRanges(df$start, df$end), strand = df$strand)
bins <- disjoin(gr)
gr_start <- gr
end(gr_start) <- start(gr)
df$firstBin <- subjectHits(findOverlaps(gr_start,bins))
df$nBins <- countOverlaps(gr,bins) #seq.int(from = 10, length.out = 6)
saveRDS(df, element_annotation_file)

bin_anno <- as.data.frame(bins)[,c(2,4,5)]
bin_anno$coverage <- countOverlaps(bins,gr)
saveRDS(bin_anno, bin_annotation_file)

