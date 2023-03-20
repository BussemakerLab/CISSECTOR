args <- commandArgs(trailingOnly=T)
cell <- args[[1]]

library(GenomicRanges)
library(rtracklayer)

dir_path <- "/rigel/hblab/users/vdf2104/processed/SuRE_42-45"

pos_track_file <- paste0(dir_path,"/",cell,"_fwd_coef_track.bw")
neg_track_file <- paste0(dir_path,"/",cell,"_rev_coef_track.bw")

#get all the grp names
grp_names <- dir(dir_path, full.names = F,recursive = F, pattern = "^chr")

#get chrs
chrs <- sapply(strsplit(grp_names,"_"),"[[",1)
chr_levels <- unique(chrs)

#construct the filenames
anno_files <- paste(dir_path, grp_names, "bin_annotation.rds", sep = "/")
coef_files <- paste(dir_path, grp_names, paste0("bin_coefs_",cell,"_spatial.rds"), sep = "/")


#initiate lists
pos_track <- vector("list", length = length(chr_levels))
names(pos_track) <- chr_levels
neg_track <- pos_track

for(chr in chr_levels){
  n <- which(chrs == chr)
  df <- do.call("rbind", lapply(n, function(i){
    anno <- readRDS(anno_files[i])
    coefs <- readRDS(coef_files[i])
    cbind(anno,coefs)
  }))
  is_pos <- df$strand == "+"
  df_pos <- df[is_pos,]
  df_neg <- df[!is_pos,]
  rm(df)
  
  pos_track[[chr]] <- coverage(IRanges(start = df_pos$start, width = df_pos$width), weight = df_pos$coefs/sqrt(df_pos$width))

  
  neg_track[[chr]] <- coverage(IRanges(start = df_neg$start, width = df_neg$width), weight = df_neg$coefs/sqrt(df_neg$width))

  
  rm(df_pos,df_neg,is_pos)
}


pos_track <- RleList(pos_track, compress = F)
names(pos_track) <- chr_levels
export(pos_track, pos_track_file, format = "bw")

neg_track <- RleList(neg_track, compress = F)
names(neg_track) <- chr_levels
export(neg_track, neg_track_file, format = "bw")
