chrN <- as.numeric(commandArgs(trailingOnly=T)[[1]])
chr <- readLines("/rigel/hblab/users/vdf2104/git/SuRE/Vince/newData/chrList.txt")[chrN]

suppressMessages(library(GenomicRanges))

rawdata_path <- "/rigel/hblab/users/vdf2104/rawdata/"
output_path <- "/rigel/hblab/users/vdf2104/processed/SuRE_42-45/"
max_grp_length <- 8000000


files <- list.files(paste0(rawdata_path,"SuRE_42-45"), full.names = T, pattern = paste0(chr,"[.]"),recursive=T)
libraries <- gsub(".*_(SuRE[0-9_]{,4})/.*","\\1",files)
libraries[nchar(libraries) == 6] <- paste(libraries[nchar(libraries) == 6],"1",sep="_")
nlib <- length(files)

df <- do.call("rbind", lapply(1:length(files), function(i){
    
    x <- read.delim(files[i], header = T, stringsAsFactors = F)[,-c(1,5:7)]
    if(ncol(x) == 15) {
	x <- x[,-15]
	}
    x <- x[-grep("2|3",x$SNPvarInf),] #easier to ignore all SNPs not 0 or 1

    count_cols <- grep("SuRE",colnames(x))
    hep_cols <- grep("HEPG2",colnames(x)) 
    k_cols <- count_cols[!(count_cols %in% hep_cols)]   

    x$K562 <- rowSums(x[,k_cols])
    x$HEPG2 <- rowSums(x[,hep_cols])
    
    x <- x[,-count_cols]
    x$lib <- rep(i,nrow(x))
    return(x)
}))

df$strand <- factor(df$strand, levels = c("+","-"))
df <- df[order(df$strand,df$start),]

ranges <- IRanges(start = df$start, end = df$end)
red <- reduce(ranges)
df$grp <- ceiling(cumsum(width(red))/max_grp_length)[subjectHits(findOverlaps(ranges,red))]
grps <- unique(df$grp)
rm(ranges,red)


stats_colnames <- paste(rep(c("ref","alt"),each = nlib*2),rep(c("plus","minus"),each = nlib), libraries, sep = ".")

for(grp in grps){
    
    grp_name <- paste(chr,grp,sep = "_")
    grp_output_path <- paste0(output_path, grp_name)
    dir.create(grp_output_path)
    snp_anno_file <- paste0(grp_output_path,"/SNP_annotation.rds")
    snp_stats_file <- paste0(grp_output_path,"/SNP_stats.rds")
    snp_triplets_file <- paste0(grp_output_path,"/SNP_triplets.rds")
    bin_annotation_file <- paste0(grp_output_path,"/bin_annotation.rds")
    element_annotation_file <- paste0(grp_output_path,"/element_annotation.rds")
    
    f <- df$grp == grp
    
    df_grp <- df[f,]
    df <- df[!f,]
    rm(f)
    
    vars <- strsplit(df_grp$SNPvarInf, ",")
    snps <- data.frame(element = rep(1:nrow(df_grp), lengths(vars)),
                       idx = as.integer(unlist(strsplit(df_grp$SNPidx, ","))),
                       rsid = unlist(strsplit(df_grp$SNP_ID, ",")),
                       pos = as.numeric(unlist(strsplit(df_grp$SNPabspos, ","))),
                       var = unlist(vars),
                       base = unlist(strsplit(df_grp$SNPbaseInf, ",")))
    rm(vars)
    df_grp <- df_grp[,c(1:3,6,10:12)]
    
    #removing duplicate rows
    dups <- duplicated(snps[,1:2])
    snps <- snps[!dups,]
    
    #removing rows for SNPs with only one observed allele
    varTable <- table(snps$idx, snps$var)
    obsAlleles <- rowSums(varTable > 0)
    oneAllele <- which(obsAlleles == 1)
    snps <- snps[!snps$idx %in% rownames(varTable)[oneAllele],]
    varTable <- varTable[-oneAllele,]
    rm(dups, obsAlleles, oneAllele)
    

    
    snp_stats <- list(N = matrix(data = table(snps$idx, df_grp$lib[snps$element], df_grp$strand[snps$element], snps$var), nrow = nrow(varTable), ncol = 4*nlib, dimnames = list(NULL, stats_colnames)),
                      K562 = matrix(data = tapply(df_grp$K562[snps$element], list(snps$idx, df_grp$lib[snps$element], df_grp$strand[snps$element], snps$var), mean), nrow = nrow(varTable), ncol = 4*nlib, dimnames = list(NULL, stats_colnames)),
                      HEPG2 = matrix(data = tapply(df_grp$HEPG2[snps$element], list(snps$idx, df_grp$lib[snps$element], df_grp$strand[snps$element], snps$var), mean), nrow = nrow(varTable), ncol = 4*nlib, dimnames = list(NULL, stats_colnames)),
                      iPCR = matrix(data = tapply(df_grp$iPCR[snps$element], list(snps$idx, df_grp$lib[snps$element], df_grp$strand[snps$element], snps$var), mean), nrow = nrow(varTable), ncol = 4*nlib, dimnames = list(NULL, stats_colnames)))
    saveRDS(snp_stats, snp_stats_file)
    rm(snp_stats)
    
    #compiling SNP annotation
    snp_alleles <- snps[,-1]
    snp_alleles <- snp_alleles[!duplicated(snp_alleles),]
    snp_alleles <- snp_alleles[order(as.numeric(snp_alleles$idx), snp_alleles$var),]
    n <- seq(1,nrow(snp_alleles),2)
    
    snp_anno <- data.frame(snp_alleles[n,1:3], 
                           base.ref = snp_alleles$base[n], 
                           base.alt = snp_alleles$base[n+1],
                           n.ref = varTable[,1],
                           n.alt = varTable[,2])
    rownames(snp_anno) <- NULL
    saveRDS(snp_anno,snp_anno_file)
    
    new_snp_index <- rep(NA,max(snp_anno$idx)+1)
    new_snp_index[snp_anno$idx + 1] <- 1:nrow(snp_anno)
    
    snp_triplets <- data.frame(i = snps$element, j = new_snp_index[snps$idx + 1], var = snps$var)
    saveRDS(snp_triplets,snp_triplets_file)
    rm(snp_alleles,new_snp_index,snps,varTable,n, snp_triplets,snp_anno)
    
    gr <- GRanges(chr, IRanges(df_grp$start, df_grp$end), strand = df_grp$strand)
    bins <- disjoin(gr)
    gr_start <- gr
    end(gr_start) <- start(gr)
    df_grp$firstBin <- subjectHits(findOverlaps(gr_start,bins))
    df_grp$nBins <- countOverlaps(gr,bins) #seq.int(from = 10, length.out = 6)
    saveRDS(df_grp, element_annotation_file)
    rm(df_grp, gr_start)
    
    bin_anno <- as.data.frame(bins)[,c(2,4,5)]
    bin_anno$coverage <- countOverlaps(bins,gr)
    saveRDS(bin_anno, bin_annotation_file)
    rm(bins,bin_anno,gr)
    
}

cat("\n",chr,"complete")

