library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(mygene)
library(scales)

datafile <- "Documents/SuRE/Vince/SuRE49_predictions_dup_removed.rds"
#datafile <- "Documents/SuRE/Vince/SuRE49_predictions.rds"
fragments <- readRDS(datafile)
fragments

genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
tss <- resize(genes, width = 1, fix = "start")
tss <- resize(tss, width = 2000, fix = "center")

tss$co <- countOverlaps(tss, fragments)
tss <- tss[order(tss$co, decreasing = T)]
tss <- tss[tss$co > 0]

fo <- findOverlaps(tss, fragments)
table(table(subjectHits(fo)))

tss$meanY <- tapply(fragments$SuRE23_49_B1_T1[subjectHits(fo)], queryHits(fo), sum)/tss$co
plot.ecdf(tss$meanY)
tss <- tss[order(tss$meanY, decreasing = T)]

genenames <- getGenes(tss$gene_id, fields=c("symbol", "name"))
tss$symbol <- genenames$symbol

rm(genenames)

ex <- grep("VTRNA|HIST", tss$symbol)


#binning
#find overlaps within tss regions
fo <- findOverlaps(fragments, tss, type = "within")
#create list of all fragments overlapping each tss
frag_by_tss <- tapply(fragments[queryHits(fo)], subjectHits(fo), c)
#make 100bp bins for each tss
bins_by_tss <- list()
for(i in 1:length(tss)){
  bins_by_tss[[i]] <- GRanges(seqnames = seqnames(tss[i]), ranges = IRanges(seq(start(tss[i]),end(tss[i]),50),width = 50), strand = strand(tss[i]))
}
# bins_by_tss <- lapply(tss, function(x){
#   GRanges(seqnames = seqnames(x), ranges = IRanges(seq(start(x),end(x),100),width = 100), strand = strand(x))
# })
#for each tss, find overlaps for start and end and combine
fo_start <- mapply(function(x,y) {
  subjectHits(findOverlaps(resize(x,1,"start"),y))
}, x = frag_by_tss, y = bins_by_tss)
fo_end <- mapply(function(x,y) {
  subjectHits(findOverlaps(resize(x,1,"end"),y))
}, x = frag_by_tss, y = bins_by_tss)
fo_bins <- mapply(function(x,y) paste(x,y,sep="_"), x = fo_start, y = fo_end)

#calculate count for each bin
bin_count <- lapply(fo_bins,table)

#for each bin in each tss calculate means
bin_means_obs <- mapply(function(x,y){
  tapply(x$SuRE23_49_B1_T1, y, mean)
}, x = frag_by_tss, y = fo_bins)
bin_means_obs2 <- mapply(function(x,y){
  tapply(x$SuRE23_49_B1_T1, y, sum)/tapply(x$iPCR, y, sum)
}, x = frag_by_tss, y = fo_bins)
bin_means_old <- mapply(function(x,y){
  tapply(exp(x$old), y, mean)
}, x = frag_by_tss, y = fo_bins)
bin_means_new <- mapply(function(x,y){
  tapply(exp(x$new), y, mean)
}, x = frag_by_tss, y = fo_bins)
bin_means_new2 <- mapply(function(x,y){
  tapply(x$new, y, mean)
}, x = frag_by_tss, y = fo_bins)


#fits
df <- as.data.frame(fragments)
fit_old <- glm(SuRE23_49_B1_T1 ~ old, family = poisson(), data = df)
fit_new <- glm(SuRE23_49_B1_T1 ~ new, family = poisson(), data = df)

#plots
bins_df <- data.frame(gene = rep(tss$symbol,lengths(bin_count)), n = unlist(bin_count), obs = unlist(bin_means_obs), old = unlist(bin_means_old), new = unlist(bin_means_new), new2 = unlist(bin_means_new2), obs2 = unlist(bin_means_obs2))

# par(mar = c(5,4,4,1)+.1, mfrow = c(1,1))
# xlim <- range(c(log2(bins_df$old), log2(bins_df$new)))
# cor(log2(bins_df$old), log2(bins_df$obs+1))
# plot(log2(bins_df$old), log2(bins_df$obs),xlab = "log2(mean SuRE23 predicted)", ylab = "log2(mean SuRE49 observed)", pch = 16, cex = sqrt(bins_df$n)/4, col = alpha("red",alpha = .2), xlim = xlim)
# abline(coef = coefficients(fit_old), col = "blue",lwd=1)
# text(0, 10, bquote(R^2 == .(round(cor(log2(bins_df$old), log2(bins_df$obs+1))^2, 2))))
# 
# cor(log2(bins_df$new), log2(bins_df$obs+1))
# plot(log2(bins_df$new), log2(bins_df$obs),xlab = "log2(mean SuRE42-45 predicted)", ylab = "log2(mean SuRE49 observed)", pch = 16, cex = sqrt(bins_df$n)/4, col = alpha("red",alpha = .2), xlim = xlim)
# abline(coef = coefficients(fit_new), col = "blue",lwd=1)
# text(0, 10, bquote(R^2 == .(round(cor(log2(bins_df$new), log2(bins_df$obs+1))^2, 2))))
# 
# par(pty="s")
# cor(bins_df$new2, log(bins_df$obs+1))
# xlim <- range(bins_df$new2)
# plot(bins_df$new2, log(bins_df$obs),xlab = "mean predicted log expression", ylab = "log(mean observed expression)", pch = 16, cex = sqrt(bins_df$n)/4, col = alpha("red",alpha = .2), xlim = xlim)
# abline(coef = coefficients(fit_new), col = "blue",lwd=1)
# text(0, 10, bquote(R^2 == .(round(cor(log2(bins_df$new), log2(bins_df$obs+1))^2, 2))))
# 
# #Figure 4.10 version
# library(RColorBrewer)
# n <- 50
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# 
# 
# plot(log2(bins_df$new), log2(bins_df$obs),xlab = "log2(mean predicted cDNA count)", ylab = "log2(mean observed cDNA count)", pch = 16, cex = sqrt(bins_df$n)/4, col = alpha(col_vector[as.integer(bins_df$gene)],alpha = .4), xlim = xlim)
# abline(coef = coefficients(fit_new), col = "blue",lwd=1)
# abline(coef = c(2.72,.7), col = "black",lwd=1)
# legend("topleft",bty = "n", legend = c(1,5,10,25,50,100), pch = 16, pt.cex = sqrt(c(1,5,10,25,50,100))/4, col = alpha("red",alpha = .2))
# 
# 
# 
# 
# ###
# 
# r2s_old <- tapply(1:nrow(bins_df), bins_df$gene, function(i){
#   cor(log2(bins_df$old[i]), log2(bins_df$obs[i]+1))^2
# })
# r2s_new <- tapply(1:nrow(bins_df), bins_df$gene, function(i){
#   cor(log2(bins_df$new[i]), log2(bins_df$obs[i]+1))^2
# })
# plot.ecdf(r2s_old, xlab = bquote(R^2))
# plot.ecdf(r2s_new, add=T,col="red")
# abline(v = c(cor(log2(bins_df$old), log2(bins_df$obs+1))^2, cor(log2(bins_df$new), log2(bins_df$obs+1))^2), col = c("black","red"))
# 
# 
# par(mfrow = c(5,5), mar = c(0,0,2,0))
# ylim = range(log2(bins_df$obs[bins_df$obs>0]))
# tapply(1:nrow(bins_df), bins_df$gene, function(i){
#   plot(log2(bins_df$old[i]), log2(bins_df$obs[i]), xlab = "", ylab = "", pch = 16, cex = sqrt(bins_df$n)/2, col = alpha("black",alpha = .5), xlim = xlim,ylim = ylim, main = bins_df$gene[i[1]], xaxt = "n", yaxt = "n")
#   
#   points(log2(bins_df$new[i]), log2(bins_df$obs[i]), pch = 16, cex = sqrt(bins_df$n)/2, col = alpha("red",alpha = .5))
#   abline(coef = coefficients(fit_new), col = "black",lwd=1)
# })
# 
# # Figure 4.11 version
# par(mfrow = c(5,5), mar = c(0,0,2,0))
# ylim = range(log2(bins_df$obs[bins_df$obs>0]))
# tapply(1:nrow(bins_df), bins_df$gene, function(i){
#   plot(log2(bins_df$new[i]), log2(bins_df$obs[i]), xlab = "", ylab = "", pch = 16, cex = sqrt(bins_df$n)/2, col = alpha("red",alpha = .5), xlim = xlim,ylim = ylim, main = bins_df$gene[i[1]], xaxt = "n", yaxt = "n")
#   abline(coef = coefficients(fit_new), col = "blue",lwd=1)
#   abline(coef = c(2.72,.7), col = "black",lwd=1)
# })
# 
# 
# 
# 
# par(mar = c(5,4,4,1)+.1, mfrow = c(1,1))
# plot(tss$meanY[order(tss$symbol)],r2s_old, xlab = "mean SuRE49 cDNA count", ylab = "SuRE23 vs. SuRE49 R squared")
# plot(tss$meanY[order(tss$symbol)],r2s_new, xlab = "mean SuRE49 cDNA count", ylab = "SuRE42-45 vs. SuRE49 R squared")
# 
# 

####new versions

#removing the redundant promoters
bins_df <- bins_df[!bins_df$gene %in% c("H1-4","H3C4","AKNAD1"),]

par(pty="s",mar = c(5,4,4,1)+.1, mfrow = c(1,1))
cor(bins_df$new2, log(bins_df$obs2+1))
cor.test(exp(bins_df$new2), bins_df$obs2)
cor.test(exp(bins_df$new2), bins_df$obs2)$p.value
xlim <- exp(c(-2.5,6.5))
f <- bins_df$obs2 > 0

plot(exp(bins_df$new2[f]), bins_df$obs2[f],xlab = "predicted expression", ylab = "mean observed expression", pch = 16, cex = sqrt(bins_df$n[f])/2, col = alpha("black",alpha = .2), xlim = xlim, ylim=xlim*(mean(df$SuRE23_49_B1_T1)/mean(df$iPCR))/mean(exp(df$new)), log="xy", axes = F, frame.plot = T)
axis(side = 1, at = c(0.1,1,10,100), labels = parse(text=paste("10^",-1:2, sep="")))
axis(side = 2, at = c(0.1,1,10,100), labels = parse(text=paste("10^",-1:2, sep="")))
curve(x*(mean(df$SuRE23_49_B1_T1)/mean(df$iPCR))/mean(exp(df$new)), from = 10^-2, to = 10^3, add=T, col = "blue")
legend(x = .1, y = 500, legend = c(1,5,10,20,40), pch = 16, pt.cex = sqrt(c(1,5,10,20,40))/2, col = alpha("black",alpha = .2))
text(1, 400, bquote(r == .(round(cor(exp(bins_df$new2), bins_df$obs2), 2))))
text(1, 300, bquote(p < 10^-16))



cors <- tapply(1:nrow(bins_df),as.character(bins_df$gene), function(n) {
  df <- bins_df[n,]
  df <- df[df$obs2 > 0,]
  return(cor(df$new2, log(df$obs2)))
})

cors2 <- tapply(1:nrow(bins_df),as.character(bins_df$gene), function(n) {
  df <- bins_df[n,]
  return(cor(exp(df$new2), df$obs2))
})

plot.ecdf(cors2, main = "", xlab = "Pearson correlation")



#plot(sapply(names(cors), function(x) tss$meanY[tss$symbol==x]),cors, pch = 16, col = "black", xlab = "mean observed expression", ylab = "correlation",log="x")

g <- c("HARS1", "H2BC8", "TMCO6","H3C6","VTRNA1-3","TMEM167B")
par(pty="s",mar = c(3,2,2,1), mfrow = c(3,2))
for(i in g){
  D <- bins_df[bins_df$gene == i,]
  fD <- D$obs2 > 0
  
  plot(exp(D$new2[fD]), D$obs2[fD],xlab = "", ylab = "", pch = 16, cex = sqrt(D$n[fD])/2, col = alpha("black",alpha = .5), xlim = xlim, ylim=xlim*(mean(df$SuRE23_49_B1_T1)/mean(df$iPCR))/mean(exp(df$new)), log="xy", axes = F, frame.plot = T)
  axis(side = 1, at = c(0.1,1,10,100), labels = parse(text=paste("10^",-1:2, sep="")))
  axis(side = 2, at = c(0.1,1,10,100), labels = parse(text=paste("10^",-1:2, sep="")))
  curve(x*(mean(df$SuRE23_49_B1_T1)/mean(df$iPCR))/mean(exp(df$new)), from = 10^-2, to = 10^3, add=T, col = "blue")
  text(1,400,i,cex=2)
}
