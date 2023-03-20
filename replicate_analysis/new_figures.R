library(glmnet)
library(methods)

data_path <- "/Users/vince/Documents/SuRE2/processed/SuRE_42-45/"
grp_path <- paste0(data_path,"chr17_1")

element_annotation_file <- paste0(grp_path,"/element_annotation.rds")
element_annotation <- readRDS(element_annotation_file)

y <- rowSums(element_annotation[,5:7])
z <- tapply(y, list(element_annotation$iPCR, element_annotation$lib), mean)

par(mar = c(5,4,4,1),pty="s")
matplot(z[1:10,], type = "b", xlim = c(0,10), ylim = c(0,max(z[1:10,])), xlab = "iPCR (input) count", ylab = "mean barcode count")
abline(h=0,lty=2,col = "grey")
abline(v=0,lty=2,col = "grey")



bin_annotation_file <- paste0(grp_path,"/bin_annotation.rds")
coef_file <- paste0(grp_path,"/bin_coefs.rds")
bin_annotation <- readRDS(bin_annotation_file)
bin_coefs <- readRDS(coef_file)
cors <- cov.wt(bin_coefs, wt = bin_annotation$width, cor= T)$cor
rownames(cors) <- c("K562 1","K562 2","K562 3","HepG2 1","HepG2 2")
colnames(cors) <- rownames(cors)
cors2 <- cors
cors2[1,1:5] <- NA
cors2[2,2:5] <- NA
cors2[3,3:5] <- NA
cors2[4,4:5] <- NA
cors2[4,5] <- NA

library(stringr)
par(mar = c(0,0,0,0),pty="s")
image(t(cors2[5:1,]),col = colorRampPalette(c("magenta","cyan"))(30),zlim = c(0,.7), xlim = c(-0.375,1.125),ylim = c(-0.375,1.125), frame.plot = F,axes=F)
text(x = rep(c(0,.25,.5,.75),each = 4), y = rev(c(0,.25,.5,.75)), labels = str_pad(round(as.vector(cors2[-1,-5]),3),width = 5,side = "right",pad = "0"))
text(y = rev(c(0,.25,.5,.75)),x = -0.125,pos=2,labels = rownames(cors2[-1,]))
text(x = c(0,.25,.5,.75),y = -0.18,pos=1,labels = rownames(cors2[-5,]), srt = 90)


plot(0,type = "n",ylim = c(-1,4),xlim=c(-1,4),frame.plot = F,axes=F)
rect(xleft = )