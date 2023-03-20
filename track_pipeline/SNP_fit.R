args <- commandArgs(trailingOnly = T)
grpN <- as.integer(args[[1]])
cell <- args[[2]]
alpha <- as.numeric(args[[3]])
logLambda <- as.numeric(args[[4]])

library(glmnet)
library(methods)

data_path <- "/rigel/hblab/users/vdf2104/processed/SuRE_42-45/"


offset_coef_file <- paste0(data_path, "optimalFitOffsetCoefficients_", cell, "_spatial.rds")
offset_coefs <- readRDS(offset_coef_file)


grp_names <- list.dirs(data_path, recursive =F, full.names=F)[seq(grpN,361,40)]

for(grp_name in grp_names){


 t0 <- T0 <- proc.time()[3]
 grp_path <- paste0(data_path, grp_name)
 element_annotation_file <- paste0(grp_path,"/element_annotation.rds")
 bin_annotation_file <- paste0(grp_path,"/bin_annotation.rds")
 offset_file <- paste0(grp_path,"/offsets_",cell,"_spatial.rds")
 bin_coef_file <- paste0(grp_path,"/bin_coefs_",cell,"_spatial.rds")
  
 if(!file.exists(offset_file)){  
 element_annotation <- readRDS(element_annotation_file)
 bin_annotation <- readRDS(bin_annotation_file)
    
    
 bin_X <- sparseMatrix(j = unlist(mapply(function(x,y) seq.int(from = x, length.out = y), 
                                            x = element_annotation$firstBin, 
                                            y = element_annotation$nBins)), 
                          p = cumsum(c(0,element_annotation$nBins)))
    
    
 X <- bin_X %*% Diagonal(x = sqrt(bin_annotation$width))
 rm(bin_annotation)
    
 y <- element_annotation[,cell]
    
 offset <- (element_annotation$end - element_annotation$start + 1)*offset_coefs$width + 
        offset_coefs$lib_iPCR[element_annotation$lib,"intercept"] + 
        offset_coefs$lib_iPCR[element_annotation$lib,"log_iPCR"]*log(element_annotation$iPCR) +
        offset_coefs$lib_iPCR[element_annotation$lib,"log_iPCR_2"]*log(element_annotation$iPCR)^2
    
 rm(element_annotation)
    
 P_bin <- ncol(bin_X)
 P_all <- ncol(X)
 N <- nrow(X) 
    
    
 cat("\nGroup", grp_name,"\n\tNumber of elements:", round(N/(10^6),2), "million\n\t\tBins:",round(P_bin/(10^6),2),"million\n")
 rm(bin_X)
    
 res <- y - exp(offset)
 l1_max <- max(abs(colSums(res*X)))
 logLambda_max <- ceiling(log(l1_max/alpha))
 lambda <- exp(rev(seq(logLambda,logLambda_max,.2)))
 rm(res)
    
 cat("\n\tPre-processing complete.\n\t\tTime:",round((proc.time()[3]-t0)/60),"minutes\n")
    
      

 t0 <- proc.time()[3]
 fit <- glmnet(X, y = y, family = "poisson", alpha = alpha, lambda = lambda/N, intercept = F, standardize = F, offset = offset)
 cat("\n\tFit complete.\n\t\tTime:",round((proc.time()[3]-t0)/60),"minutes\n")
           
 t0 <- proc.time()[3]
 offsets <- data.frame(offset1 = offset, offset2 = predict(fit, newx = X, s = lambda[length(lambda)]/N, newoffset = offset)) 
 saveRDS(offsets,offset_file)
 cat("\n\tOffsets calculated and saved.\n\t\tTime:",round((proc.time()[3]-t0)/60),"minutes\n")

 t0 <- proc.time()[3]
 bin_coefs <- fit$beta[,length(lambda)]
 saveRDS(bin_coefs,bin_coef_file)
 cat("\n\tBin coefficients saved.\n\t\tTime:",round((proc.time()[3]-t0)/60),"minutes\n")
 
 cat("\n\tTotal time:",round((proc.time()[3]-T0)/60),"minutes\n")
 }
 }
   

cat("\n\nGroup",grpN,"complete.")

