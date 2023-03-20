library(glmnet)
library(methods)

data_path <- "/Users/vince/Documents/SuRE2/processed/SuRE_42-45/"
grp_path <- paste0(data_path,"chr17_1")

element_annotation_file <- paste0(grp_path,"/element_annotation.rds")
#snp_triplets_file <- paste0(grp_path,"/SNP_triplets.rds")
bin_annotation_file <- paste0(grp_path,"/bin_annotation.rds")
coef_file <- paste0(grp_path,"/bin_coefs.rds")

element_annotation <- readRDS(element_annotation_file)

lib_X <- model.matrix(~ factor(lib) + 0, data = element_annotation)
colnames(lib_X) <- paste0("SuRE", c(42:45,45,42:44),"_",rep(c(1,2),each=4))

iPCR_X <- cbind(lib_X*log(element_annotation$iPCR),lib_X*log(element_annotation$iPCR)^2)
colnames(iPCR_X) <- paste("log_iPCR",rep(1:2, each = ncol(lib_X)), colnames(lib_X), sep = "_")

width_X <- element_annotation$end - element_annotation$start + 1

unpenalized <- cbind(lib_X, width_X, iPCR_X)
Pu <- ncol(unpenalized)

rm(lib_X, width_X, iPCR_X)

bin_annotation <- readRDS(bin_annotation_file)

sqrt_width <- sqrt(bin_annotation$width)

bin_X <- sparseMatrix(j = unlist(mapply(function(x,y) seq.int(from = x, length.out = y), 
                                        x = element_annotation$firstBin, 
                                        y = element_annotation$nBins)), 
                      p = cumsum(c(0,element_annotation$nBins)))

X <- cbind(unpenalized, bin_X %*% Diagonal(x = sqrt_width))
penalty_factor <- rep(c(0,1), c(Pu,ncol(X) - Pu))
rm(bin_X,sqrt_width,unpenalized)

N <- nrow(X)
bin_coefs <- matrix(0,nrow = nrow(bin_annotation), ncol = 0)

for(rep in 4:5){
  stats_file <- paste0(grp_path, "/validation_stats_",rep,"_spatial.rds")
  y <- element_annotation[,rep+4]
  stats <- readRDS(stats_file)
  w <- which.max(stats[,"devRatio"])
  alpha <- stats[w,"alpha"]
  ll <- stats[stats[,"alpha"]==alpha,"logLambda"]
  ll <- ll[ll >= stats[w,"logLambda"]]
  lambda <- sort(exp(ll),decreasing = T)
  
  system.time({
  fit <- glmnet(X, y = y, family = "poisson", alpha = alpha, lambda = lambda/N, intercept = F, standardize = F, penalty.factor = penalty_factor)
  })
  
  bin_coefs <- cbind(bin_coefs,fit$beta[-c(1:Pu),length(lambda)])
  saveRDS(bin_coefs, coef_file)
  
  rm(fit)
}


cors <- cov.wt(bin_coefs, wt = bin_annotation$width, cor= T)$cor
