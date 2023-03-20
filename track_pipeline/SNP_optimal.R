args <- commandArgs(trailingOnly=T)
cell <- args[[1]]
alpha <- as.numeric(args[[2]])
logLambda <- as.numeric(args[[3]])

library(glmnet)
library(methods)

data_path <- "/rigel/hblab/users/vdf2104/processed/SuRE_42-45/"
grp_names <- c("chr22_1", "chr1_5", "chr3_25", "chr12_2", "chrX_9", "chr4_23", "chr11_16", "chr12_10", "chr13_7", "chrX_17")
grp_paths <- paste0(data_path,grp_names)

element_annotation_files <- paste0(grp_paths,"/element_annotation.rds")
bin_annotation_files <- paste0(grp_paths,"/bin_annotation.rds")
offset_coef_file <- paste0(data_path, "optimalFitOffsetCoefficients_",cell,"_spatial.rds")

coefLists <- list()

for(grpN in 1:length(grp_names)){

cat("\n\n",grp_names[grpN],"starting\n")

t0 <- proc.time()[3]
element_annotation <- readRDS(element_annotation_files[grpN])
#snp_triplets <- readRDS(snp_triplets_file)
bin_annotation <- readRDS(bin_annotation_files[grpN])

lib_X <- model.matrix(~ factor(lib) + 0, data = element_annotation)
colnames(lib_X) <- paste0("SuRE", c(42:45,45,42:44),"_",rep(c(1,2),each=4))

iPCR_X <- cbind(lib_X*log(element_annotation$iPCR),lib_X*log(element_annotation$iPCR)^2)
colnames(iPCR_X) <- paste("log_iPCR",rep(1:2, each = ncol(lib_X)), colnames(lib_X), sep = "_")

width <- element_annotation$end - element_annotation$start + 1

unpenalized <- cbind(width, lib_X, iPCR_X)
P_un <- ncol(unpenalized)
rm(lib_X, width, iPCR_X)

t1 <- proc.time()[3]
cat("\nTime to unpenalized matrix:", round((t1-t0)/60,1), "minutes\n")


sqrt_width <- sqrt(bin_annotation$width)
rm(bin_annotation)

bin_X <- sparseMatrix(j = unlist(mapply(function(x,y) seq.int(from = x, length.out = y), 
                                        x = element_annotation$firstBin, 
                                        y = element_annotation$nBins)), 
                      p = cumsum(c(0,element_annotation$nBins)))
P_bin <- ncol(bin_X)

t3 <- proc.time()[3]
cat("\nTime to bin matrix:", round((t3-t1)/60,1), "minutes\n")

X <- cbind(unpenalized,  bin_X %*% Diagonal(x = sqrt_width))
rm(bin_X, sqrt_width)

t4 <- proc.time()[3]
cat("\nTime to full design matrix:", round((t4-t3)/60,1), "minutes\n")

N <- nrow(X)
P_all <- ncol(X)
cat("\nNumber of elements:", N,"\n\tBins:",P_bin)

penalty_factor <- rep(c(0,1), c(P_un, P_all - P_un))

fit_null <- glm(y ~ . + 0, family = poisson(), data = data.frame(y=element_annotation[,cell], unpenalized))

t5 <- proc.time()[3]
cat("\nTime to null model fit:", round((t5-t4)/60,1), "minutes\n")

pred_null <- predict(fit_null, newdata = data.frame(unpenalized), type = "response")
rm(fit_null)

res <- element_annotation[,cell] - pred_null
l1_max <- max(abs(colSums(res*X[,-(1:ncol(unpenalized))])))
rm(pred_null, res)

t6 <- proc.time()[3]
cat("\nTime to L1 max:", round((t6-t5)/60,1), "minutes\n")


logLambda_max <- ceiling(log(l1_max/alpha))
lambda <- exp(rev(seq(logLambda,logLambda_max,.2)))
  
fit <- glmnet(X, y = element_annotation[,cell], family = "poisson", alpha = alpha, lambda = lambda/N, intercept = F, standardize = F, penalty.factor = penalty_factor)

t7 <- proc.time()[3]
cat("\nTime to full fit:", round((t7-t6)/60,1), "minutes\n")

coefs <- fit$beta[,length(lambda)]
rm(fit)
coefLists[[grpN]] <- list(width = coefs[1], 
                 lib_iPCR = matrix(coefs[2:P_un], 
                                   nrow = 8, 
                                   ncol = 3, 
                                   dimnames = list(names(coefs[2:9]), 
                                                   c("intercept","log_iPCR","log_iPCR_2"))))


t7 <- proc.time()[3]
cat("\nTime to full fit:", round((t7-t6)/60,1), "minutes\n")

cat("\nTotal time:", round((t7-t0)/60,1), "minutes\n")
}

coefList <- list()
coefList$width <- mean(sapply(coefLists, function(x) x$width))
lib_iPCR <- lapply(coefLists, function(x) x$lib_iPCR)
coefList$lib_iPCR <- Reduce("+",lib_iPCR)/length(grp_names)

saveRDS(coefList,offset_coef_file)

