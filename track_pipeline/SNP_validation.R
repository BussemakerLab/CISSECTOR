args <- commandArgs(trailingOnly=T)
grp_name <- args[[1]]
cell <- args[[2]]

library(glmnet)
library(methods)

data_path <- "/rigel/hblab/users/vdf2104/processed/SuRE_42-45/"
grp_path <- paste0(data_path,grp_name)

element_annotation_file <- paste0(grp_path,"/element_annotation.rds")
#snp_triplets_file <- paste0(grp_path,"/SNP_triplets.rds")
bin_annotation_file <- paste0(grp_path,"/bin_annotation.rds")
stats_file <- paste0(grp_path, "/validation_stats_",cell,"_spatial.rds")

element_annotation <- readRDS(element_annotation_file)
head(element_annotation)

lib_X <- model.matrix(~ factor(lib) + 0, data = element_annotation)
colnames(lib_X) <- paste0("SuRE", c(42:45,45,42:44),"_",rep(c(1,2),each=4))

iPCR_X <- cbind(lib_X*log(element_annotation$iPCR),lib_X*log(element_annotation$iPCR)^2)
colnames(iPCR_X) <- paste("log_iPCR",rep(1:2, each = ncol(lib_X)), colnames(lib_X), sep = "_")

width_X <- element_annotation$end - element_annotation$start + 1

unpenalized <- cbind(lib_X, width_X, iPCR_X)
rm(lib_X, width_X, iPCR_X)
head(unpenalized)

#snp_triplets <- readRDS(snp_triplets_file)
#head(snp_triplets)

#snp_X <- sparseMatrix(i = snp_triplets$i, j = snp_triplets$j, x = as.numeric(snp_triplets$var) - 1.5, dims = c(nrow(element_annotation), max(snp_triplets$j)))
#rm(snp_triplets)
#dim(snp_X)

bin_annotation <- readRDS(bin_annotation_file)
head(bin_annotation)
sqrt_width <- sqrt(bin_annotation$width)
rm(bin_annotation)

bin_X <- sparseMatrix(j = unlist(mapply(function(x,y) seq.int(from = x, length.out = y), 
                                        x = element_annotation$firstBin, 
                                        y = element_annotation$nBins)), 
                      p = cumsum(c(0,element_annotation$nBins)))
dim(bin_X)

X <- cbind(unpenalized, bin_X %*% Diagonal(x = sqrt_width))
rm(bin_X,sqrt_width)
dim(X)

y <- element_annotation[,cell]
rm(element_annotation)

N <- nrow(unpenalized)
Pu <- ncol(unpenalized)

set.seed(1)
train <- sample(c(T,F), size = N, replace = T, prob = c(.9,.1))
N_train <- sum(train)

fit_null <- glm(y ~ . + 0, family = poisson(), data = data.frame(y=y[train], unpenalized[train,]))
cat("null fit complete")

pred_null <- predict(fit_null, newdata = data.frame(unpenalized), type = "response")
rm(fit_null)
length(pred_null)

logLik_null <- mean(dpois(y[!train], lambda = pred_null[!train], log = T))
logLik_sat <- mean(dpois(y[!train], lambda = y[!train], log = T))
null_deviance <- 2*(logLik_sat - logLik_null)
null_deviance

res <- y - pred_null
l1_max <- max(abs(colSums(res[train]*X[train,-(1:ncol(unpenalized))])))
rm(pred_null, res)
l1_max

penalty_factor <- rep(c(0,1), c(ncol(unpenalized),ncol(X) - ncol(unpenalized)))
u_names <- colnames(unpenalized)
rm(unpenalized)

alphas <- 2^-(1:10)

stats <- matrix(0, nrow = 0, ncol = Pu+3, dimnames = list(NULL,c("alpha", "logLambda","devRatio", u_names)))

logLambda_min <- ceiling(log(l1_max/alphas[1])) - 8

for(alpha in alphas){
  logLambda_max <- ceiling(log(l1_max/alpha))
  lambda <- exp(rev(seq(logLambda_min,logLambda_max,.2)))
  
  
  fit <- glmnet(X[train,], y = y[train], family = "poisson", alpha = alpha, lambda = lambda/N_train, intercept = F, standardize = F, penalty.factor = penalty_factor)
  cat("fit complete\n")
  
  pred_model <- predict(fit, newx = X[!train,], s = lambda/N_train, type = "response")
  dim(pred_model)
  
  logLik_model <- apply(pred_model, 2, function(m) mean(dpois(y[!train], lambda = m, log = T)))
  model_deviance <- 2*(logLik_sat - logLik_model)
  model_devRatio <- 1 - (model_deviance/null_deviance)
  stats <- rbind(stats, cbind(rep(alpha,length(lambda)),log(lambda),model_devRatio,t(as.matrix(fit$beta[1:Pu,]))))
  
  logLambda_min <- log(lambda[which.max(model_devRatio)])
  saveRDS(stats, stats_file)
  rm(fit,pred_model)
}


