#args <- commandArgs(trailingOnly=T)
grp_name <- 1

library(glmnet)
library(methods)

data_path <- "/Users/vince/Documents/SuRE2/processed/SuRE_42-45/"
grp_path <- paste0(data_path,"chr17_1")

element_annotation_file <- paste0(grp_path,"/element_annotation.rds")
#snp_triplets_file <- paste0(grp_path,"/SNP_triplets.rds")
bin_annotation_file <- paste0(grp_path,"/bin_annotation.rds")

element_annotation <- readRDS(element_annotation_file)

lib_X <- model.matrix(~ factor(lib) + 0, data = element_annotation)
colnames(lib_X) <- paste0("SuRE", c(42:45,45,42:44),"_",rep(c(1,2),each=4))

iPCR_X <- cbind(lib_X*log(element_annotation$iPCR),lib_X*log(element_annotation$iPCR)^2)
colnames(iPCR_X) <- paste("log_iPCR",rep(1:2, each = ncol(lib_X)), colnames(lib_X), sep = "_")

width_X <- element_annotation$end - element_annotation$start + 1

unpenalized <- cbind(lib_X, width_X, iPCR_X)
rm(lib_X, width_X, iPCR_X)

bin_annotation <- readRDS(bin_annotation_file)

sqrt_width <- sqrt(bin_annotation$width)
rm(bin_annotation)

bin_X <- sparseMatrix(j = unlist(mapply(function(x,y) seq.int(from = x, length.out = y), 
                                        x = element_annotation$firstBin, 
                                        y = element_annotation$nBins)), 
                      p = cumsum(c(0,element_annotation$nBins)))

X <- cbind(unpenalized, bin_X %*% Diagonal(x = sqrt_width))
rm(bin_X,sqrt_width)

N <- nrow(unpenalized)
Pu <- ncol(unpenalized)


set.seed(1)
train <- sample(c(T,F), size = N, replace = T, prob = c(.9,.1))
N_train <- sum(train)


for(rep in 1:5){
  t0 <- Sys.time()
  y <- element_annotation[,rep+4]
  stats_file <- paste0(grp_path, "/validation_stats_",rep,"_spatial.rds")
  
  
  fit_null <- glm(y ~ . + 0, family = poisson(), data = data.frame(y=y[train], unpenalized[train,]))
  
  pred_null <- predict(fit_null, newdata = data.frame(unpenalized), type = "response")
  rm(fit_null)
  
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
  #rm(unpenalized)
  
  alphas <- 2^-(1:10)
  
  stats <- matrix(0, nrow = 0, ncol = Pu+3, dimnames = list(NULL,c("alpha", "logLambda","devRatio", u_names)))
  
  logLambda_min <- ceiling(log(l1_max/alphas[1])) - 8
  
  t1 <- Sys.time()
  tdiff <- t1-t0
  units(tdiff) <- "mins"
  cat("\nPre-validation for rep",rep,"complete in",round(tdiff,2),"minutes\n")
  
  for(alpha in alphas){
    logLambda_max <- ceiling(log(l1_max/alpha))
    lambda <- exp(rev(seq(logLambda_min,logLambda_max,.2)))

    fit <- glmnet(X[train,], y = y[train], family = "poisson", alpha = alpha, lambda = lambda/N_train, intercept = F, standardize = F, penalty.factor = penalty_factor)
    
    pred_model <- predict(fit, newx = X[!train,], s = lambda/N_train, type = "response")
    dim(pred_model)
    
    logLik_model <- apply(pred_model, 2, function(m) mean(dpois(y[!train], lambda = m, log = T)))
    model_deviance <- 2*(logLik_sat - logLik_model)
    model_devRatio <- 1 - (model_deviance/null_deviance)
    stats <- rbind(stats, cbind(rep(alpha,length(lambda)),log(lambda),model_devRatio,t(as.matrix(fit$beta[1:Pu,]))))
    
    logLambda_min <- log(lambda[which.max(model_devRatio)])
    saveRDS(stats, stats_file)
    rm(fit,pred_model)
    t2 <- Sys.time()
    tdiff <- t2-t1
    units(tdiff) <- "mins"
    cat("Alpha fit complete in",round(tdiff,2),"minutes\n")
    t1 <- t2
  }
  tdiff <- t1-t0
  units(tdiff) <- "mins"
  cat("Validation for rep",rep,"complete in",round(tdiff,2),"minutes\n")
}
