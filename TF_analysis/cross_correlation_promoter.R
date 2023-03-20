#!/usr/bin/env Rscript
library(rtracklayer)
library(data.table)
library(ggplot2)
library(argparse)
library(forecast)

parser <- ArgumentParser(description='cross correlation plots')
parser$add_argument('--bw', help='bigWig file')
parser$add_argument('--plus', help='SuRE plus coefficients')
parser$add_argument('--minus', help='SuRE minus coefficients')
parser$add_argument('--out', help='output')
parser$add_argument('--bed', help='bed file with info on orientation')

args = parser$parse_args()

save(args, file='test.RData')


import_ori_bw <- function(regions, plus_bw, minus_bw, get_sense=T){
    is_minus = strand(regions)=='-'
    plus_regions = regions[is_minus]
    minus_regions = regions[!is_minus]
    if (get_sense == T){
        plus_signal = import.bw(plus_bw, which=plus_regions,
                                as = "NumericList")
        minus_signal = import.bw(minus_bw, which=minus_regions,
                                 as = "NumericList")
    } else {
        plus_signal = import.bw(minus_bw, which=plus_regions,
                                as = "NumericList")
        minus_signal = import.bw(plus_bw, which=minus_regions,
                                 as = "NumericList")
    }
    plus_m = do.call(rbind, plus_signal)
    minus_m = do.call(rbind, minus_signal)
    colnames(plus_m) = names(regions)[!is_minus]
    colnames(minus_m) = names(regions)[is_minus]

    full_m = rbind(plus_m, minus_m[,ncol(minus_m):1])
    return(full_m)
}

bed_gr = import.bed(args$bed)

glm_m = import_ori_bw(bed_gr, args$plus, args$minus)
aff_m = import_ori_bw(bed_gr, args$bw, args$bw)


aff_m = (aff_m - mean(aff_m)) / sd(aff_m)
glm_m = (glm_m - mean(glm_m)) / sd(glm_m)

cc <- sapply(-1023:1023,
             function(l) {
                 mean(aff_m[, pmax(1,1+l):pmin(1024,1024+l)] *
                      glm_m[, pmax(1,1-l):pmin(1024,1024-l)])
             })

ac <- sapply(-1023:1023,
          function(l) {
              mean(aff_m[, pmax(1,1+l):pmin(1024,1024+l)] *
                   aff_m[, pmax(1,1-l):pmin(1024,1024-l)])
          })

saveRDS(data.table(cc=cc, ac=ac), args$out)
