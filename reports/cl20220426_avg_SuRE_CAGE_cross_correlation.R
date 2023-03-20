library(rtracklayer)
library(data.table)
library(ggplot2)
library(argparse)
library(forecast)

import_ori_bw <- function(regions, plus_bw, minus_bw, get_sense=T,
                          width=1024, mean_sd_dt=NULL){
    is_minus = strand(regions)=='-'
    regions = resize(regions, width=width, fix='center')

    plus_regions = regions[!is_minus]
    minus_regions = regions[is_minus]

    seq_plus = as.character(seqnames(plus_regions))
    seq_minus = as.character(seqnames(minus_regions))

    if (get_sense == T){
        plus_signal = import.bw(plus_bw, which=plus_regions,
                                as = "NumericList")
        minus_signal = import.bw(minus_bw, which=minus_regions,
                                 as = "NumericList")
        if (!is.null(mean_sd)){
            mean_vec = c(unlist(mean_sd_dt[seq_plus, 'mean_plus']),
                         unlist(mean_sd_dt[seq_minus, 'mean_minus']))
            sd_vec = c(unlist(mean_sd_dt[seq_plus, 'sd_plus']),
                       unlist(mean_sd_dt[seq_minus, 'sd_minus']))
        }
    } else {
        plus_signal = import.bw(minus_bw, which=plus_regions,
                                as = "NumericList")
        minus_signal = import.bw(plus_bw, which=minus_regions,
                                 as = "NumericList")
        if (!is.null(mean_sd_dt)){
            mean_vec = c(unlist(mean_sd_dt[seq_plus, 'mean_minus']),
                         unlist(mean_sd_dt[seq_minus, 'mean_plus']))
            sd_vec = c(unlist(mean_sd_dt[seq_plus, 'sd_minus']),
                       unlist(mean_sd_dt[seq_minus, 'sd_plus']))
        }
    }
    plus_m = do.call(rbind, plus_signal)
    minus_m = do.call(rbind, minus_signal)
    colnames(plus_m) = names(regions)[!is_minus]
    colnames(minus_m) = names(regions)[is_minus]

    full_m = rbind(plus_m, minus_m[,ncol(minus_m):1])
    if (!is.null(mean_sd_dt)){
        for (i in 1:nrow(full_m)){
            full_m[i,] / mean_vec[i] * sd_vec[i]
        }

    }
    rownames(full_m) = c(plus_regions$name, minus_regions$name)
    return(full_m)
}

lookup = fread("/DATA/scratch/usr/c.leemans/projects/SuRE/SuRE_K562/lookup_matrix.txt.gz")
lookup[,ensembl_id:=gsub('[.].*', '', gene_id)]


gr = GRanges(lookup$seqnames, IRanges(start=lookup$pos, width=1),
             strand=lookup$strand, name=lookup$ensembl_id)



K562_plus = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/SURE_elasticNet_allele_K562_plus.bw"
K562_minus = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/SURE_elasticNet_allele_K562_minus.bw"

HEPG2_plus = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/SURE_elasticNet_allele_HEPG2_plus.bw"
HEPG2_minus = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/SURE_elasticNet_allele_HEPG2_minus.bw"


CAGE_K562_plus_1 = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/CAGE_K562_rep1_plus_ENCFF019NPZ.bigWig"
CAGE_K562_minus_1 = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/CAGE_K562_rep1_minus_ENCFF233CVF.bigWig"

CAGE_HEPG2_plus_1 = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/CAGE_HepG2_rep1_plus_ENCFF446BEY.bigWig"
CAGE_HEPG2_minus_1 = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/CAGE_HepG2_rep1_minus_ENCFF153RHC.bigWig"

CAGE_K562_plus_2 = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/CAGE_K562_rep2_plus_ENCFF550YWP.bigWig"
CAGE_K562_minus_2 = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/CAGE_K562_rep2_minus_ENCFF570ONH.bigWig"

CAGE_HEPG2_plus_2 = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/CAGE_HepG2_rep2_plus_ENCFF526ADE.bigWig"
CAGE_HEPG2_minus_2 = "/DATA/scratch/usr/c.leemans/data/tracks/hg19/CAGE_HepG2_rep2_minus_ENCFF616VEK.bigWig"


mean_sd <- function(plus, minus){
    bw_plus = BigWigFile(plus)
    bw_minus = BigWigFile(minus)

    seq_info = as.data.frame(seqinfo(BigWigFile(K562_plus)))
    mean_plus = summary(bw_plus, type='mean', as='matrix')
    mean_minus = summary(bw_minus, type='mean', as='matrix')

    sd_plus = summary(bw_plus, type='sd', as='matrix')
    sd_minus = summary(bw_minus, type='sd', as='matrix')
    result = data.table(seqname = rownames(seq_info),
                        length = seq_info$seqlengths,
                        mean_plus = mean_plus[,1],
                        mean_minus = mean_minus[,1],
                        sd_plus = sd_plus[,1],
                        sd_minus = sd_minus[,1])
    setkey(result, 'seqname')
    return(result)
}

mean_K562 = mean_sd(K562_plus, K562_minus)
mean_HEPG2 = mean_sd(HEPG2_plus, HEPG2_minus)


glm_K562_m = import_ori_bw(gr, K562_plus, K562_minus,
                           width=1400)
glm_HEPG2_m = import_ori_bw(gr, HEPG2_plus, HEPG2_minus,
                            width=1400)


CAGE_K562_m_1 = import_ori_bw(gr, CAGE_K562_plus_1, CAGE_K562_minus_1,
                              width=1400)
CAGE_HEPG2_m_1 = import_ori_bw(gr, CAGE_HEPG2_plus_1, CAGE_HEPG2_minus_1,
                               width=1400)

CAGE_K562_m_2 = import_ori_bw(gr, CAGE_K562_plus_2, CAGE_K562_minus_2,
                              width=1400)
CAGE_HEPG2_m_2 = import_ori_bw(gr, CAGE_HEPG2_plus_2, CAGE_HEPG2_minus_2,
                               width=1400)

CAGE_K562_m = (CAGE_K562_m_1 + CAGE_K562_m_2) / 2
CAGE_HEPG2_m = (CAGE_HEPG2_m_1 + CAGE_HEPG2_m_2) / 2



glm_K562_m_as = import_ori_bw(gr, K562_plus, K562_minus, get_sense=F,
                              width=1400)
glm_HEPG2_m_as = import_ori_bw(gr, HEPG2_plus, HEPG2_minus, get_sense=F,
                               width=1400)


glm_K562_m = (glm_K562_m - mean(glm_K562_m)) / sd(glm_K562_m)
glm_HEPG2_m = (glm_HEPG2_m - mean(glm_HEPG2_m)) / sd(glm_HEPG2_m)

glm_K562_m_as = (glm_K562_m_as - mean(glm_K562_m_as)) / sd(glm_K562_m_as)
glm_HEPG2_m_as = (glm_HEPG2_m_as - mean(glm_HEPG2_m_as)) / sd(glm_HEPG2_m_as)



cc_K562 <- ccf(as.numeric(t(glm_K562_m)), as.numeric(t(CAGE_K562_m)), lag.max = 500)
cc_HEPG2 <- ccf(as.numeric(t(glm_HEPG2_m)), as.numeric(t(CAGE_HEPG2_m)), lag.max = 500)
cc_K562_as <- ccf(as.numeric(t(glm_K562_m_as)), as.numeric(t(CAGE_K562_m)), lag.max = 500)
cc_HEPG2_as <- ccf(as.numeric(t(glm_HEPG2_m_as)), as.numeric(t(CAGE_HEPG2_m)), lag.max = 500)

cc_K562 = sapply(-1023:1023,
                  function(l) {
                      mean(glm_K562_m[,188:(1024+188)][, pmax(1,1+l):pmin(1024,1024+l)] *
                           CAGE_K562_m[,188:(1024+188)][, pmax(1,1-l):pmin(1024,1024-l)])
                  })

cc_HEPG2 <- sapply(-1023:1023,
                   function(l) {
                       mean(glm_HEPG2_m[,188:(1024+188)][, pmax(1,1+l):pmin(1024,1024+l)] *
                            CAGE_HEPG2_m[,188:(1024+188)][, pmax(1,1-l):pmin(1024,1024-l)])
                   })


cc_K562_as <- sapply(-1023:1023,
                     function(l) {
                         mean(glm_K562_m_as[,188:(1024+188)][, pmax(1,1+l):pmin(1024,1024+l)] *
                              CAGE_K562_m[,188:(1024+188)][, pmax(1,1-l):pmin(1024,1024-l)])
                     })

cc_HEPG2_as <- sapply(-1023:1023,
                      function(l) {
                          mean(glm_HEPG2_m_as[,188:(1024+188)][, pmax(1,1+l):pmin(1024,1024+l)] *
                               CAGE_HEPG2_m[,188:(1024+188)][, pmax(1,1-l):pmin(1024,1024-l)])
                      })



dt = data.table(lag = -1023:1023,
                cc_K562=runMean(cc_K562, 25),
                cc_HEPG2=runMean(cc_HEPG2, 25),
                cc_K562_as=runMean(cc_K562_as, 25),
                cc_HEPG2_as=runMean(cc_HEPG2_as, 25))


dt_melt = melt(dt, id.vars='lag',
               value.name='cross_correlation')

dt_melt[, celltype:=ifelse(grepl('K562', variable), 'K562', 'HEPG2')]
dt_melt[, orientation:=factor(ifelse(grepl('_as', variable), 'antisense', 'sense'),
                              levels=c('sense', 'antisense'))]



pdf('cl20220426_cross_correlation_CAGE_celltypes2.pdf', useDingbats=F, height=5)

ggplot(dt_melt, aes(x=lag, y=cross_correlation, linetype=orientation,
       color=celltype)) +
    theme_bw() +
    geom_line() +
    coord_cartesian(xlim=c(-600,600))

dev.off()
