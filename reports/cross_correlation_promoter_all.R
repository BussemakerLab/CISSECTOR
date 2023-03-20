library(reshape2)
library(ggplot2)
library(data.table)
library(gridExtra)
library(XML)
library(data.table)
library(stringr)

parse_psam <- function(file){
    data = xmlToList(xmlParse(file))
    psam_split = strsplit(data$psam, "\n")[[1]]
    floats = lapply(psam_split, function(x){
        m = str_match_all(x, "[0-9.e-]+")
        return(as.numeric(m[[1]]))
    })
    psam = do.call(cbind, floats)[-5,-1]
    rownames(psam) = c('A', 'C', 'G', 'T')
    name_match = str_match(file, ".*/(.*)_([0-9]+).xml")[2:3]
    return(list(TF=name_match[1], motif= paste(name_match, collapse='_'),
                psam=psam, nsites = as.numeric(data$meta$nsites),
                e=as.numeric(data$meta$`e-value`)))
}



file_list = list.files(paste0("/DATA/scratch/usr/c.leemans/Programs/",
                              "REDUCE_Suite/data/PSAMs/Natoli_Curated/"),
                       pattern=".xml", full.names=T)

psam_info_list = lapply(file_list, parse_psam)


psam_info = do.call(rbind.data.frame,
                    c(lapply(psam_info_list,
                             function(x){x[c("TF", "motif", "nsites", "e")]}),
                      stringsAsFactors=F))

names(psam_info_list) = psam_info$motif

psam_info_dt = data.table(psam_info)

best_tf = psam_info_dt[,.SD[order(nsites,e,decreasing=T)[1],], by='TF']


psam_list = lapply(best_tf$motif, function(x){
                      psam_info_list[[x]]$psam
                  })
names(psam_list) = best_tf$TF


cc_dir = paste0('/DATA/scratch/usr/c.leemans/projects/SuRE/',
                'SuRE_K562/promoter_affinity_20210129/cross-correlation/')

file_list = list.files(cc_dir, pattern='.rds')

cc_list = lapply(file_list, function(x){
    readRDS(paste0(cc_dir, x))
})

cc_df = do.call(rbind,
                lapply(cc_list, function(x){
                    x[,cc]
                }))

colnames(cc_df) = as.character(-1023:1023)
rownames(cc_df) = gsub('.rds', '', file_list)



cc_melt = melt(cc_df, varnames=c('TF', 'pos'), stringsAsFactors=F)


cc_melt = as.data.table(cc_melt[abs(cc_melt$pos)<500, ])

cc_melt[, correlation:="cross-\ncorrelation"]


ac_df = do.call(rbind,
                lapply(cc_list, function(x){
                    x[,ac]
                }))

colnames(ac_df) = colnames(cc_df)
rownames(ac_df) = rownames(cc_df)


ac_melt = melt(ac_df, varnames=c('TF', 'pos'), stringsAsFactors=F)

ac_melt = as.data.table(ac_melt[abs(ac_melt$pos)<500 & ac_melt$pos !=0, ])

ac_melt[, correlation:="auto-\ncorrelation"]

cc_dt = rbind(cc_melt, ac_melt)

max = cc_dt[,list(max=max(value[correlation=="cross-\ncorrelation"])), by='TF']

levels = max[order(max), TF]
cc_dt[, TF:=factor(TF, levels=levels)]



cons_list = list.files(cc_dir, pattern='.RDS')


cons_cc_list = lapply(cons_list, function(x){
    readRDS(paste0(cc_dir, x))
})


cons_cc_df = do.call(rbind,
                lapply(cons_cc_list, function(x){
                    x[,cc]
                }))
colnames(cons_cc_df) = as.character(-1023:1023)
rownames(cons_cc_df) = gsub('.RDS', '', cons_list)

cons_cc_melt = melt(cons_cc_df, varnames=c('cons_score', 'pos'),
                    stringsAsFactors=F)
cons_cc_melt$correlation = "cross-\ncorrelation"
cons_ac_df = do.call(rbind,
                lapply(cons_cc_list, function(x){
                    x[,ac]
                }))

colnames(cons_ac_df) = colnames(cons_cc_df)
rownames(cons_ac_df) = rownames(cons_cc_df)

cons_ac_melt = melt(cons_ac_df, varnames=c('cons_score', 'pos'),
                    stringsAsFactors=F)
cons_ac_melt$correlation = "auto-\ncorrelation"

cons_cc = rbind(cons_cc_melt, cons_ac_melt[cons_ac_melt$pos !=0, ])

pdf('cl20210210_cross_correlation_promoter_heatmap.pdf',
    useDingbats=F, height=10, width=10)

ggplot(cc_dt[correlation=='cross-\ncorrelation', ], aes(x=pos, y=TF, fill=value)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradient2(midpoint = 0, low='blue', high='red')

dev.off()




pdf('cl20210210_cross_correlation_promoter_individual.pdf',
    useDingbats=F, height=10, width=10)

top = rev(levels)[1:10]
plot_list = lapply(top, function(x){
    p1 = ggplot(cc_dt[TF == x, ], aes(x=pos, y=value, color=correlation)) +
        geom_line() +
        xlim(-500,500) +
        theme_bw() +
        ggtitle(x) +
        guides(color=F) +
        facet_grid(rows=vars(correlation), scales="free_y") +
        theme(axis.title=element_blank())
    p2 = ggseqlogo(psam_list[[x]]) +
        coord_fixed(ratio=1) +
        theme(text=element_text(size=10))
    return(arrangeGrob(p1,p2, heights=c(2,1)))
})
do.call(grid.arrange, c(plot_list, top='top 20 at center'))

bottom = levels[1:10]

plot_list = lapply(bottom, function(x){
    p1 = ggplot(cc_dt[TF == x, ], aes(x=pos, y=value, color=correlation)) +
        geom_line() +
        xlim(-500,500) +
        theme_bw() +
        ggtitle(x) +
        guides(color=F) +
        facet_grid(rows=vars(correlation), scales="free_y") +
        theme(axis.title=element_blank())
    p2 = ggseqlogo(psam_list[[x]]) +
        coord_fixed(ratio=1) +
        theme(text=element_text(size=10))
    return(arrangeGrob(p1,p2, heights=c(2,1)))
})
do.call(grid.arrange, c(plot_list, top='bottom 20 at center'))

plot_list = lapply(unique(cons_cc$cons_score), function(x){
    ggplot(cons_cc[cons_cc$cons_score==x, ], aes(x=pos, y=value, color=correlation)) +
        geom_line() +
        xlim(-500,500) +
        theme_bw() +
        ggtitle(x) +
        guides(color=F) +
        facet_grid(rows=vars(correlation), scales="free_y") +
        theme(axis.title=element_blank())
})
do.call(grid.arrange, plot_list)
dev.off()
