
#Per study Deseq2 for Random Effects model

#Start with the unfiltered phyloseq object per pipeline

load("~/SG_otu_tax.RData")
#setwd("~/Dropbox (Second Genome)/CRC_wBaylor/Manuscript_files/Q_SG_RE/")
library(metafor)
library(data.table)
library(DESeq2)
library(phyloseq)
library(ggplot2)

SG_tax_table <- data.frame(tax_table(SG_otu_tax))

SG_otu_tax_CRC <- subset_samples(SG_otu_tax, disease_stat%in%c("carcinoma", "control") )

#Excluding controls from the study that had only adenoma and controls
SG_otu_tax_CRC <- subset_samples(SG_otu_tax_CRC, !(Study == "Brim_V13_454"))
SG_otu_tax_CRC <- subset_samples(SG_otu_tax_CRC, !(Treatment == "Yes"))

studies = c("Zack_V4_MiSeq", "WuZhu_V3_454", "Wang_V3_454", "Chen_V13_454", "Zeller_V4_MiSeq", "Weir_V4_454", "Flemer_V34_MiSeq", "Pascual_V13_454")

study_list_deseq2 = lapply(studies, function(x, input_physeq = SG_otu_tax_CRC){
  prune_samples(get_variable(input_physeq, "Study") == x, input_physeq)
})
names(study_list_deseq2) <- studies
prevalence_norm_fun_DESeq2 = function(physeq, percent_to_keep){
  n_to_keep = percent_to_keep*nsamples(physeq)
  keep = apply(as(otu_table(physeq), "matrix"), 1, function(x) sum(x >= 1L, na.rm = TRUE) >= n_to_keep)
  ps = prune_taxa(keep, physeq)
  return(ps)
}
study_deseq2_filt = lapply(study_list_deseq2, prevalence_norm_fun_DESeq2, percent_to_keep = 0.05)
names(study_deseq2_filt) <- names(study_deseq2_filt)

#DESeq2 function:

deseq2_fun <- function(deseq_physeq){
  per_study_DeSeq = phyloseq_to_deseq2(deseq_physeq, ~disease_stat)
  gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
    if(any(x < 0)){
      return(NaN)
    }
    if(zero.propagate){
      if(any(x == 0)){
        return(0)
      }
      exp(mean(log(x), na.rm = na.rm))
    } else {
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
  }
  require("DESeq2")
  geoMeans = apply(counts(per_study_DeSeq), 1, gm_mean, na.rm=TRUE, zero.propagate = FALSE)
  per_study_DeSeq = estimateSizeFactors(per_study_DeSeq, geoMeans = geoMeans)
  per_study_DeSeq <- DESeq(per_study_DeSeq)
  deseq2_res = results(per_study_DeSeq, contrast = c("disease_stat", "carcinoma", "control"))
  return(deseq2_res)
}
per_study_deseq_SG <- lapply(study_deseq2_filt, deseq2_fun)

deseq2_sigtab_fun <- function(deseq2res_obj, SG_tax_table){
  deseq2res_obj = deseq2res_obj[order(deseq2res_obj$padj, na.last=NA, decreasing = FALSE), ]
  SG_CRC = as.data.frame(deseq2res_obj)
  SG_CRC$OTU <- rownames(SG_CRC)
  SG_tax_table$OTU <- rownames(SG_tax_table)
  SG_CRC<- merge(x = SG_CRC, y = SG_tax_table[, c("OTU", "Phylum", "Family", "Genus", "Species", "Strain") ], by = "OTU", all.x=FALSE, all.y=FALSE)
  return(SG_CRC)
}
per_study_deseq_SG_CRC <- lapply(per_study_deseq_SG, deseq2_sigtab_fun, SG_tax_table = SG_tax_table)
save(per_study_deseq_SG_CRC, file="per_study_deseq_SG_CRC.RData")

########################################################################################################################################################################################################################
#Random Effects metafor code
########################################################################################################################################################################################################################

e_CRC_SG_CRC = new.env()
load("per_study_deseq_SG_CRC.RData", envir = e_CRC_SG_CRC)
deseq_list = lapply(ls(envir = e_CRC_SG_CRC), FUN = get, envir = e_CRC_SG_CRC)
names(deseq_list) <- ls(envir = e_CRC_SG_CRC)
perstudyOTUIDs = lapply(per_study_deseq_SG_CRC, rownames)
hist(table(unlist(perstudyOTUIDs)), decreasing = TRUE)
# Appear in all studies
Reduce(intersect, perstudyOTUIDs)
OTUNumberStudies = table(unlist(perstudyOTUIDs))
(OTUNumberStudies)
length(OTUNumberStudies)
#129 OTUs occur in all 8 studies while 906 occur in 5 or more studies
library("data.table")
add_study_label = function(unfiltName, e_CRC_SG_CRC){
  df = e_CRC_SG_CRC[[unfiltName]]
  dt = NULL
  dt = data.table(df)
  dt[, Study := unfiltName]
  return(dt)
}
bigdt = rbindlist(lapply(names(e_CRC_SG_CRC$per_study_deseq_SG_CRC), add_study_label, e_CRC_SG_CRC= e_CRC_SG_CRC$per_study_deseq_SG_CRC))

#bigdt = rbindlist(lapply(names(deseq_list), add_study_label, e_CRC = e_CRC))
# dcast.data.table(bigdt, OTU ~ Study, value.var = "log2FoldChange")

setnames(bigdt, "OTU", "OTU_ID")

#install.packages("metafor")
library("metafor")
run_rma = function(OTUID, dt){
  require("metafor")
  yi = dt[OTUID]$log2FoldChange
  vi = dt[OTUID]$lfcSE
  names(yi) <- names(vi) <- dt[OTUID]$Study
  rmaRes = rma(yi, vi)
  rmaRes$slab <- dt[OTUID]$Study
  weighted = FALSE
  return(rmaRes)
}
############################################################################################################################################
#If I want to retain OTUs that occured in 5 or more studies
library("data.table")
keepOTUs <- bigdt[, list(NOTU = length(unique(Study))), by = OTU_ID][(NOTU > 4)]$OTU_ID
length(keepOTUs)
setkey(bigdt, OTU_ID)
keepdt = bigdt[keepOTUs, nomatch = 0]
setkey(keepdt, OTU_ID)
theOTUs = as.character(unique(keepdt$OTU_ID))
rmaResList_SG_CRC_filt = lapply(X = theOTUs, 
                                FUN = run_rma, dt = keepdt)
names(rmaResList_SG_CRC_filt) <- theOTUs 
save(rmaResList_SG_CRC_filt, file = "rmaResList_SG_CRC_filt.RData")
# Plot some of them
#sapply(rmaResList[1:10], forest, simplify = FALSE)

############################################################################################################################################
#If I want to run the random effects model for all samples:
library("data.table")
setkey(bigdt, OTU_ID)
theOTUs = as.character(unique(bigdt$OTU_ID))
rmaResList_SG_CRC = lapply(X = theOTUs, 
                           FUN = run_rma, dt = bigdt)
names(rmaResList_SG_CRC) <- theOTUs 
save(rmaResList_SG_CRC, file = "rmaResList_SG_CRC.RData")

############################################################################################################################################
rma2dt = function(x){
  require("data.table")
  dtlist = list(data.table(Study = "RE-Model", 
                           LogFC = as.numeric(x$b), 
                           CILB=x$ci.lb, 
                           CIUB=x$ci.ub,
                           p = x$pval,
                           tau = x$tau2,
                           SE_Tau2 = x$se.tau2,
                           QE = x$QE,
                           QEp = x$QEp,
                           I2 = x$I2,
                           H2 = x$H2),
                data.table(Study = as.character(x$slab),
                           LogFC = x$yi, 
                           CILB=x$yi - 2*sqrt(x$vi),
                           CIUB=x$yi + 2*sqrt(x$vi), 
                           p = x$pval))
  rbindlist(dtlist, use.names = TRUE, fill = TRUE)
}
#rma2dt(rmaResList[[1]])
# Convert RMA objects list into a list of data.tables
rmaResDTList_SG_CRC = lapply(rmaResList_SG_CRC, rma2dt)
rmaResDTList_SG_CRC_filt = lapply(rmaResList_SG_CRC_filt, rma2dt)

############################################################################################################################################
#OTUIDs that occur in five or more studies studies
############################################################################################################################################
rmaResDT_CRC_SG_filt <- mapply(FUN = function(dt, i){dt$OTUID <- i; return(dt)},
                               dt = rmaResDTList_SG_CRC_filt, i = names(rmaResDTList_SG_CRC_filt), SIMPLIFY = FALSE)
rmadt_CRC_SG_filt = rbindlist(rmaResDT_CRC_SG_filt)
############################################################################################################################################
rmadt_CRC_SG_filt_p = rmadt_CRC_SG_filt[(Study == "RE-Model")]
rmadt_CRC_SG_filt_p[, FDR := p.adjust(p, method = "fdr")]

rmadt_CRC_SG_filt_p <- rmadt_CRC_SG_filt_p[, list(OTUID, FDR)]
setkey(rmadt_CRC_SG_filt_p, OTUID)
setkey(rmadt_CRC_SG_filt, OTUID)
rmadt_CRC_SG_filt <- rmadt_CRC_SG_filt[rmadt_CRC_SG_filt_p]
rmadt_CRC_filt_SG_df <- data.frame(rmadt_CRC_SG_filt)
############################################################################################################################################
#Merge Taxonomy
SG_tax_table$OTU <- rownames(SG_tax_table) 
SG_tax_table$Genus <- gsub("g__", "", SG_tax_table$Genus)
SG_tax_table$Phylum <- gsub("p__", "", SG_tax_table$Phylum)
SG_tax_table$Species <- gsub("s__", "", SG_tax_table$Species)
SG_tax_table$Strain <- gsub("t__", "", SG_tax_table$Strain)

SG_tax_table <- within(SG_tax_table, Taxonomy <- paste(Phylum,Genus,Species,Strain, sep=';'))
rmadt_CRC_filt_SG_df <- merge(rmadt_CRC_filt_SG_df, SG_tax_table[, c("Taxonomy", "OTU")], by.x = "OTUID", by.y = "OTU")
#rmadt_tax$Taxonomy <- gsub("\\s","", rmadt_tax$Taxonomy)
rmadt_CRC_filt_SG_df <- data.table(rmadt_CRC_filt_SG_df)
save(rmadt_CRC_filt_SG_df, file="SG_CRC_RE.csv")

############################################################################################################################################
#Plotting
#can change color or facet wrap in the following plot
rmadt_CRC_filt_SG_df$Study <- gsub("Zack_V4_MiSeq", "Zackular_V4_MiSeq", rmadt_CRC_filt_SG_df$Study)
rmadt_CRC_filt_SG_df$Study <- as.factor(rmadt_CRC_filt_SG_df$Study)
rmadt_CRC_filt_SG_df$Study <- relevel(rmadt_CRC_filt_SG_df$Study, ref="RE-Model")

rmadt_CRC_filt_SG_df$Taxonomy <- gsub("Firmicutes;Parvimonas;97otu12932;72331", "Parvimonas micra ATCC 33270", rmadt_CRC_filt_SG_df$Taxonomy)
rmadt_CRC_filt_SG_df$Taxonomy <- gsub("Proteobacteria;unclassified;unclassified;OTU3191", "Proteobacteria unc OTU3191", rmadt_CRC_filt_SG_df$Taxonomy)
rmadt_CRC_filt_SG_df$Taxonomy <- gsub("Firmicutes;Ruminococcus;97otu15279;OTU2455", "Ruminococcus unc OTU2455", rmadt_CRC_filt_SG_df$Taxonomy)
rmadt_CRC_filt_SG_df$Taxonomy <- gsub("Firmicutes;Streptococcus;anginosus;OTU1044", "Streptococcus anginosus OTU1044", rmadt_CRC_filt_SG_df$Taxonomy)
rmadt_CRC_filt_SG_df$Taxonomy <- gsub("Bacteroidetes;Parabacteroides;distasonis;OTU553", "Parabacteroides distasonis OTU553", rmadt_CRC_filt_SG_df$Taxonomy)
rmadt_CRC_filt_SG_df$Taxonomy <- gsub("Fusobacteria;Fusobacterium;unclassified;OTU2790", "Fusobacterium unc OTU2790", rmadt_CRC_filt_SG_df$Taxonomy)

ggforest = function(x){
  require("ggplot2")
  x[, interval := CIUB - CILB]
  x[, RelConf := 1/interval]
  p = ggplot(x, aes(LogFC, Study, xmax=CIUB, xmin=CILB, color = Study)) +
    #coord_cartesian(xlim=c(-2, 2)) +
    #scale_alpha_discrete(range = c(0.2, 1)) +
    geom_vline(xintercept = 0.0, linetype=2, alpha=1) +
    geom_errorbarh(alpha=0.5, color="black", height = 0.3) + 
    geom_point(aes(size = RelConf)) +
    geom_point(data = x[(Study=="RE-Model")], color = "#6D7272", size=7) +
    scale_size(range = c(2, 5), guide=FALSE) +
    facet_wrap(~Taxonomy, shrink = FALSE) + 
    theme_bw() + 
    theme(text = element_text(size=12))
  return(p)
}

CRCpalette <- c(Brim_V13_454 = "#84CF04", Chen_V13_454 = "#01B5BB", Flemer_V34_MiSeq ="#E50E63", Pascual_V13_454 = "#8F389E", Wang_V3_454 = "#DF8236", Weir_V4_454 = "#228b22", WuZhu_V3_454 = "#F1BA2F", Zackular_V4_MiSeq = "#9F832D", Zeller_V4_MiSeq = "#154db3", "RE-Model" = "#6D7272")
colScale <- scale_colour_manual(name = "Study",values = CRCpalette)

#pforest = ggforest(rmadt_CRC_filt_SG_df[(FDR < 0.2)]) 
#pforest + colScale
#13 OTUs have FDR < 0.3
#pforest = ggforest(rmadt_tax[(FDR < 0.2)]) + ggtitle("2b")
#pforest = ggforest(rmadt_tax[1:40])
#pforest
```
#show just specific OTUIDs

```{r}
rmadt_CRC_filt_SG_df[37,5] <- 8.0
tiff("CRCFig2A.tiff", width = 15, height = 7, units = "in", res = 300,
     compression = "lzw", colortype = "true")
setkey(rmadt_CRC_filt_SG_df, OTUID)
Fig2A <- ggforest(rmadt_CRC_filt_SG_df[c("OTU1167","OTU3191","OTU2455","OTU1044", "OTU553", "OTU2790")]) +	ggtitle("2A") + theme_set(theme_bw(base_size=18)) + theme(plot.title = element_text(hjust = 0)) 
Fig2A + scale_x_continuous("LogFC", limits=c(-8,8), (breaks = c(-8, -4, 0, 4, 8))) + colScale
invisible(dev.off())

write.table(rmadt_CRC_filt_SG_df, file="RE_CRC_SG.csv")

save(rmadt_CRC_filt_SG_df, file="Fig2A_SG_RE.RData")





