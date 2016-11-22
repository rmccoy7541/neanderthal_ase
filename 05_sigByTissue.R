library(data.table)
library(magrittr)
library(ggplot2)
library(boot)
library(dplyr)
library(gplots)
library(RColorBrewer)

# loop through files that result from applying script 01a-perSiteGLMM.R to 
# data subset by tissue
tissue_list <- read.table("~/Desktop/ASE/tissues/tissueList.txt", stringsAsFactors = F)$V1

read_tissue <- function(tissue_name) {
  filename <- paste("~/Desktop/ASE/tissues/ASE_results_inla_introgressed_", tissue_name, ".txt", sep = "")
  dt <- fread(filename) %>%
    setnames(., c("CHROM", "POS", "mergeID", "EST", "CI0.025", "CI0.975"))
  dt[, TISSUE := tissue_name]
  return(dt)
}

ase_results <- do.call(rbind, lapply(tissue_list, function(x) read_tissue(x)))
ase_results[, brain_indicator := grepl("BRN", TISSUE)]

freq <- fread("~/Desktop/ASE/AF_biallelic.txt", sep = "\t") %>%
	setnames(., c("CHROM", "POS", "REF", "ALT", "RSID", "ANCESTRAL", "GLOBAL_AF", "EAS_AF", "EUR_AF"))
freq[, mergeID := paste(CHROM, POS, sep = "_")]
freq <- freq[, CHROM := NULL]
freq <- freq[, POS := NULL]

ase_results <- merge(ase_results, freq, "mergeID")
ase_results <- ase_results[ALT == ANCESTRAL | REF == ANCESTRAL]

neand <- fread("~/Desktop/ASE/neand_tag_snps_EUR.filtered.txt")
ase_results <- ase_results[mergeID %in% neand$V1]

null <- -0.04082199

ase_results$sig <- as.integer(0)
ase_results[(CI0.025 > null & REF == ANCESTRAL) | (CI0.975 < null & ALT == ANCESTRAL)]$sig <- 1
ase_results[(CI0.975 < null & REF == ANCESTRAL) | (CI0.025 > null & ALT == ANCESTRAL)]$sig <- -1

col.pan <- c('#00ccff','#ccf2ff','#f7f7f7','#e6ccff','purple')

# convert the data table to a matrix for gplots
y <- dcast(ase_results[grepl("BRN|TESTIS", TISSUE)], mergeID ~ TISSUE, value.var = "sig")
rownames <- y[, 1, with = F]
y <- y[, 2:ncol(y), with = F]
row.names(y) <- rownames$mergeID
y_imputed <- as.matrix(y)
y_imputed[is.na(y_imputed)] <- 0
remove_snps <- apply(y_imputed, 1, function(x) sum(x == 0))
y_imputed <- y_imputed[which(remove_snps != 14),]

# generate Fig. 5b (upper panel)
heatmap.2(y_imputed[rev(order(rowMeans(y_imputed))),], trace = "none", col = col.pan[c(1,3,5)], 
          density.info = "none", distfun = function(x) dist(x, method = "binary"), Rowv = F, labRow = "", 
          reorderfun = function(d, w) { reorder(d, c(13, 11, 14, 9, 7, 4, 6, 2, 12, 3, 8, 1, 5, 10), mean) })
          
sig_by_tissue <- group_by(ase_results, TISSUE) %>%
  summarise(., sig_down = sum(sig == -1), sig_up = sum(sig == 1), sig = sum(sig == 1 | sig == -1), n = n())
sig_by_tissue <- data.table(sig_by_tissue)

sig_by_tissue[, down_prop := sig_down / sig]
sig_by_tissue[, up_prop := sig_up / sig]

melt_sig_by_tissue <- melt(sig_by_tissue, id.vars = "TISSUE", measure.vars = c("down_prop", "up_prop"))
mean_other_tissues <- mean(melt_sig_by_tissue[!grepl("BRN|TESTIS", TISSUE) & variable == "down_prop"]$value)

# manually reorder to match heatmap
ordered_tissues <- c(c("BRNSPC", "BRNSPC", "BRNSNG", "BRNSNG", "BRNAMY", "BRNAMY", "BRNHPP", "BRNHPP", 
                       "BRNCTXA", "BRNCTXA", "BRNCTXB", "BRNCTXB", "BRNACC", "BRNACC", "BRNPTM", "BRNPTM", 
                       "BRNNCC", "BRNNCC", "BRNCDT", "BRNCDT", "BRNHPT", "BRNHPT", "BRNCHB", "BRNCHB", "BRNCHA", 
                       "BRNCHA", "TESTIS", "TESTIS"), 
                     melt_sig_by_tissue[!grepl("BRN|TESTIS", TISSUE)]$TISSUE)
melt_sig_by_tissue$TISSUE_order <- factor(melt_sig_by_tissue$TISSUE, levels = ordered_tissues)

# generate Fig. 5b (lower panel)
ggplot(data = melt_sig_by_tissue[grepl("BRN|TESTIS", TISSUE)], 
       aes(x = TISSUE_order, y = value, fill = variable))  +
  theme_bw() +
  geom_bar(stat = "identity", position = "stack") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = col.pan[c(1,5)]) +
  geom_hline(yintercept = mean_other_tissues, color = "red", lwd = 0.5)
  
