library(data.table)
library(dtplyr)
library(parallel)
library(boot)
library(INLA)

# read ASE read count data
dt <- fread("/net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/GTEx_MidPoint_Imputation_ASE.expression-matrixfmt-ase.tsv", header = T, verbose = T)
dt[, mergeID := paste(CHR, POS, sep = "_")]
	
# get list of high confidence Neandertal tag SNPs
neand <- fread("/net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/neand_tag_snps_EUR.filtered.txt") %>%
	setnames(., c("mergeID", "CHR", "POS", "ANC", "DER", "AA_freq", "AFR_freq", "AMR_freq", "EAS_freq", "EUR_freq", "PNG_freq", "SAS_freq", "NEAND_BASE"))
dt$neandIndicator <- FALSE
dt[which(dt$mergeID %in% neand$mergeID)]$neandIndicator <- TRUE

# exclude extreme reference ratios
dt <- dt[REF_RATIO >= 0.1 & REF_RATIO <= 0.9]

# get sex of subjects
sex <- fread("~/neanderthal_ase/2015_12_15/data/sex.txt") %>%
  setnames(., c("SUBJECT_ID", "SEX"))
dt <- merge(dt, sex, "SUBJECT_ID")

# limit to European GTEx subjects
eur <- read.table("/net/akey/vol1/home/rcmccoy/ncbi/dbGaP-8037/files/eur.list", header = F)
eur <- as.vector(eur[['V1']])
dt <- dt[SUBJECT_ID %in% eur]

dt <- dt[neandIndicator == T]

neand_genes <- unique(dt[!is.na(GENE_ID)]$GENE_ID)
rm(dt)

#####
# read in files that contain counts of divergent sites per gene
# divergent sites in the Altai genome are defined as homozygous 
# and fixed for an allele that is at less than 2% frequency in Africa

div_per_gene <- fread("~/Downloads/div_per_gene.txt") %>% 
	setnames(., c("GENE_ID", "div_sites", "gene_length")) %>%
	mutate(., div = div_sites / gene_length)
div_per_gene <- div_per_gene[GENE_ID %in% neand_genes]

gene_lengths <- div_per_gene[, c(1,3), with = F]
missense_div_per_gene <- fread("~/Downloads/missense_stop_div_per_gene.txt") %>% 
	setnames(., c("GENE_ID", "div_sites"))
missense_div_per_gene <- merge(missense_div_per_gene, gene_lengths, "GENE_ID") %>%
	mutate(., div = div_sites / gene_length)
missense_div_per_gene <- missense_div_per_gene[GENE_ID %in% neand_genes]

gtex <- fread("~/Desktop/sunjin/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct") %>%
	setnames(., "Name", "GID")

gtex$GID <- gsub("\\..*", "", gtex$GID)

tissue_names <- data.table(full_name = colnames(gtex)[-c(1:2)], 
                           abbrev = c("ADPSBQ", "ADPVSC", "ADRNLG", "ARTAORT", "ARTCRN", "ARTTBL",
			              "BLDDER", "BRNAMY", "BRNACC", "BRNCDT", "BRNCHB", "BRNCHA", "BRNCTXA",
				      "BRNCTXB", "BRNHPP", "BRNHPT", "BRNNCC", "BRNPTM", "BRNSPC",
				      "BRNSNG", "BREAST", "LCL", "FIBRBLS", "CVXECT", "CVXEND",
				      "CLNSGM", "CLNTRN", "ESPGEJ", "ESPMCS", "ESPMSL", "FLLPNT",
				      "HRTAA", "HRTLV", "KDNCTX", "LIVER", "LUNG", "SLVRYG", "MSCLSK",
				      "NERVET", "OVARY", "PNCREAS", "PTTARY", "PRSTTE", "SKINNS",
				      "SKINS", "SNTTRM", "SPLEEN", "STMACH", "TESTIS", "THYROID",
				      "UTERUS", "VAGINA", "WHLBLD"))

colnames(gtex)[-c(1:2)] <- tissue_names$abbrev

weighted_mean_div <- function(tissue_name, div_per_gene, gtex) {
	gtex_subset <- data.table(GENE_ID = gtex$GID, RPKM = gtex[[tissue_name]])
	gtex_subset <- merge(gtex_subset, div_per_gene, "GENE_ID")
	wmean <- weighted.mean(gtex_subset$div, gtex_subset$RPKM)
	return(data.table(TISSUE_ID = tissue_name, wmean = wmean, bootstrap_wmean_wrapper(gtex_subset)))
}

bootstrap_wmean_wrapper <- function(gtex_subset) {
	set.seed(123)
	boot <- do.call(c, lapply(1:1000, function(x) bootstrap_wmean(gtex_subset)))
	quant <- quantile(boot, c(0.025, 0.975))
	return(data.table(ll = quant[1], ul = quant[2]))
}

bootstrap_wmean <- function(gtex_subset) {
	boot_sample <- gtex_subset[sample(1:nrow(gtex_subset), nrow(gtex_subset), replace = T),]
	wmean <- weighted.mean(boot_sample$div, boot_sample$RPKM)
	return(wmean)
}

div_by_tissue_gene <- do.call(rbind, lapply(tissue_names$abbrev, function(x) weighted_mean_div(x, div_per_gene, gtex))) %>%
	setorder(., wmean)
missense_div_by_tissue_gene <- do.call(rbind, lapply(tissue_names$abbrev, function(x) weighted_mean_div(x, missense_div_per_gene, gtex))) %>%
	setorder(., wmean)

div_by_tissue_gene$tissue_col <- grepl("BRN", div_by_tissue_gene$TISSUE_ID)
div_by_tissue_gene$TISSUE_ID_ORDER <- factor(div_by_tissue_gene$TISSUE_ID, levels = div_by_tissue_gene$TISSUE_ID[order(div_by_tissue_gene$wmean)])
limits <- aes(x = TISSUE_ID_ORDER, ymin = ll, ymax = ul)
ggplot(data = div_by_tissue_gene, aes(x = TISSUE_ID_ORDER, y = wmean, color = tissue_col)) +
	geom_point() +
	geom_errorbar(limits) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
	xlab("Tissue") +
	ylab("RPKM-Weighted Mean Divergence Per Site")

missense_div_by_tissue_gene$tissue_col <- grepl("BRN", missense_div_by_tissue_gene$TISSUE_ID)
missense_div_by_tissue_gene$TISSUE_ID_ORDER <- factor(missense_div_by_tissue_gene$TISSUE_ID, levels = missense_div_by_tissue_gene$TISSUE_ID[order(missense_div_by_tissue_gene$wmean)])
limits <- aes(x = TISSUE_ID_ORDER, ymin = ll, ymax = ul)
ggplot(data = missense_div_by_tissue_gene, aes(x = TISSUE_ID_ORDER, y = wmean, color = tissue_col)) +
	geom_point() +
	geom_errorbar(limits) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
	xlab("Tissue") +
	ylab("RPKM-Weighted Mean Nonsyn. Divergence Per Site")
	
merged_table <- merge(byTissue, div_by_tissue_gene, "TISSUE_ID")
limits1 <- aes(x = inv.logit(mean), ymin = ll, ymax = ul)
limits2 <- aes(y = wmean, xmin = inv.logit(quant0.025), xmax = inv.logit(quant0.975))
ggplot(data = merged_table, aes(label = TISSUE_ID, x = inv.logit(mean), y = wmean, color = tissue_col)) + 
	geom_text(hjust = -0.1) + 
	geom_errorbar(limits1, alpha = 0.5) + 
	geom_errorbarh(limits2, alpha = 0.5) + 
	theme_bw() + 
	geom_point() + 
	theme(legend.position = "none") + 
	xlab("Prop. Reads Supporting Neand. Allele") + 
	ylab("RPKM-Weighted Mean Nonsyn. Div. Per Site")

cor.test(merged_table$wmean, merged_table$mean)
