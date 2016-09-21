library(data.table)
library(dplyr)
library(ggplot2)
library(INLA)

eur <- read.table("~/Desktop/ASE/eur.list", header = F)
eur <- as.vector(eur[['V1']])

neand <- fread("~/Desktop/ASE/neand_tag_snps_EUR.filtered.txt") %>%
  setnames(., c("mergeID", "CHR", "POS", "ANC", "DER", "AA_freq", "AFR_freq", "AMR_freq", "EAS_freq", "EUR_freq", "PNG_freq", "SAS_freq", "NEAND_BASE"))

sex <- fread("~/Desktop/ASE/sex.txt") %>%
  setnames(., c("SUBJECT_ID", "SEX"))

freq <- fread("~/Desktop/ASE/AF_biallelic.txt", sep = "\t") %>%
  setnames(., c("CHROM", "POS", "REF", "ALT", "RSID", "ANCESTRAL", "GLOBAL_AF", "EAS_AF", "EUR_AF"))
freq[, mergeID := paste(CHROM, POS, sep = "_")]
freq <- freq[, CHROM := NULL]
freq <- freq[, POS := NULL]

dt <- fread("~/Desktop/ASE/GTEx_MidPoint_Imputation_ASE.expression-matrixfmt-ase.introgressed.tsv")
dt <- merge(dt, sex, "SUBJECT_ID")
dt <- dt[REF_RATIO >= 0.1 & REF_RATIO <= 0.9]
dt <- dt[SUBJECT_ID %in% eur]

dt$neandIndicator <- FALSE
dt[which(dt$mergeID %in% neand$mergeID)]$neandIndicator <- TRUE
dt <- dt[neandIndicator == T]

neand_genes <- unique(dt[!is.na(GENE_ID)]$GENE_ID)

#####

div_per_gene <- fread("~/Desktop/ASE/divergence/div_per_gene.txt") %>% 
  setnames(., c("GENE_ID", "div_sites", "gene_length"))
div_per_gene[, div := div_sites / gene_length]
div_per_gene <- div_per_gene[GENE_ID %in% neand_genes]

gene_lengths <- div_per_gene[, c(1,3), with = F]

missense_div_per_gene <- fread("~/Desktop/ASE/divergence/missense_stop_div_per_gene.txt") %>% 
  setnames(., c("GENE_ID", "div_sites"))
missense_div_per_gene <- merge(missense_div_per_gene, gene_lengths, "GENE_ID")
missense_div_per_gene[, div := div_sites / gene_length]
missense_div_per_gene <- missense_div_per_gene[GENE_ID %in% neand_genes]

#####

deleterious <- fread("~/Desktop/ASE/divergence/deleterious.txt", header = F) %>% # genes with missense SNP predicted as deleterious by SIFT or PolyPhen
  setnames(., "GENE_ID")

deleterious_intersection <- fread("~/Desktop/ASE/divergence/deleterious_intersection.txt") %>% # by both SIFT and PolyPhen
  setnames(., "GENE_ID")

# remove ambiguous and low confidence ancestral allele calls, restrict to biallelic SNPs
dt_delet <- merge(dt, freq, "mergeID")
dt_delet <- dt_delet[ANCESTRAL %in% toupper(letters)]
dt_delet <- dt_delet[REF == ANCESTRAL | ALT == ANCESTRAL]
dt_delet[, DERIVED_COUNT := as.integer(NA)]
dt_delet[REF == ANCESTRAL, DERIVED_COUNT := ALT_COUNT]
dt_delet[ALT == ANCESTRAL, DERIVED_COUNT := REF_COUNT]

dt_delet[, delet := GENE_ID %in% deleterious$GENE_ID]
dt_delet[, delet_int := GENE_ID %in% deleterious_intersection$GENE_ID]

formula = DERIVED_COUNT ~ 1 + f(SUBJECT_ID, model = "iid") + f(GENE_ID, model = "iid") + f(TISSUE_ID, model = "iid") + delet
m0 <- inla(formula, data = dt_delet[neandIndicator == TRUE], family = "binomial", Ntrials = TOTAL_COUNT, quantile = c(0.005, 0.025, 0.975, 0.995))

formula = DERIVED_COUNT ~ 1 + f(SUBJECT_ID, model = "iid") + f(GENE_ID, model = "iid") + f(TISSUE_ID, model = "iid") + delet_int
m1 <- inla(formula, data = dt_delet[neandIndicator == TRUE], family = "binomial", Ntrials = TOTAL_COUNT, quantile = c(0.005, 0.025, 0.975, 0.995))

post_p <- function(model) {
  m <- model$marginals.fixed[[2]]
  lower_p <- inla.pmarginal(0, m)
  upper_p <- 1 - inla.pmarginal(0, m)
  post_pred_p <- 2 * (min(lower_p, upper_p))
  return(post_pred_p)
}

post_p(m0)
post_p(m1)

#####

gtex <- fread("~/Desktop/ASE/sunjin/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct") %>%
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


tissues <- colnames(gtex)[-(1:2)]
mean_missense_div <- function(tissues, missense_div_per_gene, gtex) {
  gtex_subset <- data.table(GENE_ID = gtex$GID, RPKM = rowMeans(gtex[, tissues, with = F]))
  gtex_subset <- merge(gtex_subset, missense_div_per_gene, "GENE_ID")
  wmean <- weighted.mean(gtex_subset$div, gtex_subset$RPKM)
  return(data.table(wmean, bootstrap_wmean_wrapper(gtex_subset)))
}
missense_summary <- rbind(data.table(tissue = "Brain", mean_missense_div(grepl("BRN", tissues), missense_div_per_gene, gtex)),
                          data.table(tissue = "Testis", mean_missense_div("TESTIS", missense_div_per_gene, gtex)), 
                          data.table(tissue = "All Tissues", mean_missense_div(tissues, missense_div_per_gene, gtex)),
                          data.table(tissue = "All Other Tissues", mean_missense_div(tissues[!grepl("BRN|TESTIS", tissues)], missense_div_per_gene, gtex)))
limits <- aes(x = tissue, ymin = ll, ymax = ul)
# generate Fig. 4c
ggplot(data = missense_summary[!grepl("All", tissue)], aes(x = tissue, y = wmean)) +
	geom_bar(stat = "identity") +
	geom_errorbar(limits, position = "identity", width = 0.2 ) +
	theme_bw() +
	ylab("RPKM-Weighted Mean Nonsyn. Divergence") +
	xlab("Tissue") +
	geom_hline(yintercept = missense_summary[tissue == "All Other Tissues"]$wmean,
                   lty = "dashed", color = "darkred")
