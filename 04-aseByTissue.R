library(data.table)
library(dtplyr)
library(magrittr)
library(parallel)
library(boot)
library(INLA)

# read ASE read count data
dt <- fread("/net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/GTEx_MidPoint_Imputation_ASE.expression-matrixfmt-ase.tsv", header = T, verbose = T)
dt[, mergeID := paste(CHR, POS, sep = "_")]
	
# get list of high confidence Neandertal tag SNPs
neand <- fread("/net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/neand_tag_snps_EUR.filtered.txt") %>%
	setnames(., c("mergeID", "CHR", "POS", "ANC", "DER", "AA_freq", "AFR_freq", "AMR_freq", "EAS_freq", "EUR_freq", "PNG_freq", "SAS_freq", "NEAND_BASE"))
dt[, neandIndicator := FALSE]
dt[mergeID %in% neand$mergeID, neandIndicator := TRUE]

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

freq <- fread("/net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/AF_biallelic.txt", sep = "\t") %>%
  setnames(., c("CHROM", "POS", "REF", "ALT", "RSID", "ANCESTRAL", "GLOBAL_AF", "EAS_AF", "EUR_AF"))
freq[, mergeID := paste(CHROM, POS, sep = "_")]
freq <- freq[, CHROM := NULL]
freq <- freq[, POS := NULL]

# remove ambiguous and low confidence ancestral allele calls, restrict to biallelic SNPs
dt <- merge(dt, freq, "mergeID")
dt <- dt[ANCESTRAL %in% toupper(letters)]
dt <- dt[REF == ANCESTRAL | ALT == ANCESTRAL]
dt[, DERIVED_COUNT := as.integer(NA)]
dt[REF == ANCESTRAL, DERIVED_COUNT := ALT_COUNT]
dt[ALT == ANCESTRAL, DERIVED_COUNT := REF_COUNT]

formula = DERIVED_COUNT ~ 1 + f(SUBJECT_ID, model = "iid") + f(GENE_ID, model = "iid") + TISSUE_ID
m0 <- inla(formula, data = dt[neandIndicator == TRUE], family = "binomial", Ntrials = TOTAL_COUNT, quantile = c(0.005, 0.025, 0.975, 0.995))
m0_results <- data.table(TISSUE_ID = rownames(m0$summary.fixed), m0$summary.fixed) %>% setorder(., mean)
m0_results[!(TISSUE_ID %in% n_samples[n < 10]$TISSUE_ID)]

# assign data to K = 5 subsets by gene for cross-validation
gene_list <- data.table(GENE_ID = unique(dt[neandIndicator == TRUE]$GENE_ID))
gene_list <- gene_list[!is.na(GENE_ID)]
set.seed(1)
gene_list[, setID := sample(rep(1:5, length = nrow(gene_list)))]

m_list <- lapply(1:5, function(x) inla(formula, data = dt[neandIndicator == TRUE & (GENE_ID %in% gene_list[setID == x]$GENE_ID)], 
                                       family = "binomial", Ntrials = TOTAL_COUNT, 
                                       quantile = c(0.005, 0.025, 0.975, 0.995)))

n_samples <- group_by(dt[neandIndicator == TRUE], TISSUE_ID) %>%
  summarise(., n = length(unique(SAMPLE_ID))) %>%
  setorder(., n)
n_samples <- data.table(n_samples)
n_samples[, TISSUE_ID := paste("TISSUE_ID", TISSUE_ID, sep = "")]

get_coef <- function(model_results_list, iteration, n_samples) {
  model_results <- m_list[[iteration]]
  model_coefs <- data.table(TISSUE_ID = rownames(model_results$summary.fixed), model_results$summary.fixed, iteration = iteration)
  model_coefs <- model_coefs[TISSUE_ID != "(Intercept)"] %>% setorder(., `0.995quant`)
  model_coefs[, tissue_order := .I]
  return(model_coefs)
}

m_list_results <- do.call(rbind, lapply(1:5, function(x) get_coef(m_list, x, n_samples)))
m_list_results[, TISSUE_ID := gsub("TISSUE_ID", "", TISSUE_ID)]
m_list_results[, color_factor := 1]
m_list_results[grepl("BRN", TISSUE_ID), color_factor := 2]
m_list_results[grepl("TESTIS", TISSUE_ID), color_factor := 3]

pdf("~/gg.pdf")
ggplot(data = m_list_results, aes(x = TISSUE_ID, y = tissue_order, color = factor(color_factor))) + 
	geom_point() +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 0), legend.position = "none") +
	ylab("Order of Upper 99% Confidence Limit of Coefficient Estimate") +
	xlab("Tissue") +
	geom_hline(yintercept = 26.5, color = "red")
dev.off()

### test on covariate-matched non-introgressed controls ###

dt[, brainIndicator := grepl("BRN", TISSUE_ID)]
dt[, testisIndicator := grepl("TESTIS", TISSUE_ID)]

formula = DERIVED_COUNT ~ 1 + f(SUBJECT_ID, model = "iid") + f(GENE_ID, model = "iid") + f(TISSUE_ID, model = "iid") + brainIndicator
m2 <- inla(formula, data = dt[neandIndicator == TRUE], family = "binomial", Ntrials = TOTAL_COUNT, quantile = c(0.005, 0.025, 0.975, 0.995))

match_features <- group_by(dt, mergeID) %>%
	summarise(., neandIndicator = unique(neandIndicator), GENE_ID = unique(GENE_ID),
		  n_subjects = length(unique(SUBJECT_ID)), n_tissues = length(unique(TISSUE_ID)))

post_p <- function(model) {
	m <- model$marginals.fixed[[2]]
	lower_p <- inla.pmarginal(0, m)
	upper_p <- 1 - inla.pmarginal(0, m)
	post_pred_p <- 2 * (min(lower_p, upper_p))
	return(post_pred_p)
}

control_regression_brain <- function(match_features, dt) {
	neand <- match_features[neandIndicator == T][sample(nrow(match_features[neandIndicator == T]))][!duplicated(GENE_ID)]
	cntrl <- match_features[neandIndicator == F][sample(nrow(match_features[neandIndicator == F]))][!duplicated(GENE_ID)]
	match_features_indep <- rbind(neand, cntrl)
	matches <- Match(Tr = match_features_indep$neandIndicator, X = match_features_indep[, 4:5, with = F], replace = FALSE)
	obs_snps <- match_features_indep[matches$index.treated,]$mergeID
	control_snps <- match_features_indep[matches$index.control,]$mergeID	
	r_subject <- cor.test(match_features_indep[matches$index.treated,]$n_subjects, match_features_indep[matches$index.control,]$n_subjects)$estimate
	r_tissue <- cor.test(match_features_indep[matches$index.treated,]$n_tissues, match_features_indep[matches$index.control,]$n_tissues)$estimate
	formula = DERIVED_COUNT ~ 1 + f(SUBJECT_ID, model = "iid") + f(GENE_ID, model = "iid") + f(TISSUE_ID, model = "iid") + brainIndicator
	m_neand <- inla(formula, data = dt[mergeID %in% obs_snps], family = "binomial", Ntrials = TOTAL_COUNT)
	m_cntrl <- inla(formula, data = dt[mergeID %in% control_snps], family = "binomial", Ntrials = TOTAL_COUNT)
	mean_n <- m_neand$summary.fixed[["mean"]][2]
	mean_c <- m_cntrl$summary.fixed[["mean"]][2]
	sd_n <- m_neand$summary.fixed[["sd"]][2]
	sd_c <- m_cntrl$summary.fixed[["sd"]][2]
	p_n <- post_p(m_neand)
	p_c <- post_p(m_cntrl)
	results <- data.table(r_subject = r_subject, r_tissue = r_tissue, mean_n = mean_n, mean_c = mean_c,
        		      sd_n = sd_n, sd_c = sd_c, p_n = p_n, p_c = p_c)
	write.table(results, "~/brain_control.txt", quote = F, row.names = F, col.names = F, append = T)
	return(results)
}

set.seed(123)
control_brain <- do.call(rbind, lapply(1:1000, function(x) control_regression_brain(match_features, dt)))

### plot ASE by tissue ###

nSamples <- group_by(dt[neandIndicator == TRUE], TISSUE_ID) %>%
	summarise(., length(unique(SAMPLE_ID))) %>%
	setnames(., c("TISSUE_ID", "nSamples")) %>%
	setorder(., nSamples)

byTissue <- data.table(cbind(rownames(m1$summary.fixed), m1$summary.fixed)) %>%
	setnames(., c("TISSUE_ID", "mean", "sd", "quant0.005", "quant0.025", "quant0.975", "quant0.995", "mode", "kld")) %>%
	mutate(., TISSUE_ID = gsub("TISSUE_ID", "", TISSUE_ID))

byTissue[1, 1] <- "ADPSBQ"
refMean <- byTissue[TISSUE_ID == "ADPSBQ"]$mean
byTissue[-1,]$mean <- byTissue[-1,]$mean + refMean
byTissue[-1,]$quant0.005 <- byTissue[-1,]$quant0.005 + refMean
byTissue[-1,]$quant0.025 <- byTissue[-1,]$quant0.025 + refMean
byTissue[-1,]$quant0.975 <- byTissue[-1,]$quant0.975 + refMean
byTissue[-1,]$quant0.995 <- byTissue[-1,]$quant0.995 + refMean

# require at least 10 samples per tissue
byTissue <- byTissue[TISSUE_ID %in% nSamples[nSamples >= 10]$TISSUE_ID]

byTissue$brainIndicator <- FALSE
byTissue[grepl("BRN", TISSUE_ID),]$brainIndicator <- TRUE

byTissue$TISSUE_ID <- factor(byTissue$TISSUE_ID, levels = byTissue$TISSUE_ID[order(byTissue$mean)])

limits <- aes(x = TISSUE_ID, ymin = exp(quant0.005) / (1 + exp(quant0.005)), ymax = exp(quant0.995) / (1 + exp(quant0.995)))
limits2 <- aes(x = TISSUE_ID, ymin = exp(quant0.025) / (1 + exp(quant0.025)), ymax = exp(quant0.975) / (1 + exp(quant0.975)))
ggplot(data = byTissue, aes(x = TISSUE_ID, y = exp(mean) / (1 + exp(mean)), color = factor(brainIndicator))) + 
	geom_point() +
	geom_errorbar(limits) +
	theme_bw() +
	theme(axis.text.x=element_text(angle = 90, hjust = 0), legend.position = "none") +
	ylab("Prop. Reads Supporting Neand. Allele (99% CI)") +
	xlab("Tissue")
