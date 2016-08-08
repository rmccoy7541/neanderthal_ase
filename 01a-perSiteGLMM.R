library(data.table)
library(dplyr)
library(parallel)
library(boot)
library(INLA)

# read ASE read count data
dt <- fread("/net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/GTEx_MidPoint_Imputation_ASE.expression-matrixfmt-ase.tsv", header = T, verbose = T)
dt$mergeID <- paste(dt$CHR, dt$POS, sep = "_")
	
# get list of high confidence Neandertal tag SNPs
neand <- fread("/net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/neand_tag_snps_EUR.filtered.txt") %>%
	setnames(., c("mergeID", "CHR", "POS", "ANC", "DER", "AA_freq", "AFR_freq", "AMR_freq", "EAS_freq", "EUR_freq", "PNG_freq", "SAS_freq", "NEAND_BASE"))
dt$neandIndicator <- FALSE
dt[mergeID %in% neand$mergeID]$neandIndicator <- TRUE

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


### compute ASE per SNP, separately considering introgressed tag SNPs and non-tag SNPs ###

fixed_effect_model <- function(snp, dt) {
	chrom <- unique(dt[mergeID == snp]$CHR[1])
	pos <- unique(dt[mergeID == snp]$POS[1])
	gene <- unique(dt[mergeID == snp]$GENE_ID)
	m1 <- glm(data = dt[mergeID == snp], formula = cbind(ALT_COUNT, REF_COUNT) ~ 1, family = binomial)
	ci95 <- unname(suppressMessages(confint(m1, level = 0.95)))
	ci99 <- unname(suppressMessages(confint(m1, level = 0.99)))
	coef <- summary(m1)$coefficients
	ctab <- data.table(chrom = chrom, pos = pos, snpid = snp, geneid = gene, est = coef[1], 
	                   ci0.005 = ci99[1], ci0.025 = ci95[1], ci0.975 = ci95[2], ci0.995 = ci99[2], 
	                   p = coef[4])
	return(ctab)
}

mixed_effect_model <- function(snp, dt, model_type, null) {
	chrom <- unique(dt[mergeID == snp]$CHR[1])
	pos <- unique(dt[mergeID == snp]$POS[1])
	gene <- unique(dt[mergeID == snp]$GENE_ID)
	if (model_type == "intercept_only"){
		formula = ALT_COUNT ~ 1
	} else if (model_type == "tissue_subject") {
		formula = ALT_COUNT ~ 1 + f(SUBJECT_ID, model = "iid") + f(TISSUE_ID, model = "iid")
	} else if (model_type == "tissue") {
		formula = ALT_COUNT ~ 1 + f(TISSUE_ID, model = "iid")
	} else if (model_type == "subject") {
		formula = ALT_COUNT ~ 1 + f(SUBJECT_ID, model = "iid")
	}
	m1 <- inla(formula, data = dt[mergeID == snp], family = "binomial", Ntrials = TOTAL_COUNT, 
	           quantile = c(0.005, 0.025, 0.975, 0.995))
	coef <- m1$summary.fixed		   
	ci95 <- c(coef[4], coef[5])
	ci99 <- c(coef[3], coef[6])
	m <- m1$marginals.fixed[[1]]
	lower_p <- inla.pmarginal(null, m)
	upper_p <- 1 - inla.pmarginal(null, m)
	post_pred_p <- 2 * (min(lower_p, upper_p))
	ctab <- data.table(chrom = chrom, pos = pos, snpid = snp, geneid = gene, est = coef[1], sd = coef[2],
	                   ci0.005 = ci99[1], ci0.025 = ci95[1], ci0.975 = ci95[2], ci0.995 = ci99[2],
	                   p = post_pred_p)
	return(ctab)
}

### wrapper function to decide which model to fit for each SNP ###

glmer_ase <- function(snp, dt, null) {
	# if one observation, fit model with only intercept
	if (nrow(dt[mergeID == snp]) == 1) {
		ctab <- mixed_effect_model(snp, dt, "intercept_only", null)
	# if > 1 levels of both subject and tissue, fit model with two random effects 
  } else if (length(unique(dt[mergeID == snp]$TISSUE_ID)) > 1 & length(unique(dt[mergeID == snp]$SUBJECT_ID)) > 1) {
		ctab <- mixed_effect_model(snp, dt, "tissue_subject", null)
	# if only one level of subject, only include tissue random effect
  } else if (length(unique(dt[mergeID == snp]$TISSUE_ID)) > 1 & length(unique(dt[mergeID == snp]$SUBJECT_ID)) == 1) {
		ctab <- mixed_effect_model(snp, dt, "tissue", null)
	# if only one level of tissue, only include subject random effect
  } else if (length(unique(dt[mergeID == snp]$TISSUE_ID)) == 1 & length(unique(dt[mergeID == snp]$SUBJECT_ID)) > 1) {
		ctab <- mixed_effect_model(snp, dt, "subject", null)
  } else {
	  stop("Unknown error.")
  }
	message(paste("Processed SNP:", snp))
	return(ctab)
}

null <- logit(median(1 - dt$REF_RATIO))

p_table_introgressed <- do.call(rbind, 
                                mclapply(unique(dt[neandIndicator == T]$mergeID), 
                                         function(x) glmer_ase(x, dt, null), mc.cores = 12))

p_table_nonintrogressed <- do.call(rbind, 
                                   mclapply(unique(dt[neandIndicator == F]$mergeID), 
					                                  function(x) glmer_ase(x, dt, null), 
					                                  mc.cores = 12))

### Benjamini–Hochberg procedure for controlling FDR ###

bh_fdr <- function(p_table, fdr) {
	setorder(p_table, p)
	p_table$rank <- 1:nrow(p_table)
	m <- nrow(p_table)
	p_table <- p_table %>%
		mutate(., bh_crit = (rank / m) * fdr)
	max_rank <- max(p_table[p < bh_crit]$rank)
	p_table$sig <- FALSE
	p_table[rank <= max_rank]$sig <- TRUE
	return(p_table)
}

fdr_sig <- function(p_table, fdr) {
	p_table_bh <- bh_fdr(p_table, fdr)
	sig_prop <- mean(p_table_bh$sig)
	return(data.table(fdr = fdr, sig_prop))
}

p_table_introgressed$neandIndicator <- TRUE
p_table_nonintrogressed$neandIndicator <- FALSE

p_table <- rbind(p_table_introgressed, p_table_nonintrogressed)
p_table_0.1 <- bh_fdr(p_table, 0.1)

fdr_results <- do.call(rbind, lapply(10^seq(-10, 0, 1), function(x) fdr_sig(p_table[neandIndicator == T], x)))

write.table(p_table_0.1, "~/p_table.txt", quote = F, row.names = F, col.names = T)

### get ancestral allele states ###

freq <- fread("/net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/AF_biallelic.txt", sep = "\t") %>%
	setnames(., c("CHROM", "POS", "REF", "ALT", "RSID", "ANCESTRAL", "GLOBAL_AF", "EAS_AF", "EUR_AF")) %>%
	mutate(., snpid = paste(CHROM, POS, sep = "_"))

freq <- freq[, CHROM := NULL]
freq <- freq[, POS := NULL]

ase_results <- merge(p_table_0.1, freq, "snpid")
# remove ambiguous and low confidence ancestral allele calls
ase_results <- ase_results[ANCESTRAL %in% toupper(letters)]
ase_results <- ase_results[REF == ANCESTRAL | ALT == ANCESTRAL]

# record direction of ASE with respect to derived (i.e. Neanderthal) allele
ase_results$direction <- as.character(NA)
ase_results[(est.mean > null & REF == ANCESTRAL) | (est.mean < null & ALT == ANCESTRAL)]$direction <- "up"
ase_results[(est.mean < null & REF == ANCESTRAL) | (est.mean > null & ALT == ANCESTRAL)]$direction <- "down"

write.table(ase_results, "~/p_table_annotated.txt", quote = F, row.names = F, col.names = T)
