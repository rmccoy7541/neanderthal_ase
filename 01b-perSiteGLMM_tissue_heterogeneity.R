library(data.table)
library(dplyr)
library(parallel)
library(boot)
library(INLA)

# read ASE read count data
dt <- fread("/net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/GTEx_MidPoint_Imputation_ASE.expression-matrixfmt-ase.tsv", header = T, verbose = T) %>%
dt[, mergeID := paste(CHR, POS, sep = "_")]

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

mixed_effect_model <- function(snp, dt, model_type) {
	if (model_type == "tissue_subject") {
		formula = ALT_COUNT ~ 1 + f(SUBJECT_ID, model = "iid") + TISSUE_ID
	} else if (model_type == "tissue") {
		formula = ALT_COUNT ~ 1 + TISSUE_ID
	}
	m1 <- inla(formula, data = dt[mergeID == snp], family = "binomial", Ntrials = TOTAL_COUNT, 
		   	     quantile = c(0.005, 0.025, 0.975, 0.995))
	return(m1)
}

null_model <- function(snp, dt, model_type) {
	if (model_type == "tissue_subject") {
		formula = ALT_COUNT ~ 1 + f(SUBJECT_ID, model = "iid")
	} else if (model_type == "tissue") {
		formula = ALT_COUNT ~ 1
	}
	m0 <- inla(formula, data = dt[mergeID == snp], family = "binomial", Ntrials = TOTAL_COUNT, 
		   	     quantile = c(0.005, 0.025, 0.975, 0.995))
	return(m0)
}

glmer_ase <- function(snp, dt) {
	chrom <- unique(dt[mergeID == snp]$CHR[1])
	pos <- unique(dt[mergeID == snp]$POS[1])
	gene <- unique(dt[mergeID == snp]$GENE_ID)
	# if one observation, fit model with only intercept
	if (length(unique(dt[mergeID == snp]$TISSUE_ID)) > 1 & length(unique(dt[mergeID == snp]$SUBJECT_ID)) > 1) {
		m1 <- mixed_effect_model(snp, dt, "tissue_subject")
		m0 <- null_model(snp, dt, "tissue_subject")
	# if only one level of subject, only include tissue random effect
	} else if (length(unique(dt[mergeID == snp]$TISSUE_ID)) > 1 & length(unique(dt[mergeID == snp]$SUBJECT_ID)) == 1) {
		m1 <- mixed_effect_model(snp, dt, "tissue")
		m0 <- null_model(snp, dt, "tissue")
	# if only one level of tissue, only include subject random effect
	} else if (length(unique(dt[mergeID == snp]$TISSUE_ID)) == 1) {
		stop("Only one tissue for this SNP.")
	} else {
		stop("Unknown error.")
	}
	b10 <- exp(m1$mlik[2] - m0$mlik[2])
	message(paste("Processed SNP:", snp, b10))
	return(data.table(snp = snp, chrom = chrom, pos = pos, gene = gene, BF = b10))
}

n_tissues <- group_by(dt[neandIndicator == T], mergeID) %>%
  summarise(., n_tissues = length(unique(TISSUE_ID)))
  
tissue_het <- do.call(rbind, lapply(n_tissues[n_tissues > 1]$mergeID, function(x) glmer_ase(x, dt[neandIndicator == T])))
