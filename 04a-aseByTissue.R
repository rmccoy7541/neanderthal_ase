library(data.table)
library(magrittr)
library(parallel)
library(boot)
library(INLA)
library(dplyr)
library(Matching)

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

formula = DERIVED_COUNT ~ -1 + f(SUBJECT_ID, model = "iid") + f(GENE_ID, model = "iid") + TISSUE_ID
m0 <- inla(formula, data = dt[neandIndicator == TRUE], family = "binomial", Ntrials = TOTAL_COUNT, quantile = c(0.005, 0.025, 0.975, 0.995))
m0_results <- data.table(TISSUE_ID = rownames(m0$summary.fixed), m0$summary.fixed)

### plot ASE by tissue ###

nSamples <- group_by(dt[neandIndicator == TRUE], TISSUE_ID) %>%
	summarise(., length(unique(SAMPLE_ID))) %>%
	setnames(., c("TISSUE_ID", "nSamples")) %>%
	setorder(., nSamples)
nSamples <- data.table(nSamples)

m0_results[, TISSUE_ID := gsub("TISSUE_ID", "", TISSUE_ID)]

# require at least 10 samples per tissue for plotting
m0_results <- m0_results[TISSUE_ID %in% nSamples[nSamples >= 10]$TISSUE_ID]
m0_results[, brainIndicator := grepl("BRN", TISSUE_ID)]
m0_results$TISSUE_ID <- factor(m0_results$TISSUE_ID, levels = m0_results$TISSUE_ID[order(m0_results$mean)])

limits <- aes(x = TISSUE_ID, ymin = exp(quant0.005) / (1 + exp(quant0.005)), ymax = exp(quant0.995) / (1 + exp(quant0.995)))
limits2 <- aes(x = TISSUE_ID, ymin = exp(quant0.025) / (1 + exp(quant0.025)), ymax = exp(quant0.975) / (1 + exp(quant0.975)))
# generate Fig. 4a
ggplot(data = m0_results, aes(x = TISSUE_ID, y = inv.logit(mean), color = factor(brainIndicator))) +
	geom_point() +
	geom_errorbar(limits) +
	theme_bw() +
	theme(axis.text.x=element_text(angle = 90, hjust = 0), legend.position = "none") +
	ylab("Prop. Reads Supporting Neand. Allele (99% CI)") +
	xlab("Tissue")