library(data.table)
library(ggplot2)
library(dtplyr)
library(boot)

tissue_list <- read.table("~/Desktop/ASE/tissues/tissueList.txt", stringsAsFactors = F)$V1

read_tissue <- function(tissue_name) {
	filename <- paste("~/Desktop/ASE/tissues/ASE_results_inla_introgressed_", tissue_name, ".txt", sep = "")
	dt <- fread(filename) %>%
		setnames(., c("CHROM", "POS", "mergeID", "EST", "CI0.025", "CI0.975")) %>%
		mutate(., TISSUE = tissue_name)
	return(dt)
}

ase_results <- do.call(rbind, lapply(tissue_list, function(x) read_tissue(x))) %>%
	mutate(., brain_indicator = grepl("BRN", TISSUE))

freq <- fread("~/Desktop/ASE/AF_biallelic.txt", sep = "\t") %>%
	setnames(., c("CHROM", "POS", "REF", "ALT", "RSID", "ANCESTRAL", "GLOBAL_AF", "EAS_AF", "EUR_AF")) %>%
	mutate(., mergeID = paste(CHROM, POS, sep = "_"))
freq <- freq[, CHROM := NULL]
freq <- freq[, POS := NULL]

ase_results <- merge(ase_results, freq, "mergeID")
ase_results <- ase_results[ALT == ANCESTRAL | REF == ANCESTRAL]
ase_results$sig <- FALSE

null <- -0.04082199

ase_results[CI0.025 > null | CI0.975 < null]$sig <- TRUE

# get frequency-filtered SNPs
ase_results <- ase_results[mergeID %in% unique(dt[neandIndicator == T]$mergeID)]

# record sigificance (set credible intervals)

by_tissue <- group_by(ase_results, TISSUE) %>%
	summarise(., prop_up   = sum((CI0.025 > null & REF == ANCESTRAL) | (CI0.975 < null & ALT == ANCESTRAL))/n(), 
				 prop_down = sum((CI0.975 < null & REF == ANCESTRAL) | (CI0.025 > null & ALT == ANCESTRAL))/n(),
				 n = n()) %>%
	mutate(., brain_indicator = grepl("BRN", TISSUE))

ggplot(data = by_tissue[n > 500], aes(x = prop_up, y = prop_down, label = TISSUE, color = factor(brain_indicator))) +
	geom_point() +
	geom_text(hjust = -0.1) +
	theme_bw() +
	theme(legend.position = "none") +
	geom_segment(x = 0, y = 0, xend = 1, yend = 1, lty = "dashed", color = "red") +
	xlab("Prop. Neand. Allele Sig. Up-Regulated") +
	ylab("Prop. Neand. Allele Sig. Down-Regulated") +
	xlim(0.05, 0.12) +
	ylim(0.05, 0.12)
