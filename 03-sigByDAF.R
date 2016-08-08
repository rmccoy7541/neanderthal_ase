library(data.table)
library(ggplot2)
library(dtplyr)

# start with p_table_0.1 and freq from 01a-perSiteGLMM.R

p_table_0.1$sort <- factor(p_table_0.1$snpid, levels = p_table_0.1[order(p_table_0.1$est.mean)]$snpid)
limits <- aes(ymin = ci0.005, ymax = ci0.995, x = snpid)
ggplot(data = p_table_0.1[neandIndicator == T], aes(x = sort, y = est.mean)) +
	geom_point(size = 0.1) +
	geom_errorbar(limits, size = 0.1) +
	theme_bw() +
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
	panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	legend.position = "none") +
	xlab("Introgressed SNP") +
	ylab(expression(hat(beta)[0])) +
	facet_grid(. ~ sig) +
	geom_hline(yintercept = null, color = "red")

ase_results <- merge(p_table_0.1, freq, "snpid")

ase_results$DERIVED_EUR_AF <- as.numeric(NA)
ase_results[REF == ANCESTRAL]$DERIVED_EUR_AF <- ase_results[REF == ANCESTRAL]$EUR_AF
ase_results[ALT == ANCESTRAL]$DERIVED_EUR_AF <- 1 - ase_results[ALT == ANCESTRAL]$EUR_AF

### stratify by allele frequency ###
propDiffBin <- function(dt, minFreq, maxFreq) {
	dtSubset <- dt[DERIVED_EUR_AF >= minFreq & DERIVED_EUR_AF < maxFreq]
	n_sig_case <- nrow(dtSubset[neandIndicator == T & sig == T])
	n_ns_case <- nrow(dtSubset[neandIndicator == T & sig == F])
	n_sig_control <- nrow(dtSubset[neandIndicator == F & sig == T])
	n_ns_control <- nrow(dtSubset[neandIndicator == F & sig == F])
	S <- 100000
	N <- rbeta(S, n_sig_case + 1, n_ns_case + 1)
	H <- rbeta(S, n_sig_control + 1, n_ns_control + 1)
	interval <- quantile(N - H, c(0.025, 0.5, 0.975))
	return(data.table(min = minFreq, max = maxFreq, 
	n_introgressed = nrow(dtSubset[neandIndicator == T]), 
	n_nonintrogressed = nrow(dtSubset[neandIndicator == F]), 
	CI0.025 = interval[1], CI0.5 = interval[2], CI0.975 = interval[3]))
}

diffByDAF <- do.call(rbind, lapply(seq(0, 1, 0.01), function(x) propDiffBin(ase_results, x, x + 0.01)))
limits <- aes(ymin = CI0.025, ymax = CI0.975, x = (min + max) / 2)

contrasts <- ggplot(data = diffByDAF[n_introgressed >= 20], aes(x = (min + max) / 2, y = CI0.5)) + 
	geom_point() + 
	geom_errorbar(limits) +
	theme_bw() +
	ylab(expression(paste(italic("p")["N"] - italic("p")["H"], " (95 % CI)"))) +
	xlab("Derived Allele Frequency (EUR)") +
	geom_hline(yintercept = 0, lty = "dashed", color = "red")
	
# plot frequency distributions

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors <- gg_color_hue(2)

ase_results$relabel <- "Non-Introgressed"
ase_results[neandIndicator == T]$relabel <- "Neand.-Introgressed"
ggplot(data = ase_results, aes(x = DERIVED_EUR_AF, fill = sig)) + 
	geom_histogram(binwidth = 0.025, position = "stack") +
	theme_bw() + 
	ylab("Number of SNPs") +
	xlab("Derived Allele Frequency (EUR)") +
	facet_grid(relabel ~ ., scales = "free") +
	theme(legend.position = "none") +
	scale_fill_manual(values = c(colors[1], colors[2]))
