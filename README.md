Descriptions of analysis scripts
---------------------

<code>01a-perSiteGLMM.R</code>: For each SNP (introgressed or nonintrogressed), fit a generalized linear mixed model to estimate the allelic effect. Then determine significance at 10% FDR. 

<code>01b-perSiteGLMM_tissue_heterogeneity.R</code>: For each SNP, fit a GLMM with and without a fixed effect of tissue. Calculate the Bayes factor to compare the two models.

<code>02-getDerivedAF.sh</code>: Extract the ancestral/derived allele calls from the 1000 genomes dataset.

<code>03-sigByDAF.R</code>: Compare the proportion of introgressed and nonintrogressed SNPs showing significant ASE, stratifying by derived allele frequency to control for power. [**Fig. 2**, **Fig. 4a**, **Fig. 4b**]

<code>04a-aseByTissue.R</code>: Fit a GLMM to the full introgressed dataset (rather than per-SNP) with tissue as a fixed effect. Compare the coefficient estimates for different tissues. [**Fig. 5a**]

<code>04b-aseByTissue_matched_control.R</code>: Fit the same model to equal-sized samples of covariate-matched non-introgressed control SNPs to further evaluate signficance of downregulation of Neanderthal alleles.

<code>05_sigByTissue.R</code>: Compare proportions of up- and down-regulated SNPs per tissue. [**Fig. 5b**]

<code>06_divergenceByTissue.R</code>: Get the expression-weighted divergence between modern human and Neanderthal gene sequences per tissue. [**Fig. 5c**]


Links to data sources
---------------------

GTEx Portal
http://www.gtexportal.org/

GTEx dbGaP Page <br/>
http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v6.p1

Altai Neanderthal VCFs <br/>
http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/

1000 Genomes VCFs <br/>
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
