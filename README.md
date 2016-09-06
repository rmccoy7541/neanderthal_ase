Descriptions of analysis scripts
---------------------

01a-perSiteGLMM.R: For each SNP (introgressed or nonintrogressed), fit a generalized linear mixed model to estimate the allelic effect.

01b-perSiteGLMM_tissue_heterogeneity.R: For each SNP, fit a GLMM with and without a fixed effect of tissue. Calculate the Bayes factor to compare the two models.

02-getDerivedAF.sh: Extract the ancestral/derived allele calls from the 1000 genomes dataset.

03-sigByDAF.R: Compare the proportion of introgressed and nonintrogressed SNPs showing significant ASE, stratifying by derived allele frequency to control for power.

04-aseByTissue.R: Fit a GLMM to the full introgressed dataset (rather than per-SNP) with tissue as a fixed effect. Compare the coefficient estimates for different tissues.

05_sigByTissue.R: Compare proportions of up- and down-regulated SNPs per tissue.

06_divergenceByTissue.R: Get the expression-weighted divergence between modern human and Neanderthal gene sequences per tissue.


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
