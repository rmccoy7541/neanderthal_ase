# to produce ASE data file, take individual files posted to GTEx Exchange and remove those with QC flags or low coverage
zcat GTEX-1117F.phased.ase.table.tsv.gz | head -1 > GTEx_MidPoint_Imputation_ASE.expression-matrixfmt-ase.tsv
zcat *.phased.ase.table.tsv.gz | sed -e1,1d | awk '{if ($22 == 0 && $23 == 0 && $24 == 0 && ($11 > 9)) print}' >> GTEx_MidPoint_Imputation_ASE.expression-matrixfmt-ase.tsv

# for European S* calls, tag SNP if the significant call is Neanderthal
# and the derived allele matches the Neanderthal allele

cat LL.callsetEUR.mr_0.99.freqs.with_snps.bed |\
  sed -e1,1d |\
  awk '{if ($22 == "neand" && $24 == "TRUE") print}' |\
  cut -f1,3-5,7-13 |\
  sort |\
  uniq > neand_tag_snps_EUR.txt

# This set is further filtered by requiring AF > 0.9 on Neanderthal-introgressed
# haplotypes and AF < 0.1 off of these haplotypes (or vice versa). This set is stored in
# tag_snps.txt.merge

LANG=C join -1 1 -2 1 tag_snps.txt.merge neand_tag_snps_EUR.txt.merge | tr ' ' '\t' > neand_tag_snps_EUR.filtered.txt
