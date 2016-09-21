#get list of variants from ASE file
cat GTEx_MidPoint_Imputation_ASE.expression-matrixfmt-ase.tsv | sed -e1,1d | awk '{print $1"\t"$2}' | uniq | sort | uniq > positions_overlap.txt
cat positions_overlap.txt | awk '{print > $1".txt"}'

#extract these records from VCF files in /net/akey/vol2/wqfu/nobackup/1KGP using positions-overlap in vcftools
#use ancestral allele info in AA field
for chrom in {1..22}
do 
  vcftools --positions-overlap $chrom.txt \
  --gzvcf /net/akey/vol2/wqfu/nobackup/1KGP/ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
  --recode \
  --recode-INFO-all \
  --out $chrom
done

#use VariantsToTable to parse the info field
#lower case indicates low confidence call of ancestral allele
# also extract continental allele frequencies
for chrom in {1..22}
do 
  java -jar -Xmx1G ~/progs/GenomeAnalysisTK.jar \
	-R reference/hs37d5.fa \
	-T VariantsToTable \
	-V $chrom.recode.vcf \
	-F CHROM \
	-F POS \
	-F REF \
	-F ALT \
	-F ID \
	-F AA \
	-F AF \
	-F EAS_AF \
	-F EUR_AF \
	--allowMissingData \
	-o $chrom.results.table;
done

for chrom in {1..22}
do
	cat $chrom.results.table | sed 's/|||//g' >> /net/akey/vol1/home/rcmccoy/neanderthal_ase/2015_12_15/data/AF.tsv
done

#use these data to match introgressed and non-introgressed variants on derived allele frequency
