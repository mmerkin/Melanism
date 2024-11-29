#! /bin/bash

input=
species=
Z=
remove=

vcftools --gzvcf $input --not-chr $Z --recode --recode-INFO-all --out ${species}

bcftools view -e 'QUAL<30' -O b ${species}.recode.vcf > ${species}_QS30.bcf
a=$(bcftools view -H ${species}_QS30.bcf | wc -l)
echo $a

bcftools view -v snps ${species}_QS30.bcf -O b > ${species}_snps.QS30.bcf
b=$(bcftools view -H ${species}_snps.QS30.bcf | wc -l)
echo $b

bcftools view -e 'INFO/DP<10' -O b ${species}_snps.QS30.bcf > ${species}_snps.QS30.DP10.bcf
c=$(bcftools view -H ${species}_snps.QS30.DP10.bcf | wc -l)
echo $c

bcftools view -e 'MAF<0.05' -O b ${species}_snps.QS30.DP10.bcf > ${species}.filtered.bcf
d=$(bcftools view -H ${species}.filtered.bcf | wc -l)
echo $d

bcftools convert -O v -o ${species}_filtered.vcf ${species}.filtered.bcf
vcftools --vcf ${species}_filtered.vcf --max-missing 0.95 --minDP 3 --mac 2 --max-alleles 2 \
--remove $remove \
--recode --recode-INFO-all --out ${species}_allfilters

vcftools --vcf ${species}_allfilters.recode.vcf --missing-indv
mawk '$5 > 0.3' out.imiss | cut -f1 > lowDP.indv
vcftools --vcf ${species}_allfilters.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out ${species}_FILTERED

echo -e "${a}\t${b}\t${c}\t${d}"

plink --vcf ${species}_FILTERED.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--vcf-half-call m --allow-no-sex \
--indep-pairwise 50 10 0.1 --out ${species}_plink


plink --vcf ${species}_FILTERED.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
--vcf-half-call m --allow-no-sex \
--extract ${species}_plink.prune.in \
--pca --out ${species}_pca
