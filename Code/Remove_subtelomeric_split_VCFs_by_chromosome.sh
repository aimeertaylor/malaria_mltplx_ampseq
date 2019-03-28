source /broad/software/scripts/useuse
reuse BEDTools
cd /seq/plasmodium/emilylav

bedtools subtract -a SNP_INDEL_Pf3D7_chrs_v3.combined.filtered_9999_Senegal.vcf -b Pf3D7_v3.subtelomeres.bed > subtelomeric/Senegal-subtelomeric.vcf
bedtools subtract -a SNP_INDEL.recalibrated.filtered.core_genome.LowMissingness.recode.vcf -b Pf3D7_v3.subtelomeres.bed > subtelomeric/FG-subtelomeric.vcf


grep "^#" SNP_INDEL_Pf3D7_chrs_v3.combined.filtered_9999_Senegal.vcf > subtelomeric/Senegal-header.vcf
grep "^#" SNP_INDEL.recalibrated.filtered.core_genome.LowMissingness.recode.vcf > subtelomeric/FG-header.vcf

cd subtelomeric

for x in Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 M76611 PFC10_API_IRAB
do
  grep "^${x}" Senegal-subtelomeric.vcf > Senegal-tmp.vcf
  cat Senegal-header.vcf Senegal-tmp.vcf > Senegal-${x}.vcf.gz
  grep "^${x}" FG-subtelomeric.vcf > FG-tmp.vcf
  cat FG-header.vcf FG-tmp.vcf > FG-${x}.vcf.gz
done
