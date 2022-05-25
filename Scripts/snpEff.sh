# Snpeff

# Build Annotation file for Pecten maximus following instructions at: http://pcingola.github.io/SnpEff/se_buildingdb/
# Run it like this to make it work:
java -jar snpEff.jar build -gff3 -v -noCheckCds -noCheckProtein Pecten_maximus.genomic 2>&1 | tee GRCh37.70.build

# Run snpeff, to see if highly divergent snps are synonymous, missense or LoF
# 1. Diagnostic loci
# Make vcf:
vcftools --vcf /prj/mar-in-gen/Pecten_M/Clean/ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --positions /prj/mar-in-gen/Pecten_M/Enrichment_new/Diagnostic_loci/positions.txt --recode --out /prj/mar-in-gen/bin/snpEff/q95
# Place the vcf file in: /prj/mar-in-gen/bin/snpEff 
# Run snpEff:
java -Xmx8g -jar snpEff.jar Pecten_maximus.genomic Diagn.recode.vcf > Diag.ann.vcf
# Extract relevant fields and add them to a final Annotation file:
grep -v '^#' Diag.ann.vcf | cut -d $'\t' -f8 | cut -d '|' -f2,3 | sed 's/|/\t/g' | sed 's/MODIFIER/NA/g' > Diag_snpEff_Annot.txt
cat /prj/mar-in-gen/Pecten_M/Enrichment_new/Ann_temp_head.txt Diag_snpEff_Annot.txt > Diag_snpEff_Annot_wH.txt
paste All_Diagnostic_SNPs_Annotations.txt Diag_snpEff_Annot_wH.txt > Final_All_Diag_Annotations.txt
