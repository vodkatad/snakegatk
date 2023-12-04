bcftools merge --missing-to-ref -m none -o mutect_paired/merged_targeted.vcf.vcf.tmp mutect_paired/CRC0282LMO-0-B.pass.vcf.gz mutect_paired/CRC1078LMO-0-B.pass.vcf.gz mutect_paired/CRC1502LMO-0-B.pass.vcf.gz mutect_paired/CRC0327LMO-0-B.pass.vcf.gz
bedtools intersect -header -u -a mutect_paired/merged_targeted.vcf.vcf.tmp -b  ../../local/share/data/Pri_Met_pairs/xgen-exome-hyb-panel-v2-targets-hg38.targets.bed > mutect_paired/merged_targeted_exons.vcf
rm mutect_paired/merged_targeted.vcf.*.tmp
bcftools annotate -I +'%CHROM:%POS:%REF:%ALT' mutect_paired/merged_targeted_exons.vcf > mutect_paired/merged_exons.vcf.id
cat mutect_paired/merged_exons.vcf.id | grep -v "^##" |  perl -ane '@gt=splice(@F,9,8); $gt=""; foreach $g (@gt) { if ($.==1) {$gt.=$g."\t";} else { @afs = split(":",$g); if ($afs[2] eq ".") {$afs[2]=0;} $gt.=$afs[2]."\t";} } chop($gt) ; print $F[2]."\t".$gt."\n";' | grep -v "," > mutect_paired/merged_exons.table_nomultiallele
        
