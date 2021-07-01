cd mutect_paired
for f in ../../${1}_clones_2/mutect_paired/*pass.vcf.gz.tbi; do ln -s $f; done
for f in ../../${1}_clones/mutect_paired/*pass.vcf.gz.tbi; do ln -s $f; done
for f in ../../${1}_clones_2/mutect_paired/*pass.vcf.gz; do ln -s $f; done
for f in ../../${1}_clones/mutect_paired/*pass.vcf.gz; do ln -s $f; done
cd ../align
for f in ../../${1}_clones/align/realigned_*bam; do ln -s $f; done
for f in ../../${1}_clones/align/realigned_*bai; do ln -s $f; done
for f in ../../${1}_clones_2/align/realigned_*bai; do ln -s $f; done
for f in ../../${1}_clones_2/align/realigned_*bam; do ln -s $f; done
cd ..

echo "snakemake  --allowed-rules all_coverage intersect_1x intersect_callable filter_platypus pass_platypus merge_vcf platypus -pn platypus/platypus_filtered.vcf.gz depth/callable_covered.bed.gz"

#wildcard_constraints:
#    sample="[a-zA-Z0-9]+-?[0-9A-Z]+-?[A-Z0-9]+-?[A-Z0-9]?"
