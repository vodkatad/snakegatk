PCGR=/mnt/trcanmed/snaketree/prj/snakegatk/local/src/pcgr/
python3 ${PCGR}/pcgr.py --pcgr_dir ${PCGR} \
--output_dir ./ \
--sample_id CRC1588LMO_nocnapy \
--genome_assembly grch38 \
--conf ${PCGR}/examples/example_COAD.toml \
--input_vcf $1 \
--tumor_site 9 \
--tumor_purity 1 \
--tumor_ploidy 3 \
--include_trials \
--assay TARGETED \
--estimate_msi_status \
--estimate_tmb \
--tumor_only \
--no_vcf_validate \
--target_size_mb 0.4 \
