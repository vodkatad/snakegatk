
 export RUN_NODE=node21 && snakemake -j 2 --use-docker align/realigned_CRC0282-07-1-B.bam align/realigned_CRC0282-07-1-E.bam &> 282_vitro_node21.sh_align.slog;
 export RUN_NODE=node21 && snakemake -j 24 --use-docker mutect_paired/CRC0282-07-1-B.pass.vcf.gz mutect_paired/CRC0282-07-1-E.pass.vcf.gz &> 282_vitro_node21.sh_mutect_paired_mice.slog;
 export RUN_NODE=node21 && snakemake -j 2 --use-docker sequenza/CRC0282-07-1-B sequenza/CRC0282-07-1-E &> 282_vitro_node21.sh_sequenza.slog;
 export RUN_NODE=node21 && snakemake -j 2 --use-docker fastqc_CRC0282-07-1-B fastqc_CRC0282-07-1-E &> 282_vitro_node21.sh_fastqc.slog;
 export RUN_NODE=node21 && snakemake -j 2 --use-docker align/CRC0282-07-1-B.bam.flagstat align/CRC0282-07-1-E.bam.flagstat &> 282_vitro_node21.sh_flagstat.slog;
 export RUN_NODE=node21 && snakemake -j 2 --use-docker align/CRC0282-07-1-B.wgsmetrics align/CRC0282-07-1-E.wgsmetrics &> 282_vitro_node21.sh_wgsmetrics.slog;
 