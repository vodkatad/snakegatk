
 export RUN_NODE=node26 && snakemake -j 2 --use-docker align/realigned_CRC0282-05-1-A.bam align/realigned_CRC0282-05-1-C.bam &> 282_vitro_node26.sh_align.slog;
 export RUN_NODE=node26 && snakemake -j 24 --use-docker mutect_paired/CRC0282-05-1-A.pass.vcf.gz mutect_paired/CRC0282-05-1-C.pass.vcf.gz &> 282_vitro_node26.sh_mutect_paired_mice.slog;
 export RUN_NODE=node26 && snakemake -j 2 --use-docker sequenza/CRC0282-05-1-A sequenza/CRC0282-05-1-C &> 282_vitro_node26.sh_sequenza.slog;
 export RUN_NODE=node26 && snakemake -j 2 --use-docker fastqc_CRC0282-05-1-A fastqc_CRC0282-05-1-C &> 282_vitro_node26.sh_fastqc.slog;
 export RUN_NODE=node26 && snakemake -j 2 --use-docker align/CRC0282-05-1-A.bam.flagstat align/CRC0282-05-1-C.bam.flagstat &> 282_vitro_node26.sh_flagstat.slog;
 export RUN_NODE=node26 && snakemake -j 2 --use-docker align/CRC0282-05-1-A.wgsmetrics align/CRC0282-05-1-C.wgsmetrics &> 282_vitro_node26.sh_wgsmetrics.slog;
 