
 export RUN_NODE=node27 && snakemake -j 2 --use-docker align/realigned_CRC0282-05-1-D.bam align/realigned_CRC0282-07-1-A.bam &> 282_vitro_node27.sh_align.slog;
 export RUN_NODE=node27 && snakemake -j 24 --use-docker mutect_paired/CRC0282-05-1-D.pass.vcf.gz mutect_paired/CRC0282-07-1-A.pass.vcf.gz &> 282_vitro_node27.sh_mutect_paired_mice.slog;
 export RUN_NODE=node27 && snakemake -j 2 --use-docker sequenza/CRC0282-05-1-D sequenza/CRC0282-07-1-A &> 282_vitro_node27.sh_sequenza.slog;
 export RUN_NODE=node27 && snakemake -j 2 --use-docker fastqc_CRC0282-05-1-D fastqc_CRC0282-07-1-A &> 282_vitro_node27.sh_fastqc.slog;
 export RUN_NODE=node27 && snakemake -j 2 --use-docker align/CRC0282-05-1-D.bam.flagstat align/CRC0282-07-1-A.bam.flagstat &> 282_vitro_node27.sh_flagstat.slog;
 export RUN_NODE=node27 && snakemake -j 2 --use-docker align/CRC0282-05-1-D.wgsmetrics align/CRC0282-07-1-A.wgsmetrics &> 282_vitro_node27.sh_wgsmetrics.slog;
 