
 export RUN_NODE=node23 && snakemake -j 2 --use-docker align/realigned_CRC0282-01-0.bam align/realigned_CRC0282-05-0.bam &> 282_vitro_node23.sh_align.slog;
 export RUN_NODE=node23 && snakemake -j 24 --use-docker mutect_paired/CRC0282-01-0.pass.vcf.gz mutect_paired/CRC0282-05-0.pass.vcf.gz &> 282_vitro_node23.sh_mutect_paired_mice.slog;
 export RUN_NODE=node23 && snakemake -j 2 --use-docker sequenza/CRC0282-01-0 sequenza/CRC0282-05-0 &> 282_vitro_node23.sh_sequenza.slog;
 export RUN_NODE=node23 && snakemake -j 2 --use-docker fastqc_CRC0282-01-0 fastqc_CRC0282-05-0 &> 282_vitro_node23.sh_fastqc.slog;
 export RUN_NODE=node23 && snakemake -j 2 --use-docker align/CRC0282-01-0.bam.flagstat align/CRC0282-05-0.bam.flagstat &> 282_vitro_node23.sh_flagstat.slog;
 export RUN_NODE=node23 && snakemake -j 2 --use-docker align/CRC0282-01-0.wgsmetrics align/CRC0282-05-0.wgsmetrics &> 282_vitro_node23.sh_wgsmetrics.slog;
 