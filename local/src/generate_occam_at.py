def main():
    #args = get_args()
    SAMPLES=['CRC1307NLH-0-B','CRC1307-02-0', 'CRC1307-08-0', 'CRC1307-09-0', 'CRC1307-02-1-A', 'CRC1307-02-1-B', 'CRC1307-02-1-E', 'CRC1307-08-1-B', 'CRC1307-08-1-D', 'CRC1307-08-1-E', 'CRC1307-09-1-B', 'CRC1307-09-1-C', 'CRC1307-09-1-E']    
    NODES=['RUN_NODE=node' + str(n) for n in range(10, 16)]
    SH=['node' + str(n) + '_at.sh' for n in range(10, 16)]
    at_command = """
export {:s} && snakemake -j 2 --use-docker align/markedDup_{:s}.bam align/markedDup_{:s}.bam &> {:s}_align.slog;
export {:s} && snakemake -j 2 --use_docker sequenza/{:s} sequenza/{:s} &> {:s}_varscan.slog;
export {:s} && snakemake -j 2 --use-docker varscan_paired/{:s}.pass.vcf.gz varscan_paired/{:s}.pass.vcf.gz &> {:s}_varscan.slog;
export {:s} && snakemake -j 2 --use-docker fastqc_{:s} fastqc_{:s} &> {:s}_fastqc.slog;
export {:s} && snakemake -j 2 --use-docker fastqc_{:s} fastqc_{:s} &> {:s}_fastqc.slog;
export {:s} && snakemake -j 2 --use-docker align/{:s}.bam.flagstat align/{:s}.bam.flagstat &> {:s}_flagstat.slog;
export {:s} && snakemake -j 2 --use-docker align/{:s}.wgsmetrics align/{:s}.wgsmetrics &> {:s}_wgsmetrics.slog;
        """
    index_sample = 0
    for n in NODES:
        with open(SH[index_sample/2], 'w') as atsh:
            atsh.write(at_command.format(*(n, SAMPLES[index_sample], SAMPLES[index_sample+1], n)*7))

            index_sample = index_sample + 2

        
# def get_args():
#     parser = ArgumentParser(description="From a sample sheet in excel format to the json samples_map for las mdam")
#     parser.add_argument("-f", "--fastq",
#                         default="./fastq/",
#                         help="the fastq dir")
#     parser.add_argument("-i", "--input",
#                         required=True,
#                         help="the tsv sample sheet")
#     parser.add_argument("-j", "--json",
#                         required=True,
#                         help="the json to be filled")
#     args = parser.parse_args()
#     return args


if __name__ == '__main__':
    main()
