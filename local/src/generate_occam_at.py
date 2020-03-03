def main():
    #args = get_args()
    #SAMPLES=['CRC1307-02-0', 'CRC1307-08-0', 'CRC1307-09-0', 'CRC1307-02-1-A', 'CRC1307-02-1-B', 'CRC1307-02-1-E', 'CRC1307-08-1-B', 'CRC1307-08-1-D', 'CRC1307-08-1-E', 'CRC1307-09-1-B', 'CRC1307-09-1-C', 'CRC1307-09-1-E']    
    SAMPLES= ['CRC1307-08-MA-A', 'CRC1307-08-MA-C', 'CRC1307-08-MA-F', 'CRC1307-08-MC-D', 'CRC1307-08-MC-E', 'CRC1307-08-MC-F', 'CRC1307-08-MI-A', 'CRC1307-08-MI-B', 'CRC1307-08-MI-F']

    NODES=['RUN_NODE=node' + str(n) for n in range(10, 14)]
    SH=['mutect_mice_node' + str(n) + '.sh' for n in range(10, 14)]
#     at_command ="""
# export {:s} && snakemake -j 2 --use-docker align/markedDup_{:s}.sorted.bam align/markedDup_{:s}.sorted.bam &> {:s}_align.slog;
# export {:s} && snakemake -j 2 --use-docker sequenza/{:s} sequenza/{:s} &> {:s}_sequenza.slog;
# export {:s} && snakemake -j 2 --use-docker varscan_paired/{:s}.pass.vcf.gz varscan_paired/{:s}.pass.vcf.gz &> {:s}_varscan.slog;
# export {:s} && snakemake -j 2 --use-docker fastqc_{:s} fastqc_{:s} &> {:s}_fastqc.slog;
# export {:s} && snakemake -j 2 --use-docker align/{:s}.bam.flagstat align/{:s}.bam.flagstat &> {:s}_flagstat.slog;
# export {:s} && snakemake -j 2 --use-docker align/{:s}.wgsmetrics align/{:s}.wgsmetrics &> {:s}_wgsmetrics.slog;
# """
    at_command ="""
    export {:s} && snakemake -j 24 --use-docker mutect_paired/{:s}.pass.vcf.gz mutect_paired/{:s}.pass.vcf.gz &> {:s}_mutect_paired_mice.slog;
    export {:s} && snakemake -j 2 --use-docker fastqc_{:s} fastqc_{:s} &> {:s}_fastqc_mice.slog;
    export {:s} && snakemake -j 2 --use-docker align/{:s}.bam.flagstat align/{:s}.bam.flagstat &> {:s}_flagstat_mice.slog;
    export {:s} && snakemake -j 2 --use-docker align/{:s}.wgsmetrics align/{:s}.wgsmetrics &> {:s}_wgsmetrics_mice.slog;
    """
# at_command="""
#export {:s} && snakemake --use-docker -j 2 sequenza/{:s} sequenza/{:s} &> {:s}_sequenza.slog;
#"""
    index_sample = 0
    for n in NODES:
        with open(SH[index_sample/2], 'w') as atsh:
            #if (index_sample+1 < len(SAMPLES)):
                #print('sequenza/'+SAMPLES[index_sample])
                #print('sequenza/'+SAMPLES[index_sample+1])
            atsh.write(at_command.format(*(n, SAMPLES[index_sample], SAMPLES[index_sample+1], SH[index_sample/2])*4))
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

# book occam 4 nodes 24h - from wednesday to friday for sequenza, then if free mutect 