def main():
    #args = get_args()
    # 12 in vitro + normal
    #SAMPLES=['CRC0282-01-0', 'CRC0282-05-0', 'CRC0282-07-0', 'CRC0282-01-1-A', 'CRC0282-01-1-B', 'CRC0282-01-1-E', 'CRC0282-05-1-A', 'CRC0282-05-1-C', 'CRC0282-05-1-D', 'CRC0282-07-1-A', 'CRC0282-07-1-B', 'CRC0282-07-1-E']
    # 9 in vivo
    #SAMPLES=['CRC0282-01-MI-A', 'CRC0282-01-MI-D']
    #SAMPLES=["CRC0282LMO-0-B", "CRC1078LMO-0-B","CRC1502LMO-0-B","CRC0327LMO-0-B"]
    SAMPLES=[ 'CRC1502-10-MA-C', 'CRC1502-10-MC-A', 'CRC1502-10-MC-C', 'CRC1502-10-MC-D', 'CRC1502-10-MI-A', 'CRC1502-10-MI-C', 'CRC1502-10-MI-G','plh']

    NODES=['RUN_NODE=node' + str(n) for n in range(32, 36)]# + ['RUN_NODE=node6','RUN_NODE=node7']
    SH=['bulk_node' + str(n) + '.sh' for n in range(32, 36)]# + ['1502_node6.sh','1502_node7.sh']
    at_command = """
 export {:s} && snakemake -j 2 --use-docker align/realigned_{:s}.bam align/realigned_{:s}.bam &> {:s}_align.slog;
 export {:s} && snakemake -j 8 --use-docker mutect_paired/{:s}.pass.vcf.gz mutect_paired/{:s}.pass.vcf.gz &> {:s}_mutect_paired.slog;
 export {:s} && snakemake -j 2 --use-docker sequenza/{:s} sequenza/{:s} &> {:s}_sequenza.slog;
 export {:s} && snakemake -j 2 --use-docker fastqc_{:s} fastqc_{:s} &> {:s}_fastqc.slog;
 export {:s} && snakemake -j 2 --use-docker align/{:s}.bam.flagstat align/{:s}.bam.flagstat &> {:s}_flagstat.slog;
 export {:s} && snakemake -j 2 --use-docker align/{:s}.wgsmetrics align/{:s}.wgsmetrics &> {:s}_wgsmetrics.slog;
 """
#export {:s} && snakemake -j 24 --use-docker platypus/platypus_filtered.vcf.gz;
#    at_command ="""
#    export {:s} && snakemake -j 2 --use-docker fastqc_{:s} fastqc_{:s} &> {:s}_fastqc_mice.slog;
#    export {:s} && snakemake -j 2 --use-docker align/{:s}.bam.flagstat align/{:s}.bam.flagstat &> {:s}_flagstat_mice.slog;
#    export {:s} && snakemake -j 2 --use-docker align/{:s}.wgsmetrics align/{:s}.wgsmetrics &> {:s}_wgsmetrics_mice.slog;
#    """
# at_command="""
#export {:s} && snakemake --use-docker -j 2 sequenza/{:s} sequenza/{:s} &> {:s}_sequenza.slog;
#"""
    index_sample = 0
    for n in NODES:
        with open(SH[int(index_sample/2)], 'w') as atsh:
            #if (index_sample+1 < len(SAMPLES)):
                #print('sequenza/'+SAMPLES[index_sample])
                #print('sequenza/'+SAMPLES[index_sample+1])
            atsh.write(at_command.format(*(n, SAMPLES[index_sample], SAMPLES[index_sample+1], SH[int(index_sample/2)])*24))
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