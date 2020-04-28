#!/usr/bin/env python

from sys import argv
import pysam


if __name__ == "__main__":

        if len(argv) != 4:
                print("No input bam/chr/coord passed as argument %d" % len(argv))
                exit()

        samfile = pysam.AlignmentFile(argv[1], "rb")
        # for each read that overlap the postion given
        target = int(argv[3])-1
        chrom = argv[2]
        #for pileupcolumn in samfile.pileup(chrom, target, target+1):
        for pileupcolumn in samfile.pileup(chrom, target, target+1, stepper='all', truncate=True, max_depth=10000):
                #print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                counts = dict(A=0, C=0, G=0, T=0);
                for base in pileupcolumn.pileups:
                    # .is_del -> the base is a deletion?
                    # .is_refskip -> the base is a N in the CIGAR string ?
                    if not base.is_del and not base.is_refskip:
                        counts[base.alignment.query_sequence[base.query_position]] += 1
                for base in counts.keys():
                    print("{}\t{}\t{}\t{}".format(pileupcolumn.pos, pileupcolumn.n, base, counts[base]))
        samfile.close()
