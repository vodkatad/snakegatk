#!/usr/bin/env python3

import argparse
import sys
import gzip

# return next bed entry as vector of chr, b, e, CN
def next_bed_entry(bedcn):
    line = bedcn.readline()
    if line != '':
        line.rstrip('\n')
        entry = line.split('\t')
        return([entry[0], int(entry[1])-1, int(entry[2]), int(entry[5])])
    return(None)

# return next germline entry with more than min_reads support as vector of chr, b, e, VAF 
def next_germ_entry(germ, min_reads=5, verbose=False):
    line = germ.readline()
    while line != '':
        line.rstrip('\n')
        entry = line.split('\t')
        if verbose:
            print("germ {}".format(entry[0]), file=sys.stderr)
        coords = entry[0].split(':')
        if int(entry[2]) >= min_reads:
            return([coords[0], int(coords[1])-1, int(coords[1]), float(entry[1])])
        else:
            line = germ.readline()
    return(None)


def chr_to_n(chrs):
    try: 
        chr_n = int(chrs)
        return chr_n
    except ValueError:
        # chrX/Y
        if chrs == 'X':
            return 23
        else: 
            return 24

# return the first entry that overlap segment and True if the list of germlines is finished
# false otherwise, with none if the current mut is preceding this segment (need to get next one at our next invocation)
# abstracting this logic (1-lookeahead) from the while makes it easier to reason about
# NOOOO
# Since mutations are punctual we have an easier task:
# the next one is inside this segment -> work on it, then advance on muts
# the next one is before this segment -> advance on muts
# the next one is after this segment -> advance on segments
# so we return: the mut + True if we need to process it then advance
# None + True if we need to advance on muts
# None + False if we need to advance on segments
# First entry of the tuple: True if muts are ended, False otherwise
last_entry = None
def get_next_overlapping_germ(segment, germ, verbose):
    if verbose:
        print("get for {}".format(segment[1]), file=sys.stderr)
    global last_entry
    if verbose:
        print("get2 for {}".format(last_entry), file=sys.stderr)
    if last_entry == None:
        last_entry = next_germ_entry(germ, verbose=verbose)
        if last_entry == None:
            return(True, None, False)

    result = None 
    advance = True
    # overlap check  a0 <= b1 && b0 <= a1;
    # https://fgiesen.wordpress.com/2011/10/16/checking-for-interval-overlap/
    # < and not <= for end excluded
    if verbose:
        print('evaluating overlap {} {}'.format(last_entry, segment), file=sys.stderr)
    if last_entry != None and last_entry[0] == segment[0]: # same chr
        if last_entry[1] < segment[2] and segment[1] < last_entry[2]: # overlap
            if verbose:
                print('ov', file=sys.stderr)
            result = last_entry
            last_entry = None # we want to advance next time
            advance = True
        elif last_entry[1] < segment[1]: # same chr before the segment, we need to advance:
            result = None
            advance = True
            last_entry = None
        else:
            # otherwise we are after this segment, keep this mut
            result = None
            advance = False
    # different chr we keep this mut or not depending on chr order
    else:
        chr_germ = last_entry[0][3:]
        chr_segm = segment[0][3:]
        chr_germ_n = chr_to_n(chr_germ)
        chr_segm_n = chr_to_n(chr_segm)
        if (chr_germ_n < chr_segm_n): # go on segm
            result = None
            advance = True
            last_entry = None
        else: # go on germ
            result = None
            advance = False
    
    if verbose:
        print('returning from get next overlap {} {}'.format(last_entry, result), file=sys.stderr)
    return(False, result, advance)

# evaluate if this germline mutation VAF is compatible with the reported CN
# TODO binomial?
def manage_overlap(gentry, segment, epsilon=0.05, verbose=False):
    cn = segment[3]
    vaf = gentry[3]
    i = 1
    for i in range(1, cn+1):
        comp = i/cn
        if vaf >= (comp-epsilon)*100 and vaf <= (comp+epsilon)*100:
            if verbose:
                print('evaluating compatibility {} {} ok'.format(cn, vaf), file=sys.stderr)
            return True
    if verbose:
        print('evaluating compatibility {} {} fail'.format(cn, vaf), file=sys.stderr)
    return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Iterate over cnvkit segments and hp germline muts to evaluate quality of integer CN inference/ploidy")

    parser.add_argument('--cn', '-c', dest='cn', action='store', help='CNVkit cn calls - has to be sorted!', required=True)
    parser.add_argument('--germ', '-g', dest='germ', action='store', help='tsv file with hypothetical germline muts - id vaf support tot, has to be sorted! ', required=True)
    parser.add_argument('--verbose', '-v', action='store_true',  help='verbose execution')

    args = parser.parse_args()
    if (args.verbose):
        print('cn\t{}'.format(args.cn), file=sys.stderr)
        print('germ\t{}'.format(args.germ), file=sys.stderr)
    
    with open(args.cn, 'r') as cn:
        cn.readline() # skip header
        with gzip.open(args.germ, 'rt') as germ:
            germ.readline() # skip header
            done = False
            mdone = False
            while not done:
                segment = next_bed_entry(cn)
                segment_ok = 0
                segment_ov = 0
                if segment is None:
                    done = True
                else:
                    while not mdone:
                        mdone, gentry, goon = get_next_overlapping_germ(segment, germ, args.verbose)
                        if gentry is not None:
                            if manage_overlap(gentry, segment, verbose=args.verbose):
                                segment_ok =  segment_ok + 1
                                # future: weight score on segments lengths 
                            segment_ov = segment_ov + 1
                        if not goon:
                            break
                        
                    print('{}\t{}\t{}\t{}\t{}'.format(segment[0],segment[1],segment[2], segment_ok, segment_ov))
            