#!/usr/bin/env python3

import argparse
import fileinput
import re
import sys
from collections import Counter

#COLUMNS = ['A','C','G','T','a','c','g','t']
patt_ins = re.compile("\\+[0-9]+")
patt_del = re.compile("\\-[0-9]+")
#patt_start_read = re.compile("\\^.") # XXX
#patt_end_read = re.compile("[$]") # XXX
#patt_ref_skip = re.compile("[<>]") # XXX
#patt_del_ref = re.compile("[*#]") # XXX
#patt_ref_bases = re.compile("[.,]") # XXX
#patt_nonRef_bases = re.compile("[ACGTacgtNn]") # XXX

DEBUG = False

def filter_and_count_iter(row_iter):
	for row in row_iter:
		chrom, pos, ref_base, depth, bases, qualities = row.strip().split('\t')

		if ref_base.upper() == 'N':
			continue

		#insertions = []
		#deletions = []
		insertions = Counter()
		deletions = Counter()

# XXX		bases_purged = purge_pileup(ref_base, bases)
		if DEBUG:	print("Here", bases, file=sys.stderr)
		L_bases = len(bases)
		ref_base_UPPER = ref_base.upper()
		ref_base_LOWER = ref_base.lower()
		#bases_purged = []
		bases_purged = {}
		bases_purged["A"] = 0
		bases_purged["C"] = 0
		bases_purged["G"] = 0
		bases_purged["T"] = 0
		bases_purged["a"] = 0
		bases_purged["c"] = 0
		bases_purged["g"] = 0
		bases_purged["t"] = 0
		#bases_purged = Counter()
		idx = 0
		while idx<L_bases:
			b = bases[idx]
			# patt_start_read
			if b == '^':
				if DEBUG: print(">", bases[idx:], "patt_start_read", file=sys.stderr)
				idx += 2
				if DEBUG: print("now:", bases[idx:], file=sys.stderr)
				continue
			# patt_end_read + patt_ref_skip + patt_del_ref
			if b in ['$', '<', '>', '*', '#']:
				if DEBUG: print(">", bases[idx:], "patt_end_read + patt_ref_skip + patt_del_ref", file=sys.stderr)
				idx += 1
				if DEBUG: print("now:", bases[idx:], file=sys.stderr)
				continue
			# patt_ins / patt_del
			if b in ['+', '-']:
				# ins ?
				m = patt_ins.match(bases, idx)
				if m:
					if DEBUG: print(">", bases[idx:], "patt_ins", file=sys.stderr)
					L = m.end()-idx
					N = int(bases[idx+1:idx+L])
					INS_bases = bases[idx+L:idx+L+N]
					#insertions.append(INS_bases)
					insertions.update([INS_bases])
					# removes also INS bases
					idx += L+N
					if DEBUG: print("now:", bases[idx:], "ins:", insertions, file=sys.stderr)
					continue
				# del ?
				m = patt_del.match(bases, idx)
				if m:
					if DEBUG: print(">", bases[idx:], "patt_del", file=sys.stderr)
					L = m.end()-idx
					N = int(bases[idx+1:idx+L])
					DEL_bases = bases[idx+L:idx+L+N]
					#deletions.append(DEL_bases)
					deletions.update([DEL_bases])
					# removes also DEL bases
					idx += L+N
					if DEBUG: print("now:", bases[idx:], file=sys.stderr)
					continue
				sys.exit("ERROR: something went wrong here: %s at this point: %s" % (bases, bases[idx:]))
			# patt_ref_bases + patt_nonRef_bases
			if b in ['A','C', 'G', 'T', 'a', 'c', 'g', 't', 'N', 'n', '.', ',']: # XXX aggiunge 0.3 sec su 1M
				if DEBUG: print(">", bases[idx:], "patt_ref_bases and patt_nonRef_bases", file=sys.stderr)
				if b == '.':
					#bases_purged.append(ref_base_UPPER)
					#bases_purged.update(ref_base_UPPER)
					bases_purged[ref_base_UPPER] += 1
				elif b == ',':
					#bases_purged.append(ref_base_LOWER)
					#bases_purged.update(ref_base_LOWER)
					bases_purged[ref_base_LOWER] += 1
				else:
					#bases_purged.append(b)
					#bases_purged.update(b)
					try:
						bases_purged[b] += 1
					except KeyError:
 						pass
				idx += 1
				if DEBUG: print("now:", bases[idx:], file=sys.stderr)
				continue

			sys.exit("ERROR: something went wrong here: %s at this point: %s" % (bases, bases[idx:]))

		counts = [
			"%s:%s" % (chrom, pos),
			ref_base_UPPER,
			bases_purged["A"],
			bases_purged["C"],
			bases_purged["G"],
			bases_purged["T"],
			bases_purged["a"],
			bases_purged["c"],
			bases_purged["g"],
			bases_purged["t"],
			#insertions and ','.join(insertions) or '.',
			#deletions and ','.join(deletions) or '.',
			insertions and ','.join([ ':'.join((map(str,i))) for i in insertions.items() ]) or '.',
			deletions and ','.join([ ':'.join((map(str,i))) for i in deletions.items() ]) or '.'
		]

		print(*counts, sep='\t')

def main():
    desc = """Parse samtools mpileup output to generate an 8-columns count file.

Return: coord, ref_base, Afwd_count, Cfwd_count, Gfwd_count, Tfwd_count, Arev_count, Crev_count, Grev_count, Trev_count, insertions_list, deletions_list

Skip 'N' positions or 'N' bases.
"""

    parser = argparse.ArgumentParser(
        description = desc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    parser.add_argument('mpileup_file', help = "samtools mpileup input file "
            "or '-' for stdin.")

    args = parser.parse_args()
    print("WARNING: filter_and_count_mpileup: 8-fields count file!", file=sys.stderr)
    row_iter = fileinput.input(args.mpileup_file)
    filter_and_count_iter(row_iter)

if __name__ == "__main__":
    main()
