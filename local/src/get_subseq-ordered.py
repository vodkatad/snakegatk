#!/usr/bin/python
import sys
from Bio.SeqIO.FastaIO import FastaIterator as FastaIterator

if len(sys.argv) != 1+2:
	sys.exit("""
	Usage: %s file.fna|- file.pos|-
	Subsequence of the sequences listed in file.pos.
	'file.pos' structure:
		seq_name start end [newSeqName]
		...
	If start > end, get the revcomp.

	Position order is kept.

	Output in STDOUT
	""" % sys.argv[0])

if sys.argv[1]=='-' and sys.argv[2]=='-':
	sys.exit("Too many STDIN input!")

def read_FASTA_file(f):
    "legge un file FASTA"
    if f == '-':
        fd = sys.stdin
    elif hasattr(f, "read"):
        fd = file
    else:
        fd = open(f)
    # iterator = FASTA.FastaReader(fd)
    iterator = FastaIterator(fd)
    for i in iterator:
        yield i
    fd.close()


fasta = read_FASTA_file(sys.argv[1])
file_pos = sys.argv[2]=='-' and sys.stdin or open(sys.argv[2])

fasta_dict = {}
for i in fasta:
	fasta_dict[i.name] = i


#positions = {}
positions = []
for x in file_pos:
	data = x.strip('\n').split('\t')
	note = (len(data)==4) and data[3] or None
	positions.append([data[0],data[1:3], note])

#for seq in fasta:
for pos in positions:
	name = pos[0]
	start,end = int(pos[1][0]),int(pos[1][1])
	note = pos[2]
	(start,end,strand) = (start<=end) and (start,end,1) or (end,start,-1)
	if name in fasta_dict:
		my_seq = fasta_dict[name]
		if start < 1 or end < 1:
			print("Invalid start or end: %s %d %d." % (name, start, end), file=sys.stderr)
			continue
		elif end > len(my_seq.seq):
			print("WARNING: The 'end' is over sequence length: %s (len=%d) %d." % (name, len(my_seq.seq), end), file=sys.stderr)
		if note:
			print(">"+note)
		else:
			print(">%s:%d-%d" % (name, start, end))
		if strand == 1:
			print(my_seq.seq[start-1:end])
		else:
			print(my_seq.seq[start-1:end].reverse_complement())
	else:
		print("WARNING: sequence '%s' not found in FASTA file. Skipping." % name, file=sys.stderr)
