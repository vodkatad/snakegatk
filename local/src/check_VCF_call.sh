set -e
set -u

if [ $# -ne 4 ]; then
	echo "Usage: $0 genome.fa file.vcf[.gz] ALN output_basename" >&2
	exit 6
fi

REF=$1
VCF=$2
ALN=$3
NAME=$4 #$(basename $ALN)

BIN_DIR=$(realpath $(dirname $0))

FLANKING=20 # flanking bases
PERC_MATCH=90
PARALLEL=6

if [ -e "$NAME.to-blat.fa" ]; then
	echo "ERROR: file $NAME.to-blat.fa already exists. Exit. (FEATURE TO BE REMOVED)" >&2
	exit 74
fi

############# function ##########
function get_VCF_call_key {
# $1: VCF
	local VCF=$1
	zcat -f $VCF | \
		awk 'BEGIN{FS="\t"; OFS="\t"}
			!/^#/ {
				chrom = $1
				pos = $2
				REF = $4
				ALT = $5
				row = $0

				key = chrom":"pos"/"REF"/"ALT

				print key, row
			}' | \
		sort -T. -S2G
}

function extract_VCF_call_coords {
# $1: VCF
# $2: flanking bases
	local VCF=$1
	local FLANK=$2
	zcat -f $VCF | \
		awk -v FL=$FLANK 'BEGIN{FS="\t"; OFS="\t"}
			!/^#/ {
				chrom = $1
				pos = $2
				REF = $4
				ALT = $5

				L = length(REF)
				start = pos-FL
				end = pos+(L-1)+FL

				key = chrom":"pos"/"REF"/"ALT
				print chrom, start, end, key
			}'
}

function get_context {
	local REF=$1
	local POS=$2
	python3 $BIN_DIR/get_subseq-ordered.py $REF $POS
}

function inject_variation {
# $1: flanking bases
# $2: FASTA
	local FLANK=$1
	local FASTA=$2
	cat $FASTA | \
		paste - - | \
		awk -v FL=$FLANK '{
				name = $1
				seq = $2

				split(name, a, "/")
				REF = a[2]
				ALT = a[3]
				L = length(REF)

				before = substr(seq, 1, FL)
				after = substr(seq, FL+L+1)
				new_seq = before ALT after

				print name
				print new_seq
			}'
}

function fastq2fasta {
# $1: FASTQ[.gz]
	zcat -f $1 | \
		awk 'NR%4==1{ print ">"substr($1, 2)} NR%4==2{print $0}'
}

function do_blat {
# $1: variations sequences
# $2: FASTQ reads
	local VAR_SEQ=$1
	local READS=$2
#	fastq2fasta $READS | \
		blat -noHead -minScore=20 -tileSize=6 -out=pslx $VAR_SEQ $READS stdout
}

function filter_blat {
# $1: blat pslx output 
# $2: percentage of match
	local BLAT_OUT=$1
	local PERC=$2

# XXX forse troppo stringente con MM==0 e indel==0....
	cat $BLAT_OUT | \
		awk -v P=$PERC 'BEGIN{FS="\t"; OFS="\t"}
			{
			score = $1
			mismatch = $2
			ins_bases_in_read = $6
			ins_bases_in_variation = $8
			read_name = $10
			Lquery = $11
			Qstart = $12+1
			Qend = $13
			Lmatch1 = $13-$12
			variation_name = $14
			Lsubject = $15
			Sstart = $16+1
			Send = $17
			Lmatch2 = $17-$16

			if ((mismatch+ins_bases_in_read+ins_bases_in_variation==0) && (score>=(P/100)*Lsubject)) print
			}'
}

function get_best_match {
# $1: pslx file
	local PSLX=$1
	sort -T. -S2G -k10,10 -k1,1rn $PSLX | \
		awk 'BEGIN{FS="\t"}
			{
				if($10!=query){
					print
					query = $10
					score = $1
					}
				else if($1==score){
					print
					}
			}'
}

function get_overlapping_reads {
# $1: ALN 
# $2: coords (BED file)
	local ALN=$1
	local COORDS=$2
#	samtools view -@10 -L $COORDS $ALN | \
#		cut -f1,10 | \
#		awk '{print ">"$1"@"NR; print $2}'
	cut -f1 $COORDS | \
		sort -u | \
		parallel --tmpdir . --colsep ' ' -j $PARALLEL "samtools view -@10 -u $ALN {1} | samtools view -@10 -L $COORDS - | cut -f1,10 | awk '{print \">\"\$1\"@\"NR; print \$2}'" :::
}
#################################

TMP_FILE=$(mktemp -p .)
# VCF calls to key (key + row)
echo get_VCF_call_key ... >&2
get_VCF_call_key $VCF > $TMP_FILE.key

# get input for 'get_context'
echo extract_VCF_call_coords ... >&2
extract_VCF_call_coords $VCF $FLANKING > $TMP_FILE.coords

# get VCF call context
echo get_context ... >&2
#head $TMP_FILE.coords | \
cat $TMP_FILE.coords | \
	get_context $REF - | \
	inject_variation $FLANKING - > $TMP_FILE.to-blat.fa

# get reads overlapping context regions
echo get_overlapping_reads ... >&2
get_overlapping_reads $ALN $TMP_FILE.coords > $TMP_FILE.reads.fa

## blat & filter blat output
echo blat ... >&2
do_blat $TMP_FILE.to-blat.fa $TMP_FILE.reads.fa | \
	filter_blat - $PERC_MATCH | \
	get_best_match - > $TMP_FILE.OK.pslx

## summarize call failed and not
(
# print VCF header
zcat -f $VCF | \
	grep '^#'
# print successfull calls
join -t '	' <(cut -f14 $TMP_FILE.OK.pslx | sort -u) $TMP_FILE.key | \
	cut -f2- | \
	sort -k1,1V -k2,2n
) | \
	gzip > $NAME.OK.vcf.gz

# stats
zcat -f $VCF | \
	grep -v '^#' | \
	awk 'END{print "Original calls:", NR}' > $NAME.VCF-call.stats
zcat -f $NAME.OK.vcf | \
	grep -v '^#' | \
	awk 'END{print "Checked calls:", NR}' >> $NAME.VCF-call.stats

# clean
rm $TMP_FILE*
