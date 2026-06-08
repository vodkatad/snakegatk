set -e
set -u

MIN_q=20 # for mpileup
MIN_Q=20 # for mpileup
PARALLEL=10 #10
MIN_DEPTH=10 # XXX TODO
MIN_SUPP=1 # XXX TODO

export MIN_q MIN_Q

if [ $# -ne 4 ]; then
	echo "Usage: $0 ref.fa VCF ALN DEP_SOURCE_DIR" >&2
	echo "Do mpileup on SNV only in the VCF, and then call variation in naive mode." >&2
	echo "Chromosomes with '_' will be discarded." >&2
	exit 6
fi

REF=$1
VCF=$2
ALN=$3
SDIR=$4 

TMPFILE=$(mktemp -p .)
KEY_LIST="$TMPFILE.key-list"
BED="$TMPFILE.bed"
MPILEUP="$TMPFILE.mpileup"
COUNTS="$TMPFILE.counts"

export ALN BED MPILEUP
############ function ###############
function get_SNVs {
# $1: VCF
	zcat -f $1 | \
		awk '(length($4)==1) && (length($5)==1) && ($1 !~ /_/)' | \
		${SDIR}/VCF_key.sh - | \
		awk -v F=$1 'BEGIN{OFS="\t"}{print $1, F}'
}

function key_list2bed {
# $1: key-list file
	awk 'BEGIN{FS="\t"; OFS="\t"}
		{
			key = $1
			n = split(key, a, ":")
			if(n!=4) {
				print "ERROR: wrong key:", $1 > "/dev/stderr"
				exit 1
			}
			chrom = a[1]
			pos = a[2]
			REF = a[3]
			ALT = a[4]

			print chrom, pos-1, pos, key
		}' $1
}

function get_common_chrs_list {
# $1: ALN
# $2: BED
	local ALN=$1
	local BED=$2
	join <(
		samtools view -H $ALN | \
			awk 'BEGIN{FS="\t"} /^@SQ/{split($2,a,":"); print a[2]}' | \
			sort
		) \
		<(
			cut -f1 $BED | \
				sort -T . -S2G -u
		)
}


function mpileup_parallel {
# $1: REF
# $2: ALN
# $3: BED
# $4: output basename
# $5: chrom
# output in <OUT>.<CHROM>
	local REF=$1
	local ALN=$2
	local BED=$3
	local OUT=$4
	local CHROM=$5
	samtools mpileup -f $REF -q $MIN_q -Q $MIN_Q -l $BED -r $CHROM $ALN > $OUT.$CHROM
}
export -f mpileup_parallel

########################################

# get SNV from VCF as key
get_SNVs $VCF | \
	sort -T . -S2G -u > $KEY_LIST

# bed file of SNVs
key_list2bed $KEY_LIST > $BED

# BED files by chrom
cut -f-3 $BED | \
	awk -v OUT=$BED '{
		C = $1
		print > OUT"."C
	}'

# mpileup on SNV positions
get_common_chrs_list $ALN $BED | \
	parallel --tmpdir . --colsep ' ' -j $PARALLEL "mpileup_parallel $REF $ALN $BED.{1} $MPILEUP {1}"

# mpileup to counts
cat $MPILEUP.* | \
	${SDIR}/mpileup2counts.py - | \
	sort -T . -S2G > $COUNTS

# match key and counts,
# then annotate SNVs (WT or MUT),
# finally calculate overlap
join -a1 -e '?' -o auto -t '	' <(cut -f4 $BED | \
			awk 'BEGIN{FS="\t"; OFS="\t"}
			{
				key = $1
				split(key, a, ":")
				coord = a[1]":"a[2]
				print coord, key
			}' $KEY_LIST | \
			sort -T . -S2G) $COUNTS | \
	sort -T . -S2G -V | \
	awk -v MIN_DEPTH=$MIN_DEPTH -v MIN_SUPP=$MIN_SUPP 'BEGIN{
		FS = "\t"
		OFS = "\t"
		bases["A"] = 4
		bases["C"] = 5
		bases["G"] = 6
		bases["T"] = 7
		print "#key", "VAF", "N_MUT", "N_WT", "Tot"
	}
	{
##chr1:14574      chr1:14574:A:G  A       16      0       38      0       1       0       0       0       .       .
		coord = $1
		key = $2
		split(key, a, ":")
		REF = a[3]
		ALT = a[4]

		if($3=="?") {
			print key, "NA", "NA", "NA", "NA"
			next
		}

		depth = 0
		for(i=4;i<=11;i++) depth+=$i

		if(depth<MIN_DEPTH) {
			print key, "NA", "NA", "NA", "NA"
			next
		}

		support = $(bases[ALT])+$(bases[ALT]+4)
		VAF = 100*support/depth

		print key, VAF, support, depth-support, depth
	}'
###################################

# clean
rm $TMPFILE $KEY_LIST $BED $BED.* $COUNTS $MPILEUP.*
