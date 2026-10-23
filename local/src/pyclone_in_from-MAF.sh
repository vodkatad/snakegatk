set -e
set -u

if [ $# -ne 3 ]; then
	echo "Usage: $0 SAMPLE sequenza_segments.txt file.maf" >&2
	echo "Filter MAF file by <SAMPLE> and return the input for pyclone." >&2
	exit 11
fi
MAX_CN=5
LOH="no" # "no": filter LOH, otherwise do not filter

SAMPLE=$1
SEQUENZA=$2
MAF=$3

############ functions ############
function sequenza2bed {
# $1: sequenza file
	awk 'BEGIN{FS="\t"; OFS="\t"}
		NR>1{
			chrom = $1
			start0 = $2-1
			end = $3
			CN_tot = $10
			A = $11
			B = $12

			if((A!="NA") && (B!="NA") && (CN_tot!="NA")) print chrom, start0, end, A, B, CN_tot
		}' $1
}

function maf2bed {
# $1: maf file
# $2: sample
	awk -v SMP=$2 'BEGIN{FS="\t"; OFS="\t"}
		NR>1 {
			chrom = $5
			start = $6
			start0 = start-1
			end = $7
			#base0 = $11
			#base1 = $12
			#base2 = $13
			sample = $16
			REF = $41
			ALT = $42

			mut_id = chrom"_"start

			if(sample==SMP) print chrom, start0, end, mut_id, sample, REF, ALT
		}' $1
}

###################################

bedtools intersect -wa -wb -a <(maf2bed $MAF $SAMPLE) -b <(sequenza2bed $SEQUENZA) | \
	awk -v MAX_CN=$MAX_CN -v LOH=$LOH 'BEGIN{FS="\t"; OFS="\t"
		print "mutation_id", \
			"sample_id", \
			"ref_counts", \
			"alt_counts", \
			"major_cn", \
			"minor_cn", \
			"normal_cn", \
			"tumour_content", \
			"error_rate"
		}
		{
			mut_id = $4
			sample = $5
			REF = $6
			ALT = $7
			A = $11
			B = $12
			CN_tot = $13

			normal_CN = 2
			if(mut_id ~ /chrY/) normal_CN = 1

			if((LOH=="no") && (B<=0)) next # no LOH
			if((CN_tot<=MAX_CN)) print mut_id, sample, REF, ALT, A, B, normal_CN, "1", "0.001"
		}'
