#!/bin/bash
set -e
set -u

if [ $# -ne 1 ]; then
	echo "Usage: $0 VCF[.gz]|-" >&2
	echo "Prepend key: 'chr:pos:REF:ALT'" >&2
	echo ">>> VCF row data must have at least 8 columns <<<" >&2
	exit 87
fi

function prepend_key {
        zcat -f $1 | \
                grep -v '^#' | \
                awk 'BEGIN{FS="\t"; OFS="\t"}
			NF<8{
				print "ERROR: VCF row data must have at least 8 columns, not", NF, "at row:", NR, "Exit." > "/dev/stderr"
				exit 1
			}
			{
				key=$1":"$2":"$4":"$5
				print key, $0
			}'
}

prepend_key $1
