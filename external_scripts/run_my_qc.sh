#!/bin/bash

# Script: run_my_fastqc.sh

# Enable "set -xv" for debugging purposes.
set -e
#set -xv

# ---------------------------------------------------------------------------
# Bootstrap script and helpers, get GIT infos
# ---------------------------------------------------------------------------

DIR="$(cd "$(dirname "$0")" && pwd)"
. ${DIR}/.helpers.sh

GIT_NAME="$(cd "$(dirname "$0")" && get_git_name)"
GIT_BRANCH="$(cd "$(dirname "$0")" && get_git_branch)"
GIT_HASH="$(cd "$(dirname "$0")" && get_git_hash)"

# ---------------------------------------------------------------------------
# Define the help
# ---------------------------------------------------------------------------

print_help()
{
    echo "Run quality control on a single sample" >&2
    echo "" >&2
    echo "Usage:   run_my_qc.sh -l FILE.star.log -f FASTQCDIR_ORG" >&2
    echo "             -b FILE.bam -g FILE.gtf -r REF.fasta -o OUTDIR" >&2
    echo "" >&2
    echo "  -h             Show this help" >&2
    echo "  -l STARLOG     The STAR log file" >&2
    echo "  -f FASTQCDIR   The directory with the fastqc reports of the original fastq files" >&2
    echo "  -b BAM         The input bam file" >&2
    echo "  -e             BAM contains single end reads" >&2
    echo "  -o OUTDIR      The folder to which qc data will be written" >&2
    echo "  -g GTF         The annotation gtf file" >&2
    echo "  -r REF         The genome reference fasta file" >&2
    echo "" >&2

}

# ---------------------------------------------------------------------------
# Parse the command line and check values
# ---------------------------------------------------------------------------

SINGLEEND=0

while getopts "l:f:b:g:r:o:he" opt; do
    case $opt in
        h)
            print_help
            exit 0
            ;;
	e)
	    SINGLEEND=1
	    ;;
	l)
	    STAR_LOG=${OPTARG}
	    ;;
	f)
	    FASTQC_ORIGINAL_DIR=${OPTARG}
	    ;;
	b)
	    BAM=${OPTARG}
	    ;;
	g)
	    GTF=${OPTARG}
	    ;;
	r)
	    REF=${OPTARG}
	    ;;
 	o)
	    OUTPUT_DIR=${OPTARG}
	    ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            ;;
    esac
done

require_var STAR_LOG -l
require_var FASTQC_ORIGINAL_DIR -f
require_var BAM -b
require_var GTF -g
require_var REF -r
require_var OUTPUT_DIR -o

mkdir -p $OUTPUT_DIR

SAMPLE=$(basename ${BAM})
SAMPLE=${SAMPLE%.bam}

print_header run_my_qc.sh SAMPLE STAR_LOG FASTQC_ORIGINAL_DIR BAM SINGLEEND GTF REF GIT_NAME GIT_BRANCH GIT_HASH

echo "Running QC" >&2
echo "" >&2

set -x
${DIR}/read_statistic_report.sh -l $STAR_LOG  -g $FASTQC_ORIGINAL_DIR -o ${OUTPUT_DIR}/read_alignment_report.tsv
if [ $SINGLEEND -eq 1 ]
then
    ${DIR}/run_RNA-SeQC.sh -i $BAM -e -t $GTF -r $REF -o ${OUTPUT_DIR}/
else
    ${DIR}/run_RNA-SeQC.sh -i $BAM -t $GTF -r $REF -o ${OUTPUT_DIR}/ 
fi
${DIR}/read_duplication.sh -i $BAM -o ${OUTPUT_DIR}/duplication/
