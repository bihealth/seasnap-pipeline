# ---------------------------------------------------------------------------
# Function require_var(VARNAME, PARAMETER_NAME)
# ---------------------------------------------------------------------------
#
# Exit with return code of 1 if the variable with the name VARNAME is not
# set.  The variable is expected to come from a parameter that is given
# for documentation/fixing purposes.

require_var()
{
    VARNAME=$1
    PARAM=$2
    test -z "${!VARNAME}" || return 0
    echo "ERROR: Missing parameter ${PARAM} ${VARNAME}" >&2
    echo "" >&2
    print_help
    exit 1
}

# ---------------------------------------------------------------------------
# Function print_header(PROG, VAR1, ...)
# ---------------------------------------------------------------------------
#
# Print program output header with the given program name and variables.

print_header()
{
    echo "Program: $1" >&2
    echo "" >&2
    shift

    echo "Parameters" >&2
    echo "" >&2

    for var in "$@"; do
        echo -e "  ${var}\t${!var}" >&2
    done

    echo "" >&2
}

create_read_file_list()
{
    FASTQ=($@)
    TMP_LIST_FASTQ=$(mktemp)
    for FILE in ${FASTQ[@]}; do
	echo $FILE >> $TMP_LIST_FASTQ
    done
    echo $TMP_LIST_FASTQ
}


get_files()
{
    INPUT_DIR=$1
    PATTERN=$2
    for f in ${INPUT_DIR}/${PATTERN}; do
	[[ -f "$f" ]] || { echo "No files with pattern ${PATTERN} found in ${INPUT_DIR}" >&2 ; exit 1; } 
    done

    files=$(ls ${INPUT_DIR}/${PATTERN})
    echo -n $files 
}

# ---------------------------------------------------------------------------
# Function get_reads_dir(STUDY_DIR, READ_SET, READ_TAG)
# ---------------------------------------------------------------------------
#
# Build the read base dir from study directory, read set name, and read
# tag name, write to stdout.  Call with $(get_read_base_dir ...)

get_reads_dir()
{
    echo "$1/samples/$2/fastq/$3"
}

# ---------------------------------------------------------------------------
# Function require_file(PATH, COMMENT)
# ---------------------------------------------------------------------------
#
# Require PATH to point ot a file or exit with code 1 otherwise.

require_file()
{
    test -f ${1} || { echo "Not a file ${1} ($2)" >&2; exit 1; }
}

# ---------------------------------------------------------------------------
# Function get_left_reads(READ_BASE_DIR)
# ---------------------------------------------------------------------------
#
# Return all left read sets, given the read base directory, by writing to
# stdout.  Exit with return code 1 if there is no such file.

get_left_reads()
{
    for x in $1/*_R1_*.f*q.gz; do
        require_file $x "left reads"
        echo $x
    done
}

# ---------------------------------------------------------------------------
# Function get_right_reads(LEFT_FILE1, ...)
# ---------------------------------------------------------------------------
#
# Return file name(s) for right reads, given the left one(s), result is
# written to stdout.

get_right_reads()
{
    echo $(echo $* | sed 's/_R1_/_R2_/g')
}

# ---------------------------------------------------------------------------
# Function get_reports_dir(STUDY_DIR, READ_SET, READ_TAG, REPORT_NAME)
# ---------------------------------------------------------------------------
#
# Build directory for the reports, given study base dir, read set name, read
# tag, and report name.

get_reports_dir()
{
    echo "$1/samples/$2/reports/$4/$3"
}

# ---------------------------------------------------------------------------
# Function get_vcf_dir(STUDY_DIR, DONOR)
# ---------------------------------------------------------------------------
#
# Build directory for the VCF given study base dir.

get_vcf_dir()
{
    echo "$1/analysis/${DONOR}/vcf"
}

# ---------------------------------------------------------------------------
# Function get_analysis_log_dir(STUDY_DIR, DONOR)
# ---------------------------------------------------------------------------
#
# Build directory for the VCF given study base dir.

get_analysis_log_dir()
{
    echo "$1/analysis/${DONOR}/log"
}

# ---------------------------------------------------------------------------
# Function get_analysis_report_dir(STUDY_DIR, DONOR)
# ---------------------------------------------------------------------------
#
# Build directory for the VCF given study base dir.

get_analysis_report_dir()
{
    echo "$1/analysis/${DONOR}/report"
}

# ---------------------------------------------------------------------------
# Function get_bam_dir(STUDY_DIR, READ_SET)
# ---------------------------------------------------------------------------
#
# Build directory for the bam, given study base dir, read set name, read
# tag, and report name.

get_bam_dir()
{
    echo "$1/samples/$2/bam"
}

# ---------------------------------------------------------------------------
# Function get_data_dir(BASE_DIR)
# ---------------------------------------------------------------------------
#
# Given the base dir to this script, return path to data directory.

get_data_dir()
{
    echo "/vol/fs01/data/kajetan/my_ngs_stock_data/"
}


# ---------------------------------------------------------------------------
# Function get_git_name()
# ---------------------------------------------------------------------------
#
# returns the git name of the current working directory

function get_git_name {
   name=$(basename $(git rev-parse --show-toplevel) 2> /dev/null) || return
   echo "${name}"
}


# ---------------------------------------------------------------------------
# Function get_git_branch()
# ---------------------------------------------------------------------------
#
# returns the git branch of the current working directory

function get_git_branch {
   ref=$(git symbolic-ref HEAD 2> /dev/null) || return
   echo "${ref#refs/heads/}"
}

# ---------------------------------------------------------------------------
# Function get_git_hash()
# ---------------------------------------------------------------------------
#
# returns the git hash of the current working directory

function get_git_hash {
    hash=$(git rev-parse HEAD 2> /dev/null) || return
    echo $hash
}

# ---------------------------------------------------------------------------
# Function join_by()
# ---------------------------------------------------------------------------
#
# returns the values of an array joined into a string by the first argument
# example: 
#   input: ids=("a" "b" "c" "d")
#   call:  join_by , ${ids[@]}
# output:  a,b,c,d
function join_by {
    local IFS="$1"; shift; echo -n "$*";
}


