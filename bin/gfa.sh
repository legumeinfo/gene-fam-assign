#!/usr/bin/env bash

# shellcheck disable=SC2154

set -o errexit -o errtrace -o nounset -o pipefail -o posix

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

# to help assign a heredoc value to a variable. The internal line return is intentional.
define(){ o=; while IFS=$'\n' read -r a; do o="$o$a"'
'; done; eval "$1=\$o"; }

define HELP_DOC <<'EOS'
NAME
  hmmsearch-to-gfa.sh  -- Given a database of gene family HMMs (a concatenation of HMMs produced by hmmbuild, then 
    compressed with hmmpress), manage the search of a collection of fasta files against the database and
    the subsequent parsing into tabular "gene family assignment" files.

SYNOPSIS
  hmmsearch-to-gfa.sh -l filepath_list

  Required:
           -l  - path to list of file paths to compressed (.gz) fasta files
           -c  - path to the config file

  Options: -h  - help
           -t  - threads (or may be set in config file; see comments below)

  VARIABLES set in config file:
    hmmdb    - Name of hmmer HMM data files (sans .h suffixes). HMM database should be in the data directory
    evalue   - E-value threshold for hmmsearch [1e-10]
    data_dir - Directory for HMM database [data]; also a reasonable place to put the list of fasta filepaths
    work_dir - Directory where all work will be done
    threads  - Value for hmmsearch --cpu. [2] Recommended to be small (2), as most of the parallelism comes
                 from running on multiple query files, determined by available threads.

AUTHORS
    Steven Cannon <steven.cannon@usda.gov>, Andrew Farmer <adf@ncgr.org>
EOS

if [ "$#" -eq 0 ]; then
  echo >&2 "$HELP_DOC" && exit 0;
fi

# Command-line interpreter
filepath_list="null"
config="null"

while getopts "l:c:t:h" opt
do
  case $opt in
    l) filepath_list=$OPTARG; echo "filepath_list: $config" ;;
    c) config=$OPTARG; echo "config: $config" ;;
    t) threads=$OPTARG; echo "threads: $threads" ;;
    h) echo >&2 "$HELP_DOC" && exit 0 ;;
    *) echo >&2 echo "$HELP_DOC" && exit 1 ;;
  esac
done

if [ "$filepath_list" == "null" ]; then
  printf "\nPlease provide the path to a list of (compressed) fasta files to be searched: -l filepath_list\n" >&2
  exit 1;
fi

if [ "$config" == "null" ]; then
  printf "\nPlease provide the path to the config file, e.g. -c config/gfa.conf\n" >&2
  exit 1;
fi

# Add shell variables from config file
# shellcheck source=/dev/null
. "${config}"

# Check for existence of third-party executables
missing_req=0
dependencies='hmmsearch'
for program in $dependencies; do
  if ! type "$program" &> /dev/null; then
    echo "Warning: executable $program is not on your PATH."
    missing_req=$((missing_req+1))
  fi
done

# Check that the bin directory is in the PATH
if ! type strip_spliceform.pl &> /dev/null; then
  printf "\nPlease add the bin directory to your PATH and try again. Try the following:\n"
  printf "\n  PATH=%s/bin:\%s\n\n" "$PWD" "$PATH"
  exit 1; 
fi

echo
echo "== Copy files into work directory, and uncompress."

WORK="$work_dir"
LIST=$(realpath "$filepath_list")
DATA=$(realpath "$data_dir")
EVALUE="$evalue"
THREADS="$threads"

NPROC=$( ( command -v nproc > /dev/null && nproc ) || getconf _NPROCESSORS_ONLN)

export NPROC=${NPROC:-1}
if [[ $NPROC -gt 1 ]]; then NPROC_PER_THREAD=$(( NPROC / THREADS )); else NPROC_PER_THREAD=1; fi
if [[ $NPROC_PER_THREAD -lt 1 ]]; then NPROC_PER_THREAD=1; fi
echo "NPROC: $NPROC"
echo "THREADS: $THREADS"
echo "NPROC_PER_THREAD: $NPROC_PER_THREAD"
echo "WORK_DIR: $WORK"

WD=$(realpath "$WORK")
mkdir -p "${WD}"
cd "${WD}" || exit

mkdir -p 00_fasta 01_hmmsearch 02_gfa

echo
echo "== Copy files into working directory, uncompressing them for access by hmmsearch"
echo
cat "$LIST" | while read -r filepath; do
  if [[ "$filepath" == \#* ]]; then
    break
  else 
    if [ -f "$filepath" ]; then
      if [[ "$filepath" =~ \.gz$ ]]; then
        file=$(basename "$filepath" .gz)
        if [ -f 00_fasta/"$file" ]; then
          echo "  Uncompressed file exists already at 00_fasta/$file"
        else
          echo "  Uncompressing input to 00_fasta/$file"
          gzip -dc "$filepath" > 00_fasta/"$file"
        fi
      else
        file=$(basename "$filepath")
        echo "  File appears to be uncompressed; use as-is: $file"
        cat "$filepath" > 00_fasta/"$file"
      fi
    else
      echo "  WARN: File $filepath was not found."
    fi
  fi
done

echo
echo "== Search fasta files against HMM database"
echo
for querypath in 00_fasta/*; do
  base=$(basename "$querypath" .faa)
  echo "hmmsearch -E $EVALUE --cpu $THREADS --tblout 01_hmmsearch/$base.hmmsearch.tbl -o /dev/null ${DATA}/${hmmdb} ${querypath}"
  hmmsearch -E "$EVALUE" --cpu "$THREADS" --tblout 01_hmmsearch/"$base".hmmsearch.tbl -o /dev/null "${DATA}"/"${hmmdb}" "${querypath}" &
  # allow to execute up to $NPROC_PER_THREAD in parallel
  if [[ $(jobs -r -p | wc -l) -ge ${NPROC_PER_THREAD} ]]; then wait -n; fi
done
wait

echo
echo "== Parse tabular hmmsearch data, picking the top match per query"
echo
for filepath in 01_hmmsearch/*.hmmsearch.tbl; do
  file=$(basename "$filepath" .hmmsearch.tbl)
  base=$(echo "$file" | perl -pe 's/(.+)\.\w+/$1/') # strips suffix like .protein_primary or .protein
  cat /dev/null > 02_gfa/"${base}"."${hmmdb}".gfa.tsv
  printf "#gene\tfamily\tprotein\te-value\tscore\n" > 02_gfa/"${base}"."${hmmdb}".gfa.tsv
  awk -v OFS="\t" '$1!~/^#/ {print $1, $3, $1, $5, $6}' "$filepath" |
    sort -k1,1 -k5nr,5nr | top_line.awk |
      strip_spliceform.pl >> 02_gfa/"${base}"."${hmmdb}".gfa.tsv &
   # allow to execute up to $NPROC in parallel
   if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
done
wait

cd "$OLDPWD" || exit

echo
echo "== Run completed. Look for results at $WORK/02_gfa/"
echo

exit 0
