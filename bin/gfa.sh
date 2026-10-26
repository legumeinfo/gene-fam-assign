#!/usr/bin/env bash

# shellcheck disable=SC2154

set -o errexit -o errtrace -o nounset -o pipefail -o posix

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

# to help assign a heredoc value to a variable. The internal line return is intentional.
define(){ o=; while IFS=$'\n' read -r a; do o="$o$a"'
'; done; eval "$1=\$o"; }

define HELP_DOC <<'EOS'
NAME
  gfa.sh  -- Given either a database of gene family HMMs (a concatenation of HMMs produced by hmmbuild,
    then compressed with hmmpress), or a collection of family multi-fasta protein files representing 
    all gene families, with defline formats like the following ...
      >Legume.fam3.00001__glyma.Wm82.gnm4.ann1.Glyma.06G281100.1
    manage the search of a collection of fasta files against the database and the subsequent parsing 
    into tabular "gene family assignment" files.

SYNOPSIS
  gfa.sh -l filepath_list

  Required:
           -l  - path to list of file paths to compressed (.gz) fasta files
           -c  - path to the config file

  Options: -h  - help
           -t  - threads (or may be set in config file; see comments below)

  VARIABLES set in config file:
    search_method - align_sample (against family protein files) or  hmmsearch (against family HMMs)
    hmmdb_path    - Name of hmmer HMM data files (sans .h suffixes). HMM database should be in the data directory
    prot_fam_path - Path to multi-fasta protein file representing all gene families
    evalue        - E-value threshold for hmmsearch [1e-10]
    consen_iden   - Minimum identity threshold for mmseqs matches [0.30]
    clust_cov     - Minimum alignment coverage for mmseqs matches [0.40]
    data_dir      - Directory for HMM database [data]; also a reasonable place to put the list of fasta filepaths
    work_dir      - Directory where all work will be done
    threads       - Value for hmmsearch --cpu. [2] Recommended to be small (2), as most of the parallelism comes
                    from running on multiple query files, determined by available threads.

AUTHOR
    Steven Cannon <steven.cannon@usda.gov>
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
dependencies='hmmsearch mmseqs'
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

LIST=$(realpath "$filepath_list")
THREADS="$threads"
WORK_DIR="$work_dir"

NPROC=$( ( command -v nproc > /dev/null && nproc ) || getconf _NPROCESSORS_ONLN)

export NPROC=${NPROC:-1}
if [[ $NPROC -gt 1 ]]; then NPROC_PER_THREAD=$(( NPROC / THREADS )); else NPROC_PER_THREAD=1; fi
if [[ $NPROC_PER_THREAD -lt 1 ]]; then NPROC_PER_THREAD=1; fi
echo "NPROC: $NPROC"
echo "THREADS: $THREADS"
echo "NPROC_PER_THREAD: $NPROC_PER_THREAD"
echo "WORK_DIR: $WORK_DIR"

WD=$(realpath "$WORK_DIR")
mkdir -p "${WD}"

run_ingest() {
  echo
  echo "== Copy files into working directory, and uncompress"
  cd "${WD}" || exit
  mkdir -p 00_fasta 01_hmmsearch 02_gfa

  query_list=()
  ct=0
  
  echo "$LIST"
  cat "$LIST"
  echo
  cat "$LIST" | while read -r filepath; do
    if [[ "$filepath" == \#* ]]; then
      break
    else 
      if [ -f "$filepath" ]; then
        if [[ "$filepath" =~ \.gz$ ]]; then
          file=$(basename "$filepath" .gz)
          query_list+=("$file")
          ct=$((ct+1))
          printf '%s\n' "$ct ${query_list[@]}"
          if [ -f 00_fasta/"$file" ]; then
            echo "  Uncompressed file exists already at 00_fasta/$file"
          else
            echo "  Uncompressing input to 00_fasta/$file"
            gzip -dc "$filepath" > 00_fasta/"$file"
          fi
        else
          file=$(basename "$filepath")
          query_list+=("$file")
          ((ct++))
          printf '%s\n' "$ct ${query_list[@]}"
          echo "  File appears to be uncompressed; use as-is: $file"
          cat "$filepath" > 00_fasta/"$file"
        fi
      else
        echo "  WARN: File $filepath was not found."
      fi
    fi
  done

  echo "Printing query_list ====="
  printf '%s\n' "$ct ${query_list[@]}"
  echo "====="
}

run_search_families() {
  cd "${WD}" || exit
  echo "Search provided annotation sets (protein) against family consensus sequences or HMMs"

  if [ -d 33_fam_match ]; then rm -rf 33_fam_match; fi
  mkdir -p 33_fam_match

  if [ "$search_method" == "align_sample" ]; then
    mkdir -p 33_mmseqs_tmp

    if [ ! -f "${prot_fam_path}" ]; then
      echo "To place the supplementary sequences into gene families, the directory ${prot_fam_path} must exist"
      echo "(containing protein sequence in each family, from a previous run of \"pandagma fam\")."
      echo "The specified prot_fam_path was: ${prot_fam_path}"
      exit 1
    fi

    SEQTYPE=1 # 3=nuc; 1=pep

    target_ext="${prot_fam_path##*.}"
    target_base=$(basename "$prot_fam_path" ".$target_ext")

    for filename in ${query_list[@]}; do
      ext="${filename##*.}"
      base=$(basename "$filename" ."$ext")

      echo "TEST: $base $ext"

      echo "  Search $base.$ext against a sampling of sequences in families in target to determine best family match"
      echo "  mmseqs easy-search 00_fasta/$filename $prot_fam_path 33_fam_match/$base.x.$target_base.m8 \\"
      echo "    33_mmseqs_tmp --search-type ${SEQTYPE} --cov-mode 0 -c ${clust_cov} -v 1"
      mmseqs easy-search "00_fasta/$filename" \
                         "$prot_fam_path" \
                         33_fam_match/"$base".x.$target_base.m8 \
                         33_mmseqs_tmp \
                         --search-type ${SEQTYPE} --cov-mode 0 -c "${clust_cov}" -v 1
    done

    # In the regex below, the targets look like GeneFamID__GeneID  ... and the __GeneID is discarded (leaving 12 columns)
    for filepath in 33_fam_match/*.m8; do
      base=$(basename "$filepath")
      < "$filepath" perl -pe 's/^(\S+)\t([^_]+)__\S+/$1\t$2/' | # sort -k1,1 -k2,2 | calc_avg_bitscore_per_fam_hit.awk |
        sort -k1,1 -k12nr,12nr | top_line.awk > 33_fam_match/"$base".top &
      if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi

      # Once done with finding the top match per query, remove the large m8 files.
      #rm $filepath
    done
    wait

    echo "== Place sequences into families based on top protein match, using the consen_iden threshold of $consen_iden."
    cat 33_fam_match/*.m8.top |
      awk -v IDEN="${consen_iden}" '$3>=IDEN {print $2 "\t" $1}' |
      sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk > 34_sup_vs_fam_consen.clust.tsv

    echo "== Derive gene family assignment (gfa) files"
    for file in 33_fam_match/*.m8.top; do
      base=$(basename "$file" .x.target_base.m8.top)
      echo "$base";
      cat "$file" |
        awk -v OFS="\t" 'BEGIN{print "#gene", "family", "protein", "e-value", "score"}
                           {print $1, $2, $1, $11, $10}' |
          strip_spliceform.pl > "33_fam_match/$base.gfa.tsv"
    done

  elif [ "$search_method" == "hmmsearch" ]; then

    for filename in ${query_list[@]}; do
      ext="${filename##*.}"
      base=$(basename "$filename" ."$ext")
      # Set --cpu 2 for hmmsearch because it is IO bound. Rather, get parallelization by running the searches in parallel.
      echo "== Search $base.$ext against a HMMer database (compressed and indexed HMMs)"
      echo "    hmmsearch -E $evalue --cpu 2 --tblout 33_fam_match/$base.hmmsearch.tbl -o /dev/null ${hmmdb_path} 00_fasta/${filename}"
      hmmsearch -E "$evalue" --cpu 2 --tblout 33_fam_match/"$base".hmmsearch.tbl -o /dev/null "${hmmdb_path}" 00_fasta/"${filename}" &
      # allow to execute up to $NPROC/2 in parallel
      if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
    done
    wait
    echo "== Parse tabular hmmsearch data, picking the top match per query"
    for filepath in 33_fam_match/*.hmmsearch.tbl; do
      file=$(basename "$filepath" .hmmsearch.tbl)
      base=$(echo "$file" | perl -pe 's/(.+)\.\w+/$1/') # strips suffix like .protein_primary or .protein
      hmmdb_name=$(basename "$hmmdb_path")
      awk -v OFS="\t" '$1!~/^#/ {print $1, $3, $1, $5, $6}' "$filepath" |
        sort -k1,1 -k5nr,5nr | top_line.awk |
          strip_spliceform.pl > 33_fam_match/"${base}"."${hmmdb_name}".gfa.tsv &
       # allow to execute up to $NPROC in parallel
       if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi

      # Once done with finding the top match per query, remove the large hmmsearch.tbl files.
      rm "$filepath"
    done
    wait

    echo "== Place sequences into families based on top protein match to HMM, using the E-value threshold of $evalue."
    cat 33_fam_match/*.gfa.tsv |
      awk -v EVAL="${evalue}" '$4<=EVAL {print $2 "\t" $1}' |
      sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk > 34_sup_vs_fam_consen.clust.tsv

  else
    echo "Unrecognized search_method: $search_method."
    echo "Expected values are \"align_sample\ or \"hmmsearch\"."
    exit 1
  fi

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 34_sup_vs_fam_consen.clust.tsv > 34_sup_vs_fam_consen.hsh.tsv
}

########################################
# Main program

run_ingest
run_search_families

export dependencies='hmmsearch famsa'
NPROC=$( command -v nproc > /dev/null && nproc )
if [[ $NPROC -lt 1 ]]; then NPROC=1; fi
export NPROC=${NPROC:-1}

echo "NPROC: $NPROC"

declare consen_iden clust_cov evalue search_method prot_fam_path hmmdb_path

exit 0

