#!/usr/bin/env bash

set -o errexit -o errtrace -o nounset -o pipefail -o posix

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

# to help assign a heredoc value to a variable. The internal line return is intentional.
define(){ o=; while IFS=$'\n' read -r a; do o="$o$a"'
'; done; eval "$1=\$o"; }

define HELP_DOC <<'EOS'
NAME
  make_reduced_fam_sample.sh

  Given a collection of gene family files, each consisting of multifasta protein sequences 
  for the family, make a single multifasta file that can be used as a target for searching 
  annotation sets against, for the purpose of assigning gene families.

  Three methods are available: "mmseqs", "random", and "hybrid".

  - The "mmseqs" method calculates a nonredundant set for each family, such that no two 
    sequences share more than a specified sequence identity. This is sensible, but slow.

  - The random method selects "num_to_keep" sequences from the gene family for use in the 
    combined target file.  This is fast, and should be fine in most cases, providing that "num_to_keep" 
    is sufficiently large that all sequence types are likely to be represented in the sample.

  - The hybrid method uses "mmseqs" for families with more than "num_to_keep" sequences, 
    and "random" for families in which the number of sequences is fewer than or equal to num_to_keep.

  The sampled sequences are concatenated into into a file of protein sequences with a defline format like so,
  where the delimiter "__" separates the family ID and the representative gene ID:
    >Legume.fam3.00001__singl.CAF01.gnm1.ann1.evm_model_Chr03_2437
    >Legume.fam3.00001__glyma.Wm82.gnm4.ann1.Glyma.06G281100.1
    >Legume.fam3.00001__tripr.MilvusB.gnm2.ann1.mRNA40501

SYNOPSIS
  make_reduced_fam_sample.sh -f path/to/gene/families [options]

  Required:
    -f  Path to gene families (each a multifasta protein file)
         "/project/legume_project/common_data/genefamilies/legume.fam3.VLMQ.sup1B_proteomes"

  Optional:
    -w  work directory [work_dir]
    -o  output filename [combined_sample.faa]
    -m  method; "mmseqs", "random", or "hybrid". [hybrid]
    -n  num_to_keep; no more than this number will be kept in each sampled family file. Use with method "random" and "hybrid". [100]
    -i  minimum identity for selecting representative sequences using mmseqs; is used for both --min-seq-id and cluster coverage [0.8]

AUTHOR
    Steven Cannon <steven.cannon@usda.gov>
EOS

if [ "$#" -eq 0 ]; then
  echo >&2 "$HELP_DOC" && exit 0;
fi

# Default flags
prot_fam_path=""
work_dir="work_dir"
method="hybrid" 
num_to_keep=100
output_filename="combined_sample.faa"
min_id="0.8"

while getopts "f:w:m:n:o:i:h" opt
do
  case $opt in
    f) prot_fam_path=$OPTARG; echo "prot_fam_path: $prot_fam_path" ;;
    w) work_dir=$OPTARG; echo "work_dir: $work_dir" ;;
    m) method=$OPTARG; echo "method: $method" ;;
    n) num_to_keep=$OPTARG; echo "num_to_keep: $num_to_keep" ;;
    o) output_filename=$OPTARG; echo "output_filename: $output_filename" ;;
    i) min_id=$OPTARG; echo "min_id: $min_id" ;;
    h) echo >&2 "$HELP_DOC" && exit 0 ;;
    *) echo >&2 echo "$HELP_DOC" && exit 1 ;;
  esac
done

# Check for existence of third-party executables
missing_req=0
dependencies='mmseqs'
for program in $dependencies; do
  if ! type "$program" &> /dev/null; then
    echo "Warning: executable $program is not on your PATH."
    missing_req=$((missing_req+1))
  fi
done

WD=$(realpath "$work_dir")
mkdir -p "${WD}"

cd "${WD}" || exit

if [ -d 31_mmseqs_sample ]; then rm -rf 31_mmseqs_sample; fi
mkdir -p 31_mmseqs_sample

mkdir -p 31_mmseqs_sample_tmp 31_mmseqs_sample

if [ ! -d "${prot_fam_path}" ]; then
  echo "The directory ${prot_fam_path} must exist, and must contain a set of multifasta files,"
  echo "each containing protein sequence for a gene family."
  echo "The specified prot_fam_path was: ${prot_fam_path}"
  exit 1
fi

# Find the total number of files
total_files=$(find "$prot_fam_path" -type f | wc -l)
count=0
rep_suffix="_rep_seq.fasta"

echo "== Find representative sequences for each gene family, using method $method. Output will go to $output_filename"
# The output filename is comparable to the file 19_palmp_sample.faa from pandagma_fam.sh ...
# which is an abbreviation of 19_pan_aug_leftover_merged_prot
cat /dev/null > "$output_filename"
if [[ "$method" == "random" ]]; then
  num_to_keep=100 # Don't search against more than this number of sequences in a family
  echo "  Prepare to search against a random sampling of $num_to_keep sequences from each family"
  echo "  (or as many as available if there are fewer than $num_to_keep in the family)."

  for filepath in "${prot_fam_path}"/*; do
    base=$(basename "$filepath")
    fasta_to_table.awk "$filepath" | shuf |
      awk -v FAM="$base" -v MAX=$num_to_keep 'ct<MAX {print ">" FAM "__" $1; print $2; ct++}' >> "$output_filename"

    ((++count))
    if (( count % 100 == 0 )) || (( count == total_files )); then
      percent=$(( count * 100 / total_files ))
      echo -ne "    Processing as random $num_to_keep: $count/$total_files ($percent%)\r"
    fi
  done
  echo ""
elif [[ $method == "mmseqs" ]]; then
  for filepath in "${prot_fam_path}"/*; do
    base=$(basename "$filepath")
    
    mmseqs easy-cluster "$filepath"  "31_mmseqs_sample/$base"  31_mmseqs_sample_tmp --min-seq-id "$min_id" -c "$min_id" --cov-mode 0 -v 1

    rep_seq="$base$rep_suffix"
    #echo "$base $rep_seq"
    fasta_to_table.awk "$rep_seq" | awk -v FAM="$base" '{print ">" FAM "__" $1; print $2}' >> "$output_filename"
  
    ((++count))
    if (( count % 100 == 0 )) || (( count == total_files )); then
      percent=$(( count * 100 / total_files ))
      echo -ne "    Processing with mmseqs at $min_id: $count/$total_files ($percent%)\r"
    fi
  done
  echo ""
elif [[ $method == "hybrid" ]]; then
  for filepath in "${prot_fam_path}"/*; do
    base=$(basename "$filepath")
    num_seqs=$(grep -c "^>" "$filepath")
    
    if [[ $num_seqs -gt $num_to_keep ]]; then
      mmseqs easy-cluster "$filepath"  "31_mmseqs_sample/$base"  31_mmseqs_sample_tmp --min-seq-id "$min_id" -c "$min_id" --cov-mode 0 -v 1

      rep_seq="$base$rep_suffix"
      echo "$base $rep_seq"
      fasta_to_table.awk "31_mmseqs_sample/$rep_seq" | awk -v FAM="$base" '{print ">" FAM "__" $1; print $2}' >> "$output_filename"
    else
      fasta_to_table.awk "$filepath" | shuf | awk -v FAM="$base" '{print ">" FAM "__" $1; print $2; ct++}' >> "$output_filename"
    fi

    ((++count))
    if (( count % 100 == 0 )) || (( count == total_files )); then
      percent=$(( count * 100 / total_files ))
      echo -ne "    Processing with mmseqs at $min_id: $count/$total_files ($percent%)\r"
    fi
  done
  echo ""
else 
  echo "The variable \"method\" must be either \"mmseqs\", \"random\", or \"hybrid\". The provided value was $method"
  exit 1
fi

exit 0

