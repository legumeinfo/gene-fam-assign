# gene-fam-assign
Use hmmsearch to assign proteomes to LIS gene families

The code in this repository manages the assignment of gene sets (proteomes) to gene families, represented
as hmmer hmm "databases" (compressed and indexed HMMs), for the LIS, SoyBase, and PeanutBase projects
(or generally, for any set of proteomes against an HMM collection).

Also see https://github.com/legumeinfo/lis_gfa/ , which accomplishes the same task, using a singularity container
rather than a conda environment. The approaches mainly differ in the hmmsearch intermediate formats and parsing approaches.
The management of parallelization is also somewhat different -- the present repo relying more on parallelization
via spawning multiple `hmmsearch` tasks and `lis_gfa` relying on the `--cpu` parameter of `hmmsearch`.  The present approach 
should be faster, considering IO limitations of `hmmsearch` -- but we have not benchmarked the two methods yet.

## Installation

The only required dependency is the [hmmer package](http://hmmer.org). Install this and put
the executables on your path however you prefer. One way:

To create a conda environment called `hmmer` from the environment.yml in this repository:

    conda env create

Then activate the conda environment:
```
    conda activate pandagma
      #or
    source activate pandagma
```

## Preparing for a run by creating the HMM database (if not yet created)

First (you're on your own for this step), generate the "hmmpress"ed HMM target. Do this by concatenating
the hmmer-format HMMs into a file, then running hmmpress on that file. Put the file into the data directory of
this repository (or if in another location, indicate the path in the config directory).

## Running the program

Modify the few lines of the config -- at least specifying the hmm database name.

Then call the driver script, `gfa.sh`, with at least the `-l lis.protein_files` parameter (the list of filepaths to the 
compressed fasta files). The (gzipped) fasta files should be local to your filesystem.

It is presumed that the work will be done on a HPC resource, called via a slurm or comparable job submission script. 
See the example: `batch_gfa.sh`

The job script accomplishes three things (apart from requesting system resources):
 - Loads the dependencies -- comprised of the hmmer package (the example script uses conda)
 - Adds scripts (in bin/) to the PATH
 - Calls the driver script, e.g. `gfa.sh -l data/lis.protein_files`

Call the batch script like so (assuming slurm and the batch script by this name):
```
  sbatch batch_gfa.sh
```
