# gene-fam-assign
Use hmmsearch to assign proteomes to LIS gene families

The code in this repository manages the assignment of gene sets (proteomes) to gene families, represented
as hmmer hmm "databases" (compressed and indexed HMMs), for the LIS, SoyBase, and PeanutBase projects
(or generally, for any set of proteomes against an HMM collection).

Also see https://github.com/legumeinfo/lis_gfa/ , which accomplishes the same task, using a singularity container
rather than a conda environment. The approaches mainly differ in the hmmsearch intermediate formats and parsing approaches.
The management of parallelization is also somewhat different -- the present repo relying more on parallelization
via spawning multiple `hmmsearch` tasks and `lis_gfa` relying on the `--cpu` parameter of `hmmsearch`.  The present approach 
should be faster, considering IO limitations of `hmmsearch`, though we have not made head-to-head comparisons of the two.

## Installation

The only required dependency is the [hmmer package](http://hmmer.org). Install this and put
the executables on your path however you prefer. One way:

To create a conda environment called `hmmer` from the environment.yml in this repository:

    conda env create

Then activate the conda environment (the command for environment activation may differ depending on the platform):
```
    conda activate hmmer
      #or
    source activate hmmer
```

## Preparing for a run by creating the HMM database (if not yet created)

A "hmmpress"ed HMM target is required. In the context of the legumeinfo/soybase/peanutbase project, such
a target is available at `/project/legume_project/common_data/genefamilies`, named as `legume.fam3.VLMQ`.
This target has the HMMs for all gene families in the "legume family3" set.

In other contexts (you have your own gene families and associated HMMs), you're on your own for generating
the "hmmpress"ed HMM target. Do this by concatenating the hmmer-format HMMs into a file, then running 
hmmpress on that file. Put the file into the data directory of this repository (or if in another 
location, indicate the path in the config directory).

## Running the program

Check the provided config file and make sure that the `hmmdb` and `data_dir` point to your HMM database.

Then call the driver script, `gfa.sh`, with at least the `-l data/lis.protein_files` and the `-c config/gfa.conf`
parameters (the list of filepaths to the list of compressed fast files to use as queries, and the config file).
The (gzipped) fasta files should be local to your filesystem.

It is presumed that the work will be done on a HPC resource, called via a slurm or comparable job submission script. 
See the example: `batch_gfa.sh`

The job script accomplishes three things (apart from requesting system resources):
 - Loads the dependencies -- comprised of the hmmer package (the example script uses conda)
 - Adds scripts (in bin/) to the PATH
 - Calls the driver script, e.g. `gfa.sh -l data/lis.protein_files -c config/gfa.conf`

Call the batch script like so (assuming slurm and the batch script by this name):
```
  sbatch batch_gfa.sh
```

## Optimizing the search speed

The greatest search speed will be achieved with a relatively large number of fasta query files and
a comparable number of threads. The hmmsearch `--cpu` parameter can also be set, but as hmmsearch
is highly IO-bound, increasing this parameter beyond 1 or 2 doesn't give net improvements in speed.
In fact, in our tests, we found no increase in speed for --cpu 2 over --cpu 1. Nevertheless, 
trusting the hmmer documentation, we recommend using `--cpu 2` (as set in the example config file), 
unless the number of data sets to be analyzed is greater than 2xNPROC; in that case, set --cpu 1.

Example: if you have 20 annotation sets to use as queries against the HMMdb, then leave `--cpu 2`
and make sure that `#SBATCH --ntasks-per-node=40` (or higher), to ensure that 20 searches run in parallel.

