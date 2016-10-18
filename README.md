# Satsuma2

Satsuma2 is an optimsed version of Satsuma, a tool to reliably align large and complex DNA sequences providing maximum sensitivity (to find all there is to find), specificity (to only find real homology) and speed (to accomodate the billions of base pairs in vertebrate genomes). Satsuma2 adresses these issues through three novel strategies:

* cross-correlation, implemented via fast Fourier transformation. 
* a match scoring scheme that eliminates almost all false hits.
* an asynchronous "battleship"-like search that enables fast whole-genome alignment.

Satsuma2 also interfaces with MizBee, a multi-scale synteny browser for exploring conservation relationships in comparative genomics data (<http://www.cs.utah.edu/~miriah/mizbee/Overview.html>).

Satsuma2 is implemented in C++ on Linux.

## Licensing
Satsuma2 is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.

## Citing Satsuma2

We plan to submit an application note that should be published during the summer of 2016. In the meantime, if you are using Satsuma2 for research that will be published before that, please contact us to discuss how you can cite the tool.

## Installation

Download the source code from <https://github.com/bjclavijo/satsuma2.git> and compile it using CMake v3.3+.  To run, Satsuma2 requires GCC v5.2+.  The binaries are generated in the bin/ directory.

NOTE: if you encounter the error "... undefined reference to `pthread_create'" during compilation, add the flag -pthread to CMakeLists.txt, i.e. change:

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -std=c++14 -O3 -w")
to:
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -pthread -std=c++14 -O3 -w")

## Quick start

1. Set the SATSUMA2_PATH environment variable to point to the directory containing the binaries.
2. Configure satsuma_run.sh to match your job submission system (LSF, PBS, SLURM etc.)
3. Run SatsumaSynteny2 with default parameters (uses one single-threaded slave);

```
./SatsumaSynteny2 -q query.fa -t target.fa -o output_dir
```

## Running Satsuma2

As Satsuma2 calls other executables (HomologyByXCorr, MergeXCorrMatches etc.), you need to set an environment variable to tell the software where to find the binaries;

```
export SATSUMA2_PATH=/path/to/binaries
```

Running SatsumaSynteny2 with no command-line arguments shows the available parameters;

```
$ ./SatsumaSynteny2

./SatsumaSynteny2: Wrapper around high-sensitivity alignments.

Available arguments:

-q<string> : query fasta sequence
-t<string> : target fasta sequence
-o<string> : output directory
-l<int> : minimum alignment length (def=0)
-t_chunk<int> : target chunk size (def=4096)
-q_chunk<int> : query chunk size (def=4096)
-sl_mem<int> : memory requirement for slaves (Gb) (def=100)
-do_refine<bool> : refinement steps (def=0)
-min_prob<double> : minimum probability to keep match (def=0.99999)
-cutoff<double> : signal cutoff (def=1.8)
-prob_table<bool> : approximate match prob using a table lookup in slaves (def=0)
-min_matches<int> : minimum matches per target to keep iterating (def=20)
-m<int> : number of jobs per block (def=4)
-slaves<int> : number of processing slaves (def=1)
-threads<int> : number of working threads per processing slave (def=1)
-km_mem<int> : memory required for kmatch (Gb) (def=100)
-km_sync<bool> : run kmatch jobs synchronously (def=1)
-seed<string> : loads seeds and runs from there (kmatch files prefix) (def=)
-min_seed_length<int> : minimum length for kmatch seeds (after collapsing) (def=24)
-max_seed_kmer_freq<int> : maximum frequency for kmatch seed kmers (def=1)
-old_seed<string> : loads seeds and runs from there (xcorr*data) (def=)
-pixel<int> : number of blocks per pixel (def=24)
-nofilter<bool> : do not pre-filter seeds (slower runtime) (def=0)
-dups<bool> : allow for duplications in the query sequence (def=0)
-dump_cycle_matches<bool> : dump matches on each cycle (for debug/testing) (def=0)
```

Required parameters are query FASTA (-q), target FASTA (-t) and output directory (-o).  The query and target sequences are chunked (based on the -t\_chunk and -q\_chunk parameters) and KMatch is used to detect aligning regions between chunks.  The number of chunks generated depends on the length of your query and target sequences.  The amount of memory reserved for KMatch can be modified using the -km\_mem parameter which defaults to 100Gb.

SatsumaSynteny2 despatches slave processes to compare the chunks which run asynchronously.  The number of slaves, threads per slave and memory limit per slave are specified using the -slaves, -threads and -sl\_mem parameters.  The default is to run one single-threaded slave using 100Gb of memory.  Slaves can be run on a single machine or submitted via a job submission system such as LSF, PBS or SLURM and the satsuma\_run.sh file is used by SatsumaSynteny2 to start the slaves.  Before running SatsumaSynteny2, you need to modify this file to suit your HPC environment by commenting out the lines you don't need with #. You will also need to change 'QueueName' to a queue that exists on your system.  For example, to run on SLURM your file should look like this;

```
#!/bin/sh

# Script for starting Satsuma jobs on different job submission environments
# One section only should be active, ie. not commented out

# Usage: satsuma_run.sh <current_path> <kmatch_cmd> <ncpus> <mem> <job_id> <run_synchronously>
# mem should be in Gb, ie. 100Gb = 100

# no submission system, processes are run locally either synchronously or asynchronously
#if [ "$6" -eq 1 ]; then
#  eval "$2"
#else
#  eval "$2" &
#fi

##############################################################################################################
## For the sections below you will need to change the queue name (QueueName) to one existing on your system ##
##############################################################################################################

# qsub (PBS systems)
#echo "cd $1; $2" | qsub -V -qQueueName -l ncpus=$3,mem=$4G -N $5

# bsub (LSF systems)
#mem=`expr $4 + 1000`
#bsub -o ${5}.log -J $5 -n $3 -q QueueName -R "rusage[mem=$mem]" "$2"

# SLURM systems
echo "#!/bin/sh" > slurm_tmp.sh
echo srun $2 >> slurm_tmp.sh
sbatch -p QueueName -c $3 -J $5 -o ${5}.log --mem ${4}G slurm_tmp.sh

```

### Notes  

* If SatsumaSynteny2 is run without a submission system, KMatch jobs will be launched synchronously in order to keep memory requirements low.  If you have plenty of memory available you can opt to run the KMatch jobs asynchronously (-km_sync 0).  KMatch requires a lot of memory and multiple KMatch processes running at the same time may cause SatsumaSynteny2 to abort if not enough memory is available.
* The parameters -km\_mem and -sl\_mem are only applied when using a job submission system.  We strongly recommend using a job submission system to run SatsumaSynteny2 which allows more control of the resource requirements of this software.
* If the output directory is not empty, SatsumaSynteny2 will not overwrite any files but exit with an error message.  
* Idling processes self-terminate after two minutes. The overall alignments will still complete, but using fewer processes.  
* If alignment runs locally but not on the server farm, check whether processes on the farm can communicate via TCP/IP.  
* Currently, the entire sequences are loaded into RAM by each process. For comparison of large genomes, we strongly recommend to make sure that the CPUs have enough RAM available (~ the size of both genomes in bytes). 

### Parameter choice, execution and data preparation  

* The default parameters should work well for most genomes.
* SatsumaSynteny2 runs most efficiently on either multi-processor machines or on clusters that are tightly coupled (fast access to files shared by the control process and the slaves)
* Especially for larger genomes, we recommend leaving one CPU dedicated to the control process SatsumaSynteny2.
* For larger genomes (>1Gb), we recommend using one chromosome of one genome as the query sequence and the entire other genome as the target sequence, and process alignments one query chromosome at a time. 
* To include large-scale duplications in the query sequence (in addition to the target sequence), use the option –dups.
* If using the option –nofilter, the number of initial searches (-ni) should be higher than the number of processes (-n) to ensure that subsequent processes have sufficient seeds. Note that initial searches will be queued to a number of processes specified by -n.
* When many processes search a tight space, the number of pixels per CPU (-m) should be small (e.g. ‘–m 1’ as in the sample script/data set) to avoid unbalanced load (i.e. some processes get all the pixels while others are starved, since they overlap). However, a small value for –m increases inter-process communication, which should be a consideration when deploying hundreds of processes.


## Making SatsumaSynteny2 converge (a temporary note)

Given a new and more exhaustive convergence model, which is still under active development, Satsuma2 may fail to converge into a single final result, and rather enter an iteration cycle, where lots of small (or not so small) changes are made to the general alignment search strategy. Instead of hiding this behaviour under a fixed cutoff after a number of iterations, we have chosen to expose it, and allow the user to examine and choose the intermediate result that best suits the biological question.

For this reason, we have introduced the parameter -dump\_cycle\_matches, which will produce an output file on each cycle. Because these output files are not particularly large and contain the whole information of the cycle, we recommend to turn this parameter on, unless you're running on datasets where you already know that the convergence setup will work correctly. You can then examine the convergence (probably using the MatchesByFeature tool if one of your genomes is annotated) and decide which solution(s) best suits your objectives.

We are working on generalising the convergence model so it behaves well under most circumstances, but still this will always be a recommendation when starting to run SatsumaSynteny2 in new scenarios.


## Output files

**<outdir>/satsuma_summary.chained.out: final coordinates**

```
Contents:  
	target sequence name  
	first target base  
	last target base  
	query sequence name  
	first query base  
	last query base  
	identity  
	orientation  
	
EXAMPLE:
chrX	5947	6164	chrX	9153	9360	0.626728	+ 
chrX	6270	6452	chrX	9472	9654	0.576923	+
```

Note: ‘space’ in fasta names is permissible for alignment, but all spaces will be replaced with “_” in the output files.

**<outdir>/MergeXCorrMatches.chained.out: final readable alignments**

```
EXAMPLE:
Query chr24 [29727636-29727834] vs target scaffold_24 [1206-1404] + length 198 check 198
Identity (w/ indel count): 52.5253 %
-------------------------------------------------------------------------------

TCCCCACTTCTAAAGTAAACTGCACATAGGGACTTCTTTCCAAAGAGCACAGTCTGGAAAGGAGGGAAAAACAATTTTAC
       ||  |||      |  | ||     ||   ||||||| |  || ||  | |||||    || || ||||| ||
ATATATTTTTAAAATATCTATTAAAATCAAACCTATGTTCCAAATATTACGGTACGAAAAGGGAAAAATAAGAATTTCAC

WYMYMWYTTYWAAAKWWMWMTKMAMATMRRRMCTWYKTTCCAAAKAKYACRGTMYGRAAAGGRRRRAAWAASAATTTYAC AMB
-------------------------------------------------------------------------------

AGTCTATAAACCTGATAAACACTACCTCAGCCAGGTGCTCAAGGGCAACATCAAGACTCGTAAGTCATGTTGATAGTAGA
||   | ||  |||  |  ||| |||| | ||| ||| ||||||  ||  |||    || | |||||  ||||||  |
AGCAAAGAAGTCTGGCAGTCACCACCTTAACCAAGTGATCAAGGTTAATGTCACTGATCATGAGTCACATTGATATAATG

AGYMWAKAARYCTGRYARWCACYACCTYARCCARGTGMTCAAGGKYAAYRTCAMKRMTCRTRAGTCAYRTTGATAKWAKR AMB
-------------------------------------------------------------------------------

TCCTATTGATATGCTTTGCAAGGACAGAGTAATTGACA
| |     ||||| | ||    ||  |     |   |
TACCCCCCATATGATGTGATGAGAAGGGCATTTCACCT

TMCYMYYSATATGMTKTGMWRRGAMRGRSWWWTYRMCW AMB
```
**<outdir>/xcorr\_aligns\_final.out - results in binary format for downstream analysis**	

## Other Satsuma tools
Run each tool with no arguments to see available options.

### Alignment tool

* ColaAlignSatsuma: realigns global alignments in satsuma format (summary coordinates file) using Cola, an efficient implementation of a collection of sequence alignment algorithms.

### Visualisation tools

* BlockDisplaySatsuma: takes a satsuma summary file and writes displayable blocks in MizBee format, see <http://www.cs.utah.edu/~miriah/mizbee/Overview.html> for how to display this using the MizBee Synteny Browser.
* ChromosomePaint: generates a comparative chromosome view in postscript format from the MizBee file generated by BlockDisplaySatsuma.
* MicroSyntenyPlot: generates a postscript visualisation of synteny from the HomologyByXCorr binary output file (xcorr\_aligns\_final.out).

### Scaffold synteny tools

* Chromosemble: runs a pipeline that scaffolds an assembly using synteny.
* OrderOrientBySynteny: orders and orients scaffolds according to a synteny map.

### Other useful tools
* MatchesByFeature: report matches by specific features using a GFF3 file.  To show matches to exon and CDS features defined in GFF file genome.gff using match files match1 and match2 use;

```
./MatchesByFeature genome.gff exon CDS - - - - - match1 match2
```
* ReverseSatsumaOut: swaps query and target columns in the satsuma output file.
* SatsumaToFasta: generates a FASTA file using a satsuma summary file from either the query or the target genome.
* SatsumaToGFF: generates a GFF3 file from a satsuma summary file. 




