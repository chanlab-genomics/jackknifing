# Jackknifing Workflow

## Jackknifing

To build our jackknife tree we first need to produce some jackknife samples (surprise surprise!). The jackknifing step mostly makes use of `jackknife.py` found in the top level directory. In short `jackknife.py` reads in a fasta file a spits out a reduced version. `jackknife.py --help` does a pretty good job of explaining what command line arguments it's expecting, so I've taken the liberty of copying and pasting the output of the `--help` flag here
```
usage: jackknife.py [-h] --input_paths INPUT_PATHS [INPUT_PATHS ...] [--output_path OUTPUT_PATH] [-v] [--portion PORTION] [--chunk_size CHUNK_SIZE] [--threads THREADS]

Randomly removes a portion of data from a fasta file.

optional arguments:
  -h, --help            show this help message and exit
  --input_paths INPUT_PATHS [INPUT_PATHS ...]
                        A list of fasta files to perform the jacknifing process
  --output_path OUTPUT_PATH
                        A folder to write the reduced fasta files.
  -v, --verbose         Include to run the script in verbose mode. Useful for checking progress. WARN: This may cause the program to run slower.
  --portion PORTION     The portion of data to be removed. The default is a 40 percent reduction.
  --chunk_size CHUNK_SIZE
                        The size of the chunks that get randomly removed from sequences. Default is a chunk size of 100.
  --threads THREADS     The number of threads used to run the jackknife algorithm. If 0 threads are specified then it will default to os.cpu_count().
```
Here's an example use of `jackknife.py`
```
python3 jackknife.py --input_path ~/AEH.fasta --output_path ~/jn_yeast --portion=40 --chunk_size=100 --threads=4
```
If you're going through this workflow by yourself, chances are you've got hundereds of genomes to jackknife across many different samples. To see how this process could be automated check `batch\array_jn\jn_array.sh` and `batch\array_jn\submit_jn_wait.sh`. The `batch\array_jn\jn_array.sh` file is a PBS job which produces a jackknife sample for each fasta file found in `ARRAY_TARGET` and saves the outputs to `OUTPUT_DIR`. You might want to change around the number of cpus, memory and walltime for the job file depending on how big/how many fasta you have. The `batch\array_jn\submit_jn_wait.sh` will automatically submit multiple `batch\array_jn\jn_array.sh` jobs with different sample values.

## Jellyfish

The next step in this workflow is to use jellyfish as well as some custom purpose python scripts to perform some analysis on our jackknifed samples. This means you will either have to install the `jellyfish` program on your computer or in you're working on a hpc, you could try running `module load jellyfish/*version*` or just `module load jellyfish` to load jellyfish (to check if they have jellyfish at all run `module avail` to see if they have have any version of jellyfish for you to use). Once jellyfish has been setup you will need to run the following string of commands on each jacknified fasta file
```

k=*kmer-count* # eg k=21
s=10000000000

file=*jackknifed fasta file* # eg file=~/AEH_red_40.fasta
jellyfish count -m $k -s $s -t $CPUS_PER_TASK -o $file.$k.jf $file && \
jellyfish dump -ct $file.$k.jf | sort -k1,1 | python2 jf_scripts/Kmers_2_NumbericRepresentation.py -o $file.${k}mer.nkc.gz && \
python2 jf_scripts/Composition_of_InputSeqs.py --fasta $file --freq $file.CharFreq && \
touch $file.done &
```
Note that these commands require a lot of memory and therefore are best run as batch jobs. Again doing this for each individual jackknifed file can be very tedious. To see how this process can be automated, see `batch\jellyfish_jobs\submit_jellyfish.sh` and `batch\jellyfish_jobs\run_jellyfish.sh`. The `batch\jellyfish_jobs\run_jellyfish.sh` is a PBS job file which runs the above string of commands on every fasta file found within `ARRAY_TARGET`. The `batch\jellyfish_jobs\submit_jellyfish.sh` bash script simply submits the `batch\jellyfish_jobs\run_jellyfish.sh` batch script for different samples. A successful completion of these commands for a fasta file it should have produced a `.*kmer-count*.jf` `.*kmer-count*mer.nkc.gz`, `.CharFreq` and a `.done` file. For example here's what a directory containing just `~/AEH_red_40.fasta` with `*kmer-count*=21` would look like
```
AEH_red_40.fasta               
AEH_red_40.fasta.21.jf         
AEH_red_40.fasta.21mer.nkc.gz  
AEH_red_40.fasta.CharFreq      
AEH_red_40.fasta.done
```