# seq-analysis-tools
Collection of simple tools to help with data acquisition, data reuse and sequence analysis

# Iterative_SRA_mapper
Python script to help with efficient and parallelized downloading, conversion, and mapping steps of publicly available RNA-Seq data. Allows users to generate count tables for species, based on a coding sequence (cds) file and ncbi SRA datasets (limited to paired end RNA-Seq data for now). Uses multiprocessing, sra-toolkit, and kallisto to efficiently determine and store expression values for organism of interest.  

options:

  -h, --help            show this help message and exit
  
  --in IN_FILE          Input file containing SRA IDs, one per line
  
  --out OUT_DIR         Output directory to store mapping results
  
  --n NUM_CORES         Number of threads being used, will be rounded down to
                        multiple of 4, default = half of available threads
                        
  --cds CDS             Path to cds file to generate kallisto index (.idx) for
                        read mapping. Can rerun with existent index
                        
  --sra_cache SRA_CACHE_DIR
                        Path to SRA cache directory

Requirements:
1. sra toolkit (https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) 
2. kallisto (https://github.com/pachterlab/kallisto)
3. Python packages: os, argparse, subprocess, multiprocessing

The external programs need to be in PATH (executable from command line). On Ubuntu, installation through apt is recommended. SRA toolkit configuration might be needed ("vdb-config -i") before usage.

How to use:
1. Find cds file for species of interest, used as target of kallisto mapping (--cds)
2. In NCBI Sequence Read Archive (https://www.ncbi.nlm.nih.gov/sra), find data for species of interest, filter to RNA, paired, fastq (only tested on Illumina data)
3. Download Accession list for use as input file (--in)
4. Make sure that cache of sra-toolkit is set and has sufficient space (--sra-cache, default="~/.ncbi/public/sra")
5. Pick output directory (--out) and number of threads to use (--n)
6. run script with python3

Example command: 
python3 Iterative_SRA_align.py \
--in /path/to/SRA_file \
--out /path/to/output_dir/ \
--cds /path/to/cds.fa \
--sra_cache /path/to/sra_cache \
--n 16

Output will be one subdir for each SRA accession, containing the count table ("abundance.tsv") and a kallisto log file.
