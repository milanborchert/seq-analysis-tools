import os
import argparse
import subprocess
from multiprocessing import Pool

def delete_directory(path):
    """Recursively delete a directory and all its contents."""
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(path)

def remove_stale_lock_files(output_directory):
    """Remove stale lock files from the output directory."""
    for file_name in os.listdir(output_directory):
        if file_name.endswith(".lock"):
            file_path = os.path.join(output_directory, file_name)
            os.remove(file_path)
            print(f"Removed stale lock file: {file_path}")

def clear_sra_cache(sra_ids, sra_cache_dir):
    """Clear SRA cache directory of files related to given SRA IDs."""
    sra_cache_dir = os.path.expanduser(sra_cache_dir)
    
    for sra_id in sra_ids:
        files_to_remove = [
            os.path.join(sra_cache_dir, f"{sra_id}.sra"),
            os.path.join(sra_cache_dir, f"{sra_id}.sra.tmp"),
            os.path.join(sra_cache_dir, f"{sra_id}.sra.prf")
        ]
        
        print("Removing the following SRA cache files:")
        for file_path in files_to_remove:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"Removed file: {file_path}")

def download_and_convert_SRA(SRA_ID, output_directory, genome_index, sra_cache_dir):
    """Download, convert, align, and delete unnecessary files for a single SRA entry"""

    # Check and retrieve the subdirectory
    sra_dir, skip = check_and_prepare_subdir(SRA_ID, output_directory)
    if skip:
        return

    # Check if lock file exists for this ID
    lock_file = os.path.join(output_directory, f"{SRA_ID}.lock")
    if os.path.exists(lock_file):
        print(f"Skipping SRA ID {SRA_ID}. Already being processed by another instance.")
        return

    # Create lock file
    open(lock_file, 'a').close()

    # Check if FASTQ files already exist for this ID
    fastq_file1 = os.path.join(sra_dir, f"{SRA_ID}_1.fastq")
    fastq_file2 = os.path.join(sra_dir, f"{SRA_ID}_2.fastq")
    if os.path.exists(fastq_file1) and os.path.exists(fastq_file2):
        print(f"Skipping SRA ID {SRA_ID}. FASTQ files already exist.")
    else:
        # Use prefetch to download SRA files
        prefetch_cmd = ["prefetch", SRA_ID]
        print(f"Running prefetch for SRA ID: {SRA_ID}")
        try:
            subprocess.run(prefetch_cmd, cwd=sra_cache_dir, check=True)
            print(f"Prefetch completed for SRA ID: {SRA_ID}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading SRA ID {SRA_ID}: {e}")
            os.remove(lock_file)  # Remove lock file
            return

        # Use fasterq-dump to convert SRA files to FASTQ format
        fasterq_dump_cmd = ["fasterq-dump", "--outdir", sra_dir, "--threads", "4", "--split-files", SRA_ID]
        print(f"Running fasterq-dump for SRA ID: {SRA_ID}")
        try:
            subprocess.run(fasterq_dump_cmd, cwd=sra_dir, check=True)
            print(f"Fasterq-dump completed for SRA ID: {SRA_ID}")
        except subprocess.CalledProcessError as e:
            print(f"Error converting SRA ID {SRA_ID} to FASTQ: {e}")
            os.remove(lock_file)  # Remove lock file
            # Delete SRA files if fasterq-dump fails
            sra_file = os.path.join(sra_cache_dir, f"{SRA_ID}.sra")
            if os.path.exists(sra_file):
                os.remove(sra_file)
            return

    print("Split FASTQ files created:", fastq_file1, fastq_file2)

    # Check if kallisto output files already exist
    abundance_file = os.path.join(sra_dir, "abundance.tsv")
    if os.path.exists(abundance_file):
        print(f"Skipping kallisto for SRA ID {SRA_ID}. Abundance file already exists.")
    else:
        # Use kallisto to perform the alignment
        kallisto_cmd = ["kallisto", "quant", "--index", genome_index, "--output-dir", sra_dir, "--threads", "4", fastq_file1, fastq_file2]
        print(f"Running kallisto for SRA ID: {SRA_ID}")
        try:
            subprocess.run(kallisto_cmd, cwd=sra_dir, check=True)
            print(f"Kallisto completed for SRA ID: {SRA_ID}")
        except subprocess.CalledProcessError as e:
            print(f"Error running kallisto for SRA ID {SRA_ID}: {e}")
            # Remove lock file after processing
            os.remove(lock_file)
            # Delete SRA and FASTQ files if kallisto fails
            if os.path.exists(abundance_file):
                os.remove(abundance_file)
            if os.path.exists(fastq_file1):
                os.remove(fastq_file1)
            if os.path.exists(fastq_file2):
                os.remove(fastq_file2)
            sra_file = os.path.join(sra_cache_dir, f"{SRA_ID}.sra")
            if os.path.exists(sra_file):
                os.remove(sra_file)
            return

    # Delete unnecessary files after alignment
    try:
        os.remove(fastq_file1)  # Remove FASTQ file 1
        print("Deleted FASTQ file 1:", fastq_file1)
        os.remove(fastq_file2)  # Remove FASTQ file 2
        print("Deleted FASTQ file 2:", fastq_file2)
        sra_cache_file = os.path.join(sra_cache_dir, f"{SRA_ID}.sra")
        if os.path.exists(sra_cache_file):
            os.remove(sra_cache_file)  # Remove SRA cache file
            print("Deleted SRA cache file:", sra_cache_file)
    except Exception as e:
        print(f"Error deleting files for SRA ID {SRA_ID}: {e}")

    # Remove lock file after processing
    os.remove(lock_file)
    print("Removed lock file:", lock_file)

def check_and_prepare_subdir(SRA_ID, output_directory):
    """Check and prepare the subdirectory for the given SRA ID."""
    sra_dir = os.path.join(output_directory, SRA_ID)
    abundance_file = os.path.join(sra_dir, "abundance.tsv")
    
    if os.path.isdir(sra_dir):
        if os.path.isfile(abundance_file):
            print(f"Subdirectory for SRA ID {SRA_ID} exists and contains the abundance file. Skipping this ID.")
            return sra_dir, True
        else:
            print(f"Subdirectory for SRA ID {SRA_ID} exists but does not contain the abundance file. Deleting the subdirectory and its contents.")
            try:
                delete_directory(sra_dir)
                print(f"Deleted subdirectory {sra_dir} and its contents.")
            except Exception as e:
                print(f"Failed to delete subdirectory {sra_dir}: {e}")
                raise
    else:
        print(f"No subdirectory found for SRA ID {SRA_ID}.")
    
    os.makedirs(sra_dir)
    print(f"Created directory: {sra_dir}")
    return sra_dir, False

def main(args):
    input_file = args.in_file
    output_directory = args.out_dir
    cds = args.cds
    sra_cache_dir = os.path.expanduser(args.sra_cache_dir)

    if args.num_cores is None:
        args.num_cores = os.cpu_count() // 2
        print(f"Defaulting to half of available cores: {args.num_cores} cores.")
    
    num_cores = args.num_cores
    num_process = num_cores // 4
    print(f"{num_cores} threads were allocated, resulting in {num_process} processes with 4 threads each being spawned.")
    print(f"Actual thread usage: {int(num_process) * 4}")

    # Create output directory if it doesn't exist (moved earlier in the code)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print("Created output directory:", output_directory)

    # Read SRA IDs from input file
    with open(input_file, "r") as f:
        sra_ids = [line.strip() for line in f]

    # Remove stale lock files from output directory
    remove_stale_lock_files(output_directory)
    
    # Clear SRA cache directory
    clear_sra_cache(sra_ids, sra_cache_dir)

    # Generate kallisto index file if it doesn't exist
    genome_index = os.path.splitext(cds)[0] + ".idx"
    print(f"Checking if index file exists at: {genome_index}")

    # Ensure the index file exists and isn't being overwritten
    if os.path.isfile(genome_index):
         print(f"Index file already exists: {genome_index}, skipping kallisto index creation.")
    else:
        try:
            print(f"Index file not found, generating kallisto index at: {genome_index}")
            cmd_index = ["kallisto", "index", "-i", genome_index, cds, "--make-unique"]
            subprocess.run(cmd_index, check=True)
            print(f"Kallisto index file generated: {genome_index}")
        except subprocess.CalledProcessError as e:
            print(f"Error generating kallisto index file: {e}")
            return

    # Process each SRA ID using multiprocessing Pool
    with Pool(num_process) as pool:
        pool.starmap(download_and_convert_SRA, [(sra_id, output_directory, genome_index, sra_cache_dir) for sra_id in sra_ids])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download, convert, align, and delete SRA files.")
    parser.add_argument("--in", dest="in_file", required=True, help="Input file containing SRA IDs, one per line")
    parser.add_argument("--out", dest="out_dir", required=True, help="Output directory to store alignment results")
    parser.add_argument("--n", dest="num_cores", type=int, help="Number of threads being used, will be rounded down to multiple of 4, default = half of avail threads")
    parser.add_argument("--cds", dest="cds", required=True, help="Path to cds file to generate kallisto index (.idx) for read mapping. Can rerun with existent index")
    parser.add_argument("--sra_cache", dest="sra_cache_dir", default="~/.ncbi/public/sra", help="Path to SRA cache directory")
    args = parser.parse_args()

    main(args)

