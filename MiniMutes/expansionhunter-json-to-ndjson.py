# Shoutout Rashin Alabri for the original script
# From https://gist.github.com/rashidalabri/5218b6510376e78c0a006cbf9b563b60

#!/usr/bin/env python3
# SBATCH --job-name=ExpansionHunter JSON to NDJSON
# SBATCH --partition=<insert partition name>
# SBATCH --mail-type=END,FAIL
# SBATCH --mail-user=<insert email>
# SBATCH --cpus-per-task=16
# SBATCH --mem=16G
# SBATCH --time=24:00:00

# --- Configuration

# Adjust this to the number of cores available
CORES = 4

# Directory where NDJSONs should be placed
OUTPUT_DIR = "data/ndjson/"

# Wildcard path to ExpansionHunter JSON outputs. Each file's name is assumed to be the sample name
INPUT_FILES = "exams/*.json"

# Path to JQ command line tool binary
# Leave as the default if the tool is available in your environment
# Note: in HPC systems where you need to run `module load jq` to access the binary, simply
# load it into your environment and then use the output of `which jq` as the path to the binary
JQ_BIN = "jq"

# --- End

import os
from glob import glob
from pathlib import Path


def convert(sample_id, in_path, out_path):
    """Runs the JQ command line tool to convert a single JSON to NDJSON.
    The JQ syntax first accesses the `LocusResults` key in the output.
    Then, filters out the `Variants`. It includes only the sample ID,
    genotype, motif and reference region for each variant. Each entry
    in the ND-JSON file is a JSON object suffixed by a new line--hence
    the newline delimited format.
    """
    cmd = f"cat {in_path} | {JQ_BIN} -c '.LocusResults[] | .Variants[] | {{ sample: \"{sample_id}\", genotype: .Genotype, motif: .RepeatUnit, region: .ReferenceRegion }}' > {out_path}"
    return os.system(cmd)


if __name__ == "__main__":
    # first, find all samples present in the input directory
    samples = []
    for p in glob(INPUT_FILES):
        # extract sample name from path
        sample_id = p.split("/")[-1].split(".")[0]
        samples.append((sample_id, p))

    # for each core, chunk the list of inputs and run the conversion
    chunk_size = len(samples) // CORES + 1
    for i in range(CORES):
        pid = os.fork()
        if pid > 0:
            # we are in the parent process, proceed to launch next process
            continue
        else:
            # we are in the child process, start converting inputs
            start = i * chunk_size
            batch = samples[start : start + chunk_size]

            print(
                f"[Core {i + 1}] Beginning to process {len(batch)} samples", flush=True
            )

            for j, (sample_id, in_path) in enumerate(batch):
                if j % 10 == 0:
                    print(f"[Core {i + 1}] Progress {j}/{len(batch)} samples", flush=True)
                out_path = Path(OUTPUT_DIR) / (sample_id + ".ndjson")
                convert(sample_id, in_path, str(out_path))

            # when finished, exit the child process
            print(f"[Core {i + 1}] Finished", flush=True)
            exit()

    # make the parent process wait for all child proccesses
    for _ in range(CORES):
        os.wait()

    print(f"All input files converted. Exiting.", flush=True)
