#!/bin/bash

# In terminal: bash submit_job.sh

#Define log-file
log_file="/home/pcpt3/HXLinker/HX_linkers.log"
exec > "${log_file}" 2>&1

#Run first python script to create folders and files
echo "Running python3 orca_eval.py"
source activate env4

INPUT="/home/pcpt3/HXLinker/conformers.xyz"
ARRAY_SIZE=1296
CHRG=0 
OUTPUT_FOLDER="/rds/user/pcpt3/hpc-work/HXLinkers/runall"

# Define the range of indices and the base directory
BASE_DIR="$OUTPUT_FOLDER/ORCA_job_"
MAX_INDEX=1296  # Adjust this to the maximum job number
MISSING_ARRAY=()

# Loop through all possible indices to check whether it has already been run before
# This prevents calculating some again 
for i in $(seq 0 $((MAX_INDEX - 1))); do
    FILE_PATH="${BASE_DIR}${i}/atoms-dft.xyz"
    if [ ! -f "$FILE_PATH" ]; then
        MISSING_ARRAY+=($i)  # Add the missing index to the array
    fi
done

# Create the SLURM array specification
if [ ${#MISSING_ARRAY[@]} -gt 0 ]; then
    ARRAY_SIZE=${#MISSING_ARRAY[@]}
    SLURM_ARRAY="--array=$(IFS=,; echo "${MISSING_ARRAY[*]}")"
    echo "SLURM array: $SLURM_ARRAY"
else
    echo "No missing indices found."
fi

#sbatch $SLURM_ARRAY  --export=INPUT=${INPUT},CHRG=${CHRG},OUTPUT_FOLDER=${OUTPUT_FOLDER} orca_eval.peta4-icelake
