#!/bin/bash
#SBATCH -J ORCA_FF
#SBATCH -A THEORYRIG-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time 12:00:00
##SBATCH --time 1:00:00
##SBATCH --qos=intr

export OMP_STACKSIZE=512M

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
workdir="$SLURM_SUBMIT_DIR"

source /home/pcpt3/.bashrc
conda init
conda activate env3-p11
python --version

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment
module load python3
#module load openmpi/4.1.1/gcc/mnop75he

module unload -f intel-oneapi-mpi
module load openmpi/4.1.8/intel/xrukf3ol

#export PATH=/rds/user/pcpt3/hpc-work/orca:$PATH
#export LD_LIBRARY_PATH=/rds/user/pcpt3/hpc-work/orca:$LD_LIBRARY_PATH

export PATH=/rds/user/pcpt3/hpc-work/orca_6_1_0:$PATH
export LD_LIBRARY_PATH=/rds/user/pcpt3/hpc-work/orca_6_1_0:$LD_LIBRARY_PATH

#export OPENMPI_HOME=/usr/local/software/spack/spack-views/rhel8-icelake-20211027_2/openmpi-4.1.1/gcc-11.2.0/mnop75hedc4hru2c22d6oqlu7opdc4jg 
export OPENMPI_HOME=/usr/local/software/spack/spack-views/rocky8-icelake-20220710/openmpi-4.1.8/intel-2021.6.0/xrukf3ol2ktl7hd6j5ga6y7ft3sz7o3r
export PATH=$OPENMPI_HOME/bin:$PATH
export LD_LIBRARY_PATH=$OPENMPI_HOME/lib:$LD_LIBRARY_PATH

#Create folder for ORCAjobs
export SLURM_SUBMIT_DIR

dir_path=$(pwd)
#mkdir -p "$dir_path"
#cd "$dir_path"
   
# Creating local scratch folder 
export scratchlocation=/rds/user/pcpt3/hpc-work/Orca-tmp
tdir=$(mktemp -d $scratchlocation/ORCAjob__$SLURM_JOB_ID-XXXX)

# Copy only the necessary stuff in submit directory to scratch directory. Add more (e.g. .gbw) here if needed.
cp  $dir_path/*.inp $tdir/
cp  $dir_path/*.xyz $tdir/

# Creating nodefile in scratch
echo $SLURM_NODELIST > $tdir/$job.nodes

# cd to scratch
cd $tdir

rm -f ${dir_path}/ORCA.out
echo "Job execution start: $(date)" >>  ${dir_path}/ORCA.out
echo "Shared library path: $LD_LIBRARY_PATH" >>  ${dir_path}/ORCA.out
echo "Slurm Job ID is: ${SLURM_JOB_ID}" >>  ${dir_path}/ORCA.out
echo "Slurm Job name is: ${SLURM_JOB_NAME}" >>  ${dir_path}/ORCA.out
echo $SLURM_NODELIST >> ${dir_path}/ORCA.out
echo "Now in directory: $(pwd)"

# Execute the application
/rds/user/pcpt3/hpc-work/orca_6_1_0/orca ORCA.inp > ORCA.out
cp $tdir/*.{xyz,gbw,out,nodes,final,system,new,info,mol} $dir_path
rm -rf $tdir
