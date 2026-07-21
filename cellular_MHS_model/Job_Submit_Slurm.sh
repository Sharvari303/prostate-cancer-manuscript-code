#!/bin/bash
#SBATCH --job-name=PCaHybridModel
#SBATCH --output=slurm_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=20:00:00
#SBATCH --partition=RM-shared
#SBATCH --account=mcb200052p

# Load Python 2.7 module
module load python/2.7.18

# Activate virtual environment
source /ocean/projects/mcb200052p/sskemkar/PCa_Hybrid_Model/venv_pca/bin/activate

# Add CopasiSE to PATH
export PATH=/ocean/projects/mcb200052p/sskemkar/PCa_Hybrid_Model/Copasi_Executables/Linux_64Bit:$PATH

cd /ocean/projects/mcb200052p/sskemkar/PCa_Hybrid_Model

echo "Job started on $(hostname) at $(date)"
echo " "

python combined_ode_boolean.py -p Model_Input.json

echo " "
echo "Job ended on $(hostname) at $(date)"
