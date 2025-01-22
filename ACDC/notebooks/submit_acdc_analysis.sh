#!/bin/bash
#SBATCH --job-name=acdc_analysis      # Job name
#SBATCH --output=acdc_analysis.out    # Output log
#SBATCH --error=acdc_analysis.err     # Error log
#SBATCH --mail-type=ALL               # Mail notifications (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mahima.arunkumar.1@gmail.com  
#SBATCH --partition=lrz-dgx-a100-80x8
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Number of tasks (use 1 task)
#SBATCH --cpus-per-task=16            # Number of CPU cores 
#SBATCH --gres=gpu:1                  # Request 1 GPU
#SBATCH --mem=128G                    # Memory allocation 
#SBATCH --time=48:00:00               # Time limit

# Print a message to indicate the job has started
echo "Starting the SLURM job: $(date)"
 

# Activate Python environment
module load python/3.10 
source ~/envs/acdc/bin/activate

#Set up PATH
export PATH=$PATH:/dss/dsshome1/0F/di93quv/.local/bin
export PYTHONPATH=$PYTHONPATH:/dss/dsshome1/0F/di93quv/Systems_biomedicine/acdc

# Load required modules
module load cuda/11.4   # Load CUDA if using GPU
module load gcc/10.2    # Load GCC for compatibility if needed

# Print a message after loading modules
echo "Modules loaded: $(date)"
echo "Installing Python dependencies with pip..."

# Install libraries
pip install --user pandas numpy scikit-learn scipy seaborn matplotlib phenograph umap-learn scanpy anndata joblib plotly

# Install ACDC package locally from its path
echo "Installing ACDC package..."
pip install -e /dss/dsshome1/0F/di93quv/Systems_biomedicine/acdc
echo "Dependencies installed: $(date)"


# Run Python script
echo "Running Python script..."

# Run the script
python ./run_acdc_analysis.py
echo "Python script finished: $(date)"
