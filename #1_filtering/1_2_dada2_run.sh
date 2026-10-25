#!/bin/bash
#SBATCH --job-name=2_dada2
#SBATCH --output=2_dada2_%j.out
#SBATCH --error=2_dada2_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=Node4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G

set -euo pipefail

source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate dada2_env

WORK="/storage2/flower_eDNA/analysis/5_20260819_new_pipe"
cd "$WORK"

echo "[INFO] $(date) starting DADA2"
echo "[INFO] Rscript: $(which Rscript)"

Rscript 2_dada2.R

echo "[DONE] $(date) finished"
