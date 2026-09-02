#!/bin/bash
#SBATCH --job-name=nf-test-suite
#SBATCH --partition=hsph
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=/n/home04/smaharjan/biobakery-nextflow/test/results/slurm_%j.out
#SBATCH --error=/n/home04/smaharjan/biobakery-nextflow/test/results/slurm_%j.err

set -euo pipefail

# Load environment
source /n/lab_storage/huttenhower_lab/tools/hutlab/src/hutlabrc_rocky8.sh
module use /n/lab_storage/huttenhower_lab/tools/hutlab/src/modules_rocky8
module load jdk/21.0.2-fasrc01

# The suite runs about a dozen nextflow drivers at once, each its own JVM,
# which is what the memory and cpu requests above are for. The drivers only
# submit and wait; the real work runs in the SLURM jobs they spawn.

# Confirm nextflow is reachable
/n/lab_storage/huttenhower_lab/tools/nextflow/24.10.4/bin/nextflow -version

# Run the test suite from the repo root
cd /n/home04/smaharjan/biobakery-nextflow
bash test/run_tests.sh
