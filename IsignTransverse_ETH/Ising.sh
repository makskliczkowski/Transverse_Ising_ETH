#!/bin/bash
#SBATCH --time=167:59:59
#SBATCH --mem=60G
#SBATCH --cpus-per-task=24
#SBATCH --job-name=AGP

export MKL_SERIAL=yes
export CPATH="/usr/include/hdf5/serial"
export LD_LIBRARY_PATH="${MKLROOT}/lib/intel64"

module purge
module load Armadillo/9.900.1-foss-2020a
module load intel/2021b

./Ising.o -L 8 -Ls 1 -Ln 7 -g 0.025 -gn 3 -gs 0.025 -h 0.2 -hn 15 -hs 0.2 -th 24 -m 0 -r 1 -k 0 -p 1 -x 1 -b 0 ./results/>& log_AGP.txt