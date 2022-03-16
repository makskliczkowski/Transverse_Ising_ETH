#!/bin/bash
#SBATCH --job-name=AGP
#SBATCH --output=logs/log-%j-%a.out
#SBATCH --error=errors/log-%j-%a.err
#SBATCH --time=168:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=10

export MKL_SERIAL=yes
export CPATH="/usr/include/hdf5/serial"
export LD_LIBRARY_PATH="${MKLROOT}/lib/intel64"

module purge
module load Armadillo/9.900.1-foss-2020a
module load intel/2021b
module load imkl/2021.4.0
module load OpenBLAS/0.3.18-GCC-11.2.0

#INPUTS TO SCRIPT ARE: operator, dg, h_min

num=3
if [ $# != $num ]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: operator, dg, h_min"
  exit 1
fi

# set number of realisations, array from L=8 to L=16
	R_ARR=(10 10 10 10 10 10 10 1 1)
	r=${R_ARR[`expr ${1}-8`]}

# set helper vars 
s_max=7

# set input parameters: to multiply use \*, cause * means 'all files'
# also to have floating-point use =$echo("scale=num_of_digits; {expresion}" | bc)
	site=`expr $SLURM_ARRAY_TASK_ID % $s_max`
	step=`expr $SLURM_ARRAY_TASK_ID / $s_max`
	g=$(echo $2 $step | awk '{printf "%.3f", $1 + $1 * $2}')

# output filename
opName=""
if [ $1 == 0 ]; then
	opName="AGP_SigmaZ_j=${site}"
elif [ $1 == 1 ]; then
	opName="AGP_SigmaX_j=${site}"
elif [ $1 == 2 ]; then
	opName="AGP_SigmaZ_q=${site}"
elif [ $1 == 3 ]; then
	opName="AGP_SigmaX_q=${site}"
else
	opName="AGP_H_q=${site}"
fi

	filename="./run_logs/${opName}_g=${g}_hmin=${3}"

#print all variables to see if all correct
	echo "h_min=${3}"
	echo "g=${g}"
	echo "site=${site}"
	echo "operator=${1}"
	echo "name=${filename}"

./Ising.o -L 8 -Ls 1 -Ln 7 -g $g -h $3 -hs 0.01 -hn 300 -w 0.01 -th 10 -m 0 -r 1 -op $1 -fun 1 -s $site -b 0 >& ${filename}.log