#!/bin/bash
#SBATCH --mail-user=rafal.swietek@ijs.si
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --job-name=timeEvol
#SBATCH --output=logs/log-%j-%a.out
#SBATCH --error=errors/log-%j-%a.err
#SBATCH --time=168:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=20

export MKL_SERIAL=yes
export CPATH="/usr/include/hdf5/serial"
export LD_LIBRARY_PATH="${MKLROOT}/lib/intel64"

module purge
module load Armadillo/9.900.1-foss-2020a
module load intel/2021b

#INPUTS TO SCRIPT ARE: L,  h,  operator

num=3
if [ $# != $num ]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: L,  h,  operator"
  exit 1
fi
# set number of realisations, array from L=8 to L=16
	R_ARR=(10 10 10 10 10 10 10 1 1)
	r=${R_ARR[`expr ${1}-8`]}

# set input parameters: to multiply use \*, cause * means 'all files'
# also to have floating-point use =$echo("scale=num_of_digits; {expresion}" | bc)
	site=`expr $SLURM_ARRAY_TASK_ID % ${1}`
	step=`expr $SLURM_ARRAY_TASK_ID / ${1}`
	g=$(echo 0.2 $step | awk '{printf "%.3f", $1 + $1 * $2}')

# output filename
	filename="timeEvol_L=${1}_h=${2}_g=${g}_site=${site}_op=${3}"

#print all variables to see if all correct
	echo "L=${1}"
	echo "h=${2}"
	echo "g=${g}"
	echo "site=${site}"
	echo "operator=${3}"
	echo "name=${filename}"
	echo "realisations=${r}"

./Ising.o -L $1 -g $g -h $2 -w 0.01 -th 20 -m 0 -r $r -op $3 -fun 0 -s $site -b 0 >& ${filename}.log