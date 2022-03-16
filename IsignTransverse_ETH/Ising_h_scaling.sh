#!/bin/bash
#SBATCH --job-name=timeEvol
#SBATCH --output=logs/log-%j-%a.out
#SBATCH --error=errors/log-%j-%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=90G
#SBATCH --cpus-per-task=20

export MKL_SERIAL=yes
export CPATH="/usr/include/hdf5/serial"
export LD_LIBRARY_PATH="${MKLROOT}/lib/intel64"

module purge
module load Armadillo/9.900.1-foss-2020a
module load intel/2021b
module load imkl/2021.4.0
module load OpenBLAS/0.3.18-GCC-11.2.0

#INPUTS TO SCRIPT ARE: L,  h,  operator, dg

num=4
if [ $# != $num ]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: L,  h,  operator, dg"
  exit 1
fi
# set number of realisations, array from L=8 to L=16
	R_ARR=(10 10 10 10 10 10 10 1 1)
	r=${R_ARR[`expr ${1}-8`]}

# set helper vars 
s_max=0;
if [ $3 -gt 1 ]
then	
	echo "chose k-space operators"
	s_max=$(echo $1 | awk '{printf "%d",  $1 / 2. + $1 % 2}')
else
	echo $3
	s_max=$1
	echo "chose local operators"
fi
echo "smax=${s_max}"
# set input parameters: to multiply use \*, cause * means 'all files'
# also to have floating-point use =$echo("scale=num_of_digits; {expresion}" | bc)
	site=`expr $1 / 2`
	step=`expr $SLURM_ARRAY_TASK_ID`
	g=$(echo $4 $step | awk '{printf "%.3f", $1 + $1 * $2}')

# output filename
opName=0
if [ $3 == 0 ]; then
	opName="SigmaZ_j=${site}_spectral"
elif [ $3 == 1 ]; then
	opName="SigmaX_j=${site}_spectral"
elif [ $3 == 2 ]; then
	opName="SigmaZ_q=${site}_spectral"
elif [ $3 == 3 ]; then
	opName="SigmaX_q=${site}_spectral"
else
	opName="H_q=${site}_spectral"
fi

	filename="${opName}_L=${1}_h=${2}_g=${g}"

#print all variables to see if all correct
	echo "L=${1}"
	echo "h=${2}"
	echo "g=${g}"
	echo "site=${site}"
	echo "operator=${3}"
	echo "name=${filename}"
	echo "realisations=${r}"

./Ising.o -L $1 -g $g -h $2 -w 0.01 -th 20 -m 0 -r $r -op $3 -fun 0 -s $site -b 0 >& ${filename}.log