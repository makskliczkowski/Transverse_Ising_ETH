#!/bin/bash
#SBATCH --job-name=spectrals
#SBATCH --output=logs/spectral_log-%A-%a.out
#SBATCH --mem-bind=verbose,local

#INPUTS TO SCRIPT ARE: L,  h/g,  operator, dg/dh, thread_num

num=5
if [[ $# != $num ]]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: L,  h/g,  operator, dg/dh, thread_num"
  exit 1
fi
thread_num=$5
# set number of realisations, array from L=8 to L=16
	R_ARR=(20 20 20 20 20 10 10 10 10)
	r=${R_ARR[`expr ${1}-8`]}

# set helper vars 
s_max=0;
if [[ $3 -gt 2 ]]; then	
	echo "chose k-space operators"
	s_max=$(echo $1 | awk '{printf "%d",  $1 / 2. + $1 % 2}')
else
	echo $3
	s_max=$1
	echo "chose local operators"
fi

module purge
export OMP_NUM_THREADS=$thread_num

module load OpenBLAS/0.3.18-GCC-11.2.0
module load HDF5/1.12.0-gompi-2021a
module load intel/2022.00

echo "smax=${s_max}"

h_vs_g=1
# set input parameters: to multiply use \*, cause * means 'all files'
# also to have floating-point use =$echo("scale=num_of_digits; {expresion}" | bc)
	site=`expr $SLURM_ARRAY_TASK_ID % $s_max`
	step=`expr $SLURM_ARRAY_TASK_ID / $s_max`
if [[ h_vs_g == 0 ]]; then
	g=$(echo $4 $step | awk '{printf "%.3f", $1 + $1 * $2}')
	h=$2
else
	h=$(echo $4 $step | awk '{printf "%.3f", $1 + $1 * $2}')
	g=$2
fi;

# output filename
opName=0
if [[ $3 == 0 ]]; then
	opName="SigmaZ_j=${site}_spectral"
elif [[ $3 == 1 ]]; then
	opName="SigmaX_j=${site}_spectral"
elif [[ $3 == 2 ]]; then
	opName="H_j=${site}_spectral"
elif [[ $3 == 3 ]]; then
	opName="SigmaZ_q=${site}_spectral"
elif [[ $3 == 4 ]]; then
	opName="SigmaX_q=${site}_spectral"
elif [[ $3 == 5 ]]; then
	opName="H_q=${site}_spectral"
elif [[ $3 == 6 ]]; then
	opName="TFIM_LIOM_plus_n=${site}_spectral"
elif [[ $3 == 7 ]]; then
	opName="TFIM_LIOM_minus_n=${site}_spectral"
else 
	opName="damn"
fi

	filename="run_logs/${opName}_L=${1}_h=${h}_g=${g}"

#print all variables to see if all correct
	echo "L=${1}"
	echo "h=${h}"
	echo "g=${g}"
	echo "site=${site}"
	echo "operator=${3}"
	echo "name=${filename}"
	echo "realisations=${r}"

./Ising.o -L $1 -g $g -h $h -w 0.01 -th $thread_num -m 0 -r $r -op $3 -fun 1 -s $site -b 0 >& ${filename}.log