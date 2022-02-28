#!/bin/bash
#SBATCH --job-name=g_scaling
#SBATCH --output=logs/g_scale_log-%j-%a.out
#SBATCH --cpus-per-task=16

module purge
#unset MKL_SERIAL
#export CPATH="/usr/include/hdf5/serial"
#export LD_LIBRARY_PATH="${MKLROOT}/lib/intel64"
export OMP_NUM_THREADS=16

#module load Armadillo/9.900.1-foss-2020a
module load OpenBLAS/0.3.18-GCC-11.2.0
module load foss/2021b
#module load intel/2021a

#INPUTS TO SCRIPT ARE: L,  h,  dg, fun, operator , site

operator=0; site=0;
num=4	 #minimum number of required input
suffix=""
if [[ $# -lt 4 ]]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: L,  h, dg, fun"
  exit 1
elif [[ $# -eq 4 ]]; then
  operator=0; site=0;
elif [[ $# -eq 5 ]]; then
  operator=$5; site=0;
  suffix="_op=${5}";
else
  operator=$5;  site=$6;  
  suffix="_op=${5}_site=${6}";
fi
# set number of realisations, array from L=8 to L=16
	R_ARR=(1000 1000 1000 1000 1000 500 200 50 10)
	r=${R_ARR[`expr ${1}-8`]}

# set input parameters: to multiply use \*, cause * means 'all files'
# also to have floating-point use =$echo("scale=num_of_digits; {expresion}" | bc)
	g=$(echo $3 $SLURM_ARRAY_TASK_ID | awk '{printf "%.3f", $1 + $1 * $2}')

# output filename
funName=""
if [[ $4 -eq 0 ]]; then
	funName="Spectrals"
elif [[ $4 -eq 1 ]]; then
	funName="AGP_vg"
elif [[ $4 -eq 2 ]]; then
	funName="AGP_vh"
elif [[ $4 -eq 3 ]]; then
	funName="TFIM_LIOMs"
elif [[ $4 -eq 4 ]]; then
	funName="RelaxationTimes"
else
	funName="Entropy"
fi

	filename="run_logs/${funName}_L=${1}_h=${2}_g=${g}${suffix}"

#print all variables to see if all correct
	echo "L=${1}"
	echo "h=${2}"
	echo "g=${g}"
	echo "site=${site}"
	echo "operator=${operator}"
	echo "name=${filename}"
	echo "realisations=${r}"
#g++ -std=c++2a main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising_${funName}_L=${1}_h=${2}_g=${g}${suffix}.o\
# -DARMA_DONT_USE_WRAPPER -I/home/rswietek/LIBRARIES_CPP/armadillo-10.8.2/include -llapack -lopenblas -fopenmp -lpthread -lm -lstdc++fs -fomit-frame-pointer -Ofast >& compile_${funName}_L=${1}_h=${2}_g=${g}${suffix}.log
 
./Ising.o -L $1 -g $g -h $2 -w 0.01 -th 16 -m 0 -w 0.01 -r $r -op $operator -fun $4 -s $site -b 0 >& ${filename}.log