#!/bin/bash
#SBATCH --output=logs/h_and_g_scale_log-%j-%a.out

#INPUTS TO SCRIPT ARE: L,  dh,  dg, fun, operator , site, thread_num

#set boolean value
ch=1


operator=$5; site=$6; thread_num=$7;
num=4	 #minimum number of required input
suffix="_op=${5}_site=${6}_id=${SLURM_JOB_ID}";
if [[ $# -lt 4 ]]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: L,  dh, dg, fun"
  exit 1
fi
module purge
export OMP_NUM_THREADS=$thread_num

module load OpenBLAS/0.3.18-GCC-11.2.0
module load HDF5/1.12.0-gompi-2021a
module load intel/2022.00

# set number of realisations, array from L=8 to L=16
	#R_ARR=(600 600 400 400 300 200 200 100 100)										# ch = 0
	#R_ARR=(600 600 400 400 300 300 200 100 100 100 100 100 50 50 50 50 50) 	# ch = 1
	R_ARR=(200 200 200 200 200 100 100 50 50) # diagonalizaiton
	r=${R_ARR[`expr ${1}-8`]}

# set input parameters: to multiply use \*, cause * means 'all files'
# also to have floating-point use =$echo("scale=num_of_digits; {expresion}" | bc)
	_max=30 # max number of points on g-axis (and h-axis)
    g_step=`expr $SLURM_ARRAY_TASK_ID / $_max`
    g=$(echo $2 $g_step | awk '{printf "%.3f", $1 + $1 * $2}')

    h_step=`expr $SLURM_ARRAY_TASK_ID % $_max`
    h=$(echo $3 $h_step | awk '{printf "%.3f", $1 + $1 * $2}')
    
# output filename
funName=""
if [[ $4 -eq 0 ]]; then
	funName="Diagonalize"
elif [[ $4 -eq 1 ]]; then
	funName="Spectrals"
elif [[ $4 -eq 2 ]]; then
	funName="EntropyEvolution"
	if [[ $ch == 1 ]]; then
		funName="${funName}_lanczos";
	fi
elif [[ $4 -eq 3 ]]; then
	funName="SFF"
elif [[ $4 -eq 4 ]]; then
	funName="RelaxationTimes"
elif [[ $4 -eq 5 ]]; then
	funName="Benchmark"
elif [[ $4 -eq 6 ]]; then
	funName="AGP"
else
	funName="Other"
fi
filename="run_logs/${funName}_L=${1}_h=${h}_g=${g}${suffix}"
#print all variables to see if all correct
	echo "L=${1}"
	echo "h=${h}"
	echo "g=${g}"
	echo "site=${site}"
	echo "operator=${operator}"
	echo "name=${filename}"
	echo "realisations=${r}"
#g++ -std=c++2a main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising_${funName}_L=${1}_h=${2}_g=${g}${suffix}.o\
# -DARMA_DONT_USE_WRAPPER -I/home/rswietek/LIBRARIES_CPP/armadillo-10.8.2/include -llapack -lopenblas -fopenmp -lpthread -lm -lstdc++fs -fomit-frame-pointer -Ofast >& compile_${funName}_L=${1}_h=${2}_g=${g}${suffix}.log
 
./Ising.o -L $1 -g $g -h $h -th $thread_num -m 0 -w 0.1 -r $r\
 -op $operator -fun $4 -s $site -b 0 -ch $ch "${@:8}" >& ${filename}.log
