#!/bin/bash
#SBATCH --output=logs/g_scale_log-%j-%a.out

#--------------------- MODULES
module purge

module load OpenBLAS/0.3.18-GCC-11.2.0
module load HDF5/1.12.0-gompi-2021a
module load intel/2022.00

# to avoid illegal instructions
#--constraint=rack-6

#INPUTS TO SCRIPT ARE: L,  x, y, dz, fun, operator , site, thread_num, jobid

#set boolean value
ch=0
L=$1;
x=$2
y=$3
dz=$4	
fun=$5
operator=$6; 
site=$7; 
thread_num=$8;
jobid=$9;
export OMP_NUM_THREADS=$thread_num

echo $L, $x, $y, $dz, $fun, $operator, $site, $thread_num
num=5	 #minimum number of required input
if [[ $# -lt 5 ]]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: L, x, y, dz, fun"
  exit 1
fi

# set number of realisations, array from L=8 to L=16
	R_ARR=(200 200 200 100 100 40 40 20 10);	r=${R_ARR[`expr ${L}-8`]}
	#R_ARR=(200 200 200 100 100 50 50 20 10);	r=${R_ARR[`expr ${L}-12`]}
	#r = 1;
# set input parameters: to multiply use \*, cause * means 'all files'
# also to have floating-point use =$echo("scale=num_of_digits; {expresion}" | bc)
	#par0=$dz;
	#par0=0.1
	par=$(echo $dz $SLURM_ARRAY_TASK_ID $par0 | awk '{printf "%.3f", $3 + $1 * $2}')
	#g=$par;	h=$x;		J=$y;
	#g=$y;		h=$par;	J=$x;
	#g=$x;		h=$y;		J=$par;
	
	wx=0.3; # -- disorder with longest thouless time in ergodic regime
	g=$x;		h=$y;		wx=$par;	J=0.05;

	seed=$(echo $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT | awk '{printf "%d", $1 + $2 / $3 * $1}')  
#--------------------- output filename
funName=""
if [[ $fun -eq 0 ]]; then
	funName="Diagonalize"
elif [[ $fun -eq 1 ]]; then
	funName="Spectrals"
elif [[ $fun -eq 2 ]]; then
	funName="EntropyEvolution"
	if [[ $ch == 1 ]]; then
		funName="${funName}_lanczos";
	fi
elif [[ $fun -eq 3 ]]; then
	funName="SFF"
elif [[ $fun -eq 4 ]]; then
	funName="RelaxationTimes"
elif [[ $fun -eq 5 ]]; then
	funName="Benchmark"
elif [[ $fun -eq 6 ]]; then
	funName="AGP"
else
	funName="Other"
fi
suffix="_op=${operator}_site=${site}_seed=${seed}";
filename="run_logs/${funName}_L=${L}_J=${J}_h=${h}_g=${g}_w=${wx}${suffix}"

#--------------------- print all variables to see if all correct
	echo "L=${L}"
	echo "J=${J}"
	echo "h=${h}"
	echo "g=${g}"
	echo "site=${site}"
	echo "operator=${operator}"
	echo "name=${filename}"
	echo "realisations=${r}"
	echo "seed=${seed}"
#g++ -std=c++2a main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising_${funName}_L=${1}_h=${2}_g=${g}${suffix}.o\
# -DARMA_DONT_USE_WRAPPER -I/home/rswietek/LIBRARIES_CPP/armadillo-10.8.2/include -llapack -lopenblas -fopenmp -lpthread -lm -lstdc++fs -fomit-frame-pointer -Ofast >& compile_${funName}_L=${1}_h=${2}_g=${g}${suffix}.log
 
 #--------------------- run
./Ising.o "${@:10}" -L $L -J $J -g $g -h $h -th $thread_num -m 0 -w $wx -r $r\
 -op $operator -fun $fun -s $site -b 0 -ch $ch -seed $seed -jobid $jobid >& ${filename}.log
