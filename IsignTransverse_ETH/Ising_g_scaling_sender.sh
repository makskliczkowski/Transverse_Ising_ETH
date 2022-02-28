#!/bin/bash

h=$1	
fun=$2
num=2	 #minimum number of required input
echo $fun, $1
operator=0; site=0
if [[ $# -lt 2 ]]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: h,  fun, operator, site"
  exit 1
elif [[ $# == 3 ]]; then
  operator=$3; site=0
else
  operator=$3;  site=$4;  
fi

sbatch --mem=8G --time=33:00:00 --array=0-19%7 Ising_g_scaling.sh 8 $h 0.05 $fun
sbatch --mem=8G --time=33:00:00 --array=0-19%7 Ising_g_scaling.sh 9 $h 0.05 $fun
sbatch --mem=8G --time=33:00:00 --array=0-19%7 Ising_g_scaling.sh 10 $h 0.05 $fun
sbatch --mem=8G --time=33:00:00 --array=0-19%7 Ising_g_scaling.sh 11 $h 0.05 $fun
sbatch --mem=8G --time=33:00:00 --array=0-19%7 Ising_g_scaling.sh 12 $h 0.05 $fun
sbatch --mem=16G --time=48:00:00 --array=0-19%7 Ising_g_scaling.sh 13 $h 0.05 $fun
sbatch --mem=20G --time=48:00:00 --array=0-19%7 Ising_g_scaling.sh 14 $h 0.05 $fun
#sbatch --mem=40G --time=48:00:00 --array=0-19%7 Ising_g_scaling.sh 15 $h 0.05 $fun
#sbatch --mem=90G --time=168:00:00 --array=0-19%7 Ising_g_scaling.sh 16 $h 0.05 $fun

#sbatch --mem=8G --time=3:00:00 --array=0-29%15 Ising_g_scaling.sh 10 $h 0.05 $fun
#sbatch --mem=16G --time=24:00:00 --array=0-29%15 Ising_g_scaling.sh 11 $h 0.05 $fun
#sbatch --mem=16G --time=24:00:00 --array=0-29%15 Ising_g_scaling.sh 12 $h 0.05 $fun
#sbatch --mem=20G --time=24:00:00 --array=0-29%15 Ising_g_scaling.sh 13 $h 0.05 $fun
#sbatch --mem=20G --time=96:00:00 --array=0-29%15 Ising_g_scaling.sh 14 $h 0.05 $fun
#sbatch --mem=30G --time=168:00:00 --array=0-29%15 Ising_g_scaling.sh 15 $h 0.05 $fun
#sbatch --mem=120G --time=168:00:00 --array=0-29%15 Ising_g_scaling.sh 16 $h 0.05 $fun