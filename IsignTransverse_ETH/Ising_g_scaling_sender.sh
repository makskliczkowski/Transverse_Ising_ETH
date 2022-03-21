#!/bin/bash

h=$1	
fun=$2
num=2	 #minimum number of required input
echo $1, $2, $3, $4, $5
operator=0; site=0;
name=g_scaling;
if [[ $# -lt 2 ]]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: h,  fun, operator, site"
  #exit 1
elif [[ $# == 3 ]]; then
  operator=$3;
elif [[ $# == 4 ]]; then
  operator=$3; site=$4;  
elif [[ $# == 5 ]]; then
  operator=$3; site=$4;  name=$5
else
  operator=0; site=0; name=g_scaling;
fi
echo $h, $fun, $operator, $site, $name

#SBATCH --output=logs/g_scale_log-%j-%a.out
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=2G --time=3:00:00 --array=0-19%10 Ising_g_scaling.sh 8 $h 0.05 $fun $operator $site 8
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=2G --time=3:00:00 --array=0-19%10 Ising_g_scaling.sh 9 $h 0.05 $fun $operator $site 8
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-19%10 Ising_g_scaling.sh 10 $h 0.05 $fun $operator $site 8
sbatch --constraint=rack-6 --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=16 --mem=8G --time=3:00:00 --array=0-19%7 Ising_g_scaling.sh 11 $h 0.05 $fun $operator $site 16
sbatch --constraint=rack-6 --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=16 --mem=8G --time=3:00:00 --array=0-19%7 Ising_g_scaling.sh 12 $h 0.05 $fun $operator $site 16
sbatch --constraint=rack-6 --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=16 --mem=16G --time=48:00:00 --array=0-19%7 Ising_g_scaling.sh 13 $h 0.05 $fun $operator $site 16
sbatch --constraint=rack-6 --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=16 --mem=16G --time=48:00:00 --array=0-19%7 Ising_g_scaling.sh 14 $h 0.05 $fun $operator $site 16
#sbatch --constraint=rack-6 --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=32 --mem=80G --time=48:00:00 --array=0-19%10 Ising_g_scaling.sh 15 $h 0.05 $fun $operator $site 32
#sbatch --constraint=rack-6 --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=40 --mem=160G --time=168:00:00 --array=0-19%10 Ising_g_scaling.sh 16 $h 0.05 $fun $operator $site 40

#sbatch --mem=8G --time=3:00:00 --array=0-29%15 Ising_g_scaling.sh 10 $h 0.05 $fun
#sbatch --mem=16G --time=24:00:00 --array=0-29%15 Ising_g_scaling.sh 11 $h 0.05 $fun
#sbatch --mem=16G --time=24:00:00 --array=0-29%15 Ising_g_scaling.sh 12 $h 0.05 $fun
#sbatch --mem=20G --time=24:00:00 --array=0-29%15 Ising_g_scaling.sh 13 $h 0.05 $fun
#sbatch --mem=20G --time=96:00:00 --array=0-29%15 Ising_g_scaling.sh 14 $h 0.05 $fun
#sbatch --mem=30G --time=168:00:00 --array=0-29%15 Ising_g_scaling.sh 15 $h 0.05 $fun
#sbatch --mem=120G --time=168:00:00 --array=0-29%15 Ising_g_scaling.sh 16 $h 0.05 $fun
