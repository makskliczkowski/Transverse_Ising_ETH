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
else
  operator=$3; site=$4;  name=$5
fi
echo $h, $fun, $operator, $site, $name, "${@:6}"

sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=1G --time=3:00:00 Ising_g_scaling.sh 10 $h 0.0 $fun $operator $site 8 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=1G --time=3:00:00 Ising_g_scaling.sh 12 $h 0.0 $fun $operator $site 8 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=2G --time=3:00:00 Ising_g_scaling.sh 14 $h 0.0 $fun $operator $site 8 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=4G --time=3:00:00 Ising_g_scaling.sh 16 $h 0.0 $fun $operator $site 24 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=8G --time=96:00:00 Ising_g_scaling.sh 18 $h 0.0 $fun $operator $site 24 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=16G --time=96:00:00 Ising_g_scaling.sh 20 $h 0.0 $fun $operator $site 40 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=16G --time=96:00:00 Ising_g_scaling.sh 22 $h 0.0 $fun $operator $site 40 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=32G --time=96:00:00 Ising_g_scaling.sh 24 $h 0.0 $fun $operator $site 40 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=60G --time=96:00:00 Ising_g_scaling.sh 26 $h 0.0 $fun $operator $site 40 "${@:6}"

exit;
#SBATCH --output=logs/g_scale_log-%j-%a.out
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=2G --time=3:00:00 --array=0-19%10 Ising_g_scaling.sh 8 $h 0.05 $fun $operator $site 8 "${@:6}"
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=2G --time=3:00:00 --array=0-19%10 Ising_g_scaling.sh 9 $h 0.05 $fun $operator $site 8 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=2G --time=3:00:00 --array=0-19%5 Ising_g_scaling.sh 10 $h 0.05 $fun $operator $site 4 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-19%5 Ising_g_scaling.sh 11 $h 0.05 $fun $operator $site 8 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=16 --mem=8G --time=48:00:00 --array=0-19%5 Ising_g_scaling.sh 12 $h 0.05 $fun $operator $site 16 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=16 --mem=16G --time=48:00:00 --array=0-19%5 Ising_g_scaling.sh 13 $h 0.05 $fun $operator $site 16 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=24 --mem=16G --time=96:00:00 --array=0-19%5 Ising_g_scaling.sh 14 $h 0.05 $fun $operator $site 24 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=24 --mem=60G --time=96:00:00 --array=0-19%5 Ising_g_scaling.sh 15 $h 0.05 $fun $operator $site 24 "${@:6}"
sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=32 --mem=120G --time=168:00:00 --array=0-19%5 Ising_g_scaling.sh 16 $h 0.05 $fun $operator $site 32 "${@:6}"

#sbatch --mem=8G --time=3:00:00 --array=0-29%15 Ising_g_scaling.sh 10 $h 0.05 $fun
#sbatch --mem=16G --time=24:00:00 --array=0-29%15 Ising_g_scaling.sh 11 $h 0.05 $fun
#sbatch --mem=16G --time=24:00:00 --array=0-29%15 Ising_g_scaling.sh 12 $h 0.05 $fun
#sbatch --mem=20G --time=24:00:00 --array=0-29%15 Ising_g_scaling.sh 13 $h 0.05 $fun
#sbatch --mem=20G --time=96:00:00 --array=0-29%15 Ising_g_scaling.sh 14 $h 0.05 $fun
#sbatch --mem=30G --time=168:00:00 --array=0-29%15 Ising_g_scaling.sh 15 $h 0.05 $fun
#sbatch --mem=120G --time=168:00:00 --array=0-29%15 Ising_g_scaling.sh 16 $h 0.05 $fun
