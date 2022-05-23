#!/bin/bash

dh=$1	
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
echo $dh, $fun, $operator, $site, $name, "${@:6}"
#sbatch --job-name="${name}_L=9" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-0 Ising_h_and_g_scaling.sh 9 $dh 0.05 $fun $operator $site 1 "${@:6}"
#sbatch --job-name="${name}_L=10" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-0 Ising_h_and_g_scaling.sh 10 $dh 0.05 $fun $operator $site 1 "${@:6}"
#sbatch --job-name="${name}_L=11" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-0 Ising_h_and_g_scaling.sh 11 $dh 0.05 $fun $operator $site 1 "${@:6}"
#sbatch --job-name="${name}_L=12" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-0 Ising_h_and_g_scaling.sh 12 $dh 0.05 $fun $operator $site 1 "${@:6}"
#sbatch --job-name="${name}_L=13" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-0 Ising_h_and_g_scaling.sh 13 $dh 0.05 $fun $operator $site 1 "${@:6}"
#exit;
sbatch --job-name="${name}_L=9" --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=2G --time=3:00:00 --array=0-29%2 Ising_h_and_g_scaling.sh 9 $dh 0.05 $fun $operator $site 8 "${@:6}"
sbatch --job-name="${name}_L=10" --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-29%2 Ising_h_and_g_scaling.sh 10 $dh 0.05 $fun $operator $site 8 "${@:6}"
sbatch --job-name="${name}_L=11" --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-29%2 Ising_h_and_g_scaling.sh 11 $dh 0.05 $fun $operator $site 8 "${@:6}"
sbatch --job-name="${name}_L=12" --output=logs/${name}_log-%A-%a.out --cpus-per-task=16 --mem=20G --time=24:00:00 --array=0-29%5 Ising_h_and_g_scaling.sh 12 $dh 0.05 $fun $operator $site 16 "${@:6}"
#sbatch --job-name="${name}_L=13" --output=logs/${name}_log-%A-%a.out --cpus-per-task=32 --mem=40G --time=168:00:00 --array=0-29%2 Ising_h_and_g_scaling.sh 13 $dh 0.05 $fun $operator $site 32 "${@:6}"
#sbatch --job-name="${name}_L=14" --output=logs/${name}_log-%A-%a.out --cpus-per-task=64 --mem=100G --time=168:00:00 --array=0-29%2 Ising_h_and_g_scaling.sh 14 $dh 0.05 $fun $operator $site 64 "${@:6}"
#sbatch --job-name="${name}_L=15" --output=logs/${name}_log-%A-%a.out --cpus-per-task=96 --mem=200G --time=300:00:00 --array=0-29%5 Ising_h_and_g_scaling.sh 15 $dh 0.05 $fun $operator $site 96 "${@:6}"
#sbatch --job-name="${name}_L=16" --output=logs/${name}_log-%A-%a.out --cpus-per-task=32 --mem=70 --time=300:00:00 --array=0-29%10 Ising_h_and_g_scaling.sh 16 $dh 0.05 $fun $operator $site 12 "${@:6}"
exit;
# #SBATCH --output=logs/g_scale_log-%j-%a.out
#sbatch --job-name="${name}_L=8" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-899%5 Ising_h_and_g_scaling.sh 8 $dh 0.05 $fun $operator $site 1 "${@:6}"
sbatch --job-name="${name}_L=9" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-899%5 Ising_h_and_g_scaling.sh 9 $dh 0.05 $fun $operator $site 1 "${@:6}"
sbatch --job-name="${name}_L=10" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-899%5 Ising_h_and_g_scaling.sh 10 $dh 0.05 $fun $operator $site 1 "${@:6}"
sbatch --job-name="${name}_L=11" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-899%5 Ising_h_and_g_scaling.sh 11 $dh 0.05 $fun $operator $site 1 "${@:6}"
sbatch --job-name="${name}_L=12" --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=2G --time=3:00:00 --array=0-899%5 Ising_h_and_g_scaling.sh 12 $dh 0.05 $fun $operator $site 4 "${@:6}"
sbatch --job-name="${name}_L=13" --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=8G --time=24:00:00 --array=0-899%2 Ising_h_and_g_scaling.sh 13 $dh 0.05 $fun $operator $site 4 "${@:6}"
sbatch --job-name="${name}_L=14" --output=logs/${name}_log-%A-%a.out --cpus-per-task=8 --mem=16G --time=168:00:00 --array=0-899%2 Ising_h_and_g_scaling.sh 14 $dh 0.05 $fun $operator $site 8 "${@:6}"
sbatch --job-name="${name}_L=15" --output=logs/${name}_log-%A-%a.out --cpus-per-task=12 --mem=50G --time=168:00:00 --array=0-899%1 Ising_h_and_g_scaling.sh 15 $dh 0.05 $fun $operator $site 12 "${@:6}"
sbatch --job-name="${name}_L=16" --output=logs/${name}_log-%A-%a.out --cpus-per-task=16 --mem=80G --time=200:00:00 --array=0-899%1 Ising_h_and_g_scaling.sh 16 $dh 0.05 $fun $operator $site 16 "${@:6}"
