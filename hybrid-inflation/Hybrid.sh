#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q
#$ -N chaotic
#$ -o job_out
#$ -e job_out
#$ -pe OpenMP 14

#source /opt/intel/bin/compilervars.sh intel64
#export OMP_NUM_THREADS=$NSLOTS

data=$(grep 'const std::string sdatadir' parameters.hpp | sed 's|//.*||' | sed 's/.*= *"\(.*\)".*/\1/')
model=$(grep 'const std::string model' model.hpp | sed 's|//.*||' | sed 's/.*= *"\(.*\)".*/\1/')

mkdir -p "$data"
mkdir -p "$data""/$model"
mkdir -p "$data""/$model/animation"

for ((i=0; i<1; i++))
do
    echo "Number $i"
	./STOLAS $i
done
