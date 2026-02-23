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

for ((i=0; i<5; i++))
do
    echo "Processing index $i"
    mkdir -p noisedata/noiselist_$i
	./noisemap $i
done
