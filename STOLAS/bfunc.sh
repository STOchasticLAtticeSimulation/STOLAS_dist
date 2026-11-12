#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q
#$ -N chaotic
#$ -o job_out
#$ -e job_out
#$ -pe OpenMP 14

ll=0
numdata=100
for ((i=ll; i<ll+numdata; i++))
do
for m in 2000.
do
   echo "calPzeta=$m"
   SIGMA=$m ./functions $i
done
done
