#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q
#$ -N chaotic
#$ -o job_out
#$ -e job_out
#$ -pe OpenMP 14


MODEL=$(grep '^#define MODEL' model.hpp | awk '{print $3}')

data=$(grep 'const std::string sdatadir' parameters.hpp \
        | sed 's|//.*||' \
        | sed 's/.*= *"\(.*\)".*/\1/')

if [ "$MODEL" -eq 0 ]; then
    model="chaotic"
elif [ "$MODEL" -eq 1 ]; then
    model="Starobinsky"
elif [ "$MODEL" -eq 2 ]; then
    model="USR"
elif [ "$MODEL" -eq 3 ]; then
    model="hybrid"
else
    echo "Unknown MODEL=$MODEL"
    exit 1
fi

mkdir -p "$data"
mkdir -p "$data/$model"
mkdir -p "$data/$model/animation"

ll=0
numdata=100
for ((i=ll; i<ll+numdata; i++))
do
for m in 100.
do
   echo "calPzeta=$m"
   SIGMA=$m ./STOLAS $i
done
done
