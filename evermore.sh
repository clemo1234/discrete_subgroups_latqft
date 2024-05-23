#!/bin/bash

beta1=$1
beta2=$2
seed=$RANDOM
D=$3
nt=$4
nx=$5
group=$4
fol=$7

export OMP_NUM_THREADS=4
export GOMP_CPU_AFFINITY=0-3

out=${fol}out_b${beta1}_g${group}_nt${nt}_nx${nx}_s${seed}.log

echo "./dym-mod-metro.cpp ./groups/my$group $D $nt $nx 0 $beta1 $beta2 $seed" > $out
date >> $out
./dym-mod-metro ./groups/$group $Z $nt $nx 0 $beta1 0 $seed >> $out
#./dym-mod-metro ./groups/my$group $D $nt $nx 0 $beta1 0 $seed >> $out
#./dym-mod-metro-savecfg ./groups/mys1080-v4 4 $nt $nx 0 $beta1 $beta2 $seed >> $out
date >> $out
