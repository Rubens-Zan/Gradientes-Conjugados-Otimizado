#!/bin/bash

METRICA="L3 L2CACHE FLOPS_DP FLOPS_AVX"
TAMANHOS="32 64 128 256 512 1000 2000 4000 8000"

echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

for n in $TAMANHOS
do
    for k in $METRICA
    do
        likwid-perfctr -C 3 -g ${k} -m -O ./cgSolver -n ${n} -k 7 -p 1 -i 150 -o outputs/output_${n}_${k}.txt > logs/lik${n}${k}.log
    	echo "Terminei" ${n} ${k}
    done
done

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor 
