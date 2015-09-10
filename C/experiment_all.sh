#!/bin/bash
batch="all"
T=10;
for N in 16 32 48 64 80 96 112 128 144 160 176 192 360 2520 27720
#for N in 16
  do
#for dt in .001
for dt in 1 .1 .01 .001 .0001 
  do
#for C in 1
for C in {1..12}
  do
for compo in "S" "LT"
  do
for synch in 0 1
  do
for order in "naive" "smart" "fractal"
  do
    cmd="./sphere_C$C -compo $compo -synch $synch -N $N -order $order -batch $batch -T $T -Tend 10 -dto 1 -dt $dt"
    echo $cmd
    eval $cmd
  done
  done
  done
  done
  done
  done
