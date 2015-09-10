#!/bin/bash
batch="C_compo_N360"
T=10;
#for N in 16 32 48 64 80 96 112 128 144 160 176 192 360 2520 27720
for N in 360
  do
#for dt in 1 .1 .01 .001 .0001 
for dt in .001
  do
#for synch in 0 1
for synch in 1
  do
#for order in "naive" "smart" "fractal"
for order in "fractal"
  do
    cmd="./sphere_C12 -compo S -synch $synch -N $N -order $order -batch ref -T 1 -Tend .1 -dto .01 -dt .0001"
    echo $cmd
    eval $cmd
for C in {1..12}
#for C in 1
  do
for compo in "S" "LT"
#for compo in "S"
  do
    cmd="./sphere_C$C -compo $compo -synch $synch -N $N -order $order -batch $batch -T $T -Tend .1 -dto .01 -dt $dt"
    echo $cmd
    eval $cmd
  done
  done
  done
  done
  done
  done
