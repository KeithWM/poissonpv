#!/bin/bash
batch="BEA"
T=1
Tend=1
synch=1
order=fractal
#for N in 16 32 48 64 80 96 112 128 144 160 176 192 360 2520 27720
for N in 4
  do
#for dt in .1 .01 .001 .0001 
for dt in 1
  do
#for compo in "S" "LT"
for compo in "S"
  do
#for C in {4..12}
for C in 1
  do
#for synch in 0 1
for synch in 1
  do
#for order in "naive" "smart" "fractal"
for order in "fractal"
  do
    cmd="./sphere_C$C -compo $compo -synch $synch -N $N -order $order -batch $batch -T $T -Tend $Tend -dto $dt -dt $dt"
    echo $cmd
    eval $cmd
  done
  done
  done
  done
  done
  done
