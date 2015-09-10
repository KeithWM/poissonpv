#!/bin/bash
batch="many"
T=1
Tend=1
dto=.1
synch=1
order=fractal
#for N in 16 32 48 64 80 96 112 128 144 160 176 192 360 2520 27720
#for N in 24 48 120 240 480 1200 2400 4800 12000
for N in 48 384 
#for N in 2048
  do
#for dt in 1 .1 .01 .001 .0001 
for dt in .0001
  do
#for compo in "S" "LT"
for compo in "S"
  do
    cmd="./sphere_C12 -compo S -synch $synch -N $N -order $order -batch ref -T 1 -Tend $Tend -dto $dto -dt .01"
    echo $cmd
#    eval $cmd
#for C in {4..12}
for C in 12
  do
#for synch in 0 1
for synch in 1
  do
#for order in "naive" "smart" "fractal"
#for order in "fractal_new"
#for order in "naive" "smart" "fractal" "fractal_new"
for order in "smart"
  do
    cmd="./sphere_C$C -compo $compo -synch $synch -N $N -order $order -batch $batch -T $T -Tend $Tend -dto $dto -dt $dt"
    echo $cmd
    eval $cmd
  done
  done
  done
  done
  done
  done
