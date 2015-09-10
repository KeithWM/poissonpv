#!/bin/bash
batch="rings"
T=5
Tend=100
dto=1
synch=1
order=fractal
#for N in 16 32 48 64 80 96 112 128 144 160 176 192 360 2520 27720
for N in 6
  do
cmd="./sphere_C1 -compo S -synch $synch -N $N -order $order -batch ref -T 1 -Tend $Tend -dto $dto -dt .00001 -energy rings"
echo $cmd
eval $cmd
for dt in .0001 .0002 .0005 .001 .002 .005 .01 .02 .05 .1
#for dt in .001
  do
for compo in "S" "LT"
#for compo in "S"
  do
#for C in {4..12}
#for C in 12
for C in 1
  do
#for synch in 0 1
for synch in 1
  do
#for order in "naive" "smart" "fractal"
for order in "fractal"
  do
    cmd="./sphere_C$C -compo $compo -synch $synch -N $N -order $order -batch $batch -T $T -Tend $Tend -dto $dto -dt $dt -energy rings"
    echo $cmd
    eval $cmd
  done
  done
  done
  done
  done
  done
