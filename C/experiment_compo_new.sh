#!/bin/bash
batch="compo_C8"
T=5
Tend=.1
dto=.01
synch=1
order=fractal
#for N in 16 32 48 64 80 96 112 128 144 160 176 192 360 2520 27720
for N in 360 
#for N in 2520
  do
cmd="./sphere_C12 -compo S -synch $synch -N $N -order $order -batch ref -T 1 -Tend $Tend -dto $dto -dt .00001"
echo $cmd
#eval $cmd
#for dt in .001 .002 .005 .01 
#for dt in .00005 .00002 .00001
for dt in .00001 .00002 .00005 .0001 .0002 .0005 .001 .002 .005 .01
#for dt in .02 .05 .1 
  do
#for compo in "Y4" "M4" "Y6" "M6"
for compo in "LT" "S"
#for compo in "S"
  do
#for C in {4..12}
for C in 8
  do
#for synch in 0 1
for synch in 1
  do
for order in "naive" "smart" "fractal" "fractal_new"
#for order in "fractal"
#for order in "smart"
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
