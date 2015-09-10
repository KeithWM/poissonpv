#!/bin/bash
batch="C64"
T=5
Tend=1
dto=.1
synch=1
order=fractal
#for dt in 1 .1 .01 .001 .0001 
for dt in .01
  do
#for compo in "S" "LT"
for compo in "S"
  do
#for C in 16 32 64 128
for C in {1..12}
#for C in 1
  do
    N=$[C*64]
    cmd="./sphere_C12 -compo S -synch $synch -N $N -order $order -batch ref -T 1 -Tend $Tend -dto $dto -dt .0001"
    echo $cmd
#    eval $cmd
#for synch in 0 1
for synch in 1
  do
#for order in "naive" "smart" "fractal"
for order in "fractal_new"
  do
    cmd="./sphere_C1 -compo $compo -synch $synch -N $N -order $order -batch $batch -T $T -Tend $Tend -dto $dto -dt $dt"
    echo $cmd
#    eval $cmd
    cmd="./sphere_C$C -compo $compo -synch $synch -N $N -order $order -batch $batch -T $T -Tend $Tend -dto $dto -dt $dt"
    echo $cmd
    eval $cmd
  done
  done
  done
  done
  done
