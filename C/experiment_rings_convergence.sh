#!/bin/bash
batch="convergence"
T=5
Tend=1000
dto=1
synch=1
order=fractal
for N in 6
  do
cmd="./sphere_C1 -compo M6 -synch $synch -N $N -order $order -batch ref -T 1 -Tend $Tend -dto $dto -dt .00001 -energy rings"
echo $cmd
eval $cmd
for compo in "S" "LT" "M4" "M6" "Y4" "Y6"
  do
for dt in .0001 .0002 .0005 .001 .002 .005 .01 .02 .05 .1
  do
for C in 1
  do
for synch in 1
  do
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
