$MATLABROOT/bin/mex -client engine -v -output test test.c 
#$MATLABROOT/bin/mex -v -largeArrayDims -f ~/matopts.sh CFLAGS="-std=c99 -v" -output test test.c 
for C in {1..1}
#for C in 16 32 64 128
  do
    echo "Compiling for $C cores"
#    $MATLABROOT/bin/mex -f ~/matopts.sh CFLAGS="-std=c99" -DTHREADS=$C -Dchar16_t=uint16_t -o sphere_C$C sphere.c readMatlab.c writeMatlab.c normal.c parallel.c spherestep.c conserved.c store.c
    $MATLABROOT/bin/mex -client engine -v CFLAGS="-std=c99" -DTHREADS=$C -Dchar16_t=uint16_t -output sphere_C$C sphere.c readMatlab.c writeMatlab.c normal.c parallel.c spherestep.c conserved.c store.c
#    $MATLABROOT/bin/mex -v CFLAGS="-v" -DTHREADS=$C -Dchar16_t=uint16_t -output test test.c
  done
