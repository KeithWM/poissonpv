#for C in {1..12}
for C in 16 32 64 128
  do
    N=$[C*64];
    for energy in "lowest" "neutral" "highest"
      do
        echo $N $energy;
#        ./ics -N $N -energy $energy
      done
  done

for N in 360 2520 27720
  do
    for energy in "lowest" "neutral" "highest"
      do
        echo $N $energy;
#        ./ics -N $N -energy $energy
      done
  done

#for N in 24 48 120 240 480 1200 2400 4800 12000
for N in 36 360 3600 
  do
    for energy in "neutral"
      do
        echo $N $energy;
        ./ics -N $N -energy $energy
      done
  done
