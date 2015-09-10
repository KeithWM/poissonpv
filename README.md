# poissonpv

By Keith Myerscough CC BY-NC-SA.

Licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International licence,
for more information, visit http://creativecommons.org/licenses/by-nc-sa/4.0/

The code developed here has been developed to be used on a single machine and later on a different machine, both Macintosh computers
the parallelization is specific to Mac's Grand Central Dispatch, and will not work on other types of computers. Nevertheless, I believe
it is good scientific practice to be open about one's code. That said, I have not made particularly much effort to make the code readable.
If you are interested in the code, but frustrated by some part that is not clear, BY ALL MEANS LET ME KNOW. I will be very happy to help
anyone who wants to study this code, you can reach me at
  keith@myerscough.nl via e-mail
  @KeithWM on twitter
  KeithWM on Github
  and other places under the same moniker

files need to be compiled from within Matlab using a command such as

cd ./C
!./compile_all_ceres.sh

and then executed as 
!./sphere_C$C -compo $compo -synch $synch -N $N -order $order -batch $batch -T $T -Tend $Tend -dto $dt -dt $dt

where the variables denoted with the cygil '$' need to be replaced for (valid) parameters, as detailed in sphere.c

