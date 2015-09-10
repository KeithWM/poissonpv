command = './sphere_C1 -energy collapse -batch convergence -Tend 50. -compo %s -N 3 -order fractal -dt %f -dto .5';

compos = {'LT','S','M4','M6','Y4','Y6'};
dts = kron([1 2 5], 10.^[-1 -2 -3 -4]);
dts = .2

for compo = compos
    for dt = dts
        unix(sprintf(command, compo{1}, dt));
    end
end