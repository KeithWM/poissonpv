command = './sphere_C1 -energy street -batch convergence -Tend 100. -compo %s -N 12 -order fractal -dt %f -dto 1';

compos = {'LT','S','M4','M6','Y4','Y6'};
dts = kron([1 2 5], 10.^[-1 -2 -3 -4]);

for compo = compos
    for dt = dts
        unix(sprintf(command, compo{1}, dt));
    end
end