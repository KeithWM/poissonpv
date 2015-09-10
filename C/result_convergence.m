clear all
close all
styles = {'-','-.','--',':'};
markers = {'s','o','^','v','+','x'};
names = {'Stable vortex ring', 'Von K\''arm\''an vortex streets', 'Collapsing vortices'};

for i=1:6
    fig = figure(i);
    set(fig,'PaperUnits','centimeters','PaperPosition',[0 0 15 9]);
    set(fig,'Position',[400 400 500 300]);
%     set(findall(fig,'-property','FontSize'),'FontSize',24)
    hold on
end
load ../goldenmap.mat
%%

base = 'convergence/rings/Tend1.00e+03/C1_compo%s_synch1_N6_orderfractal_dt%.2e';

compos = {'LT','S','Y4','M4','Y6','M6'};
dts = kron(10.^[-1 -2 -3 -4], [5 2 1]);

errors_H = zeros(length(dts), length(compos));
errors_r = zeros(length(dts), length(compos));

path = sprintf(base, 'M6', .0001);
list = dir(path);
in_ref = load(sprintf('%s/%s', path, list(end).name));
H_ref = in_ref.H(1);
x_ref = in_ref.x;
y_ref = in_ref.y;
z_ref = in_ref.z;

L = 10;

for i=1:length(dts)
    dt = dts(i);
    for j=1:length(compos)
        compo = compos{j};
        path = sprintf(base, compo, dt);
        list = dir(path);
        in = load(sprintf('%s/%s', path, list(end).name));
        
        errors_H(i,j) = norm(in.H(1:L) - H_ref, 2);
        errors_r(i,j) = sqrt(sum(sum((in.x-x_ref).^2 + (in.y-y_ref).^2 + (in.z-z_ref).^2)));
    end
end

base = '../../VaLe/worksheets/polvani_dritschel/data/pd_%s_long_h%+.2e.mat';
methods = {'rk2','rk4'};

in_ref = load(sprintf(base, 'rk4', .0001));
H_ref = in_ref.energies(1);
r_ref = in_ref.vortices;

errors_H_VaLe = zeros(length(dts), length(methods))*NaN;
errors_r_VaLe = zeros(length(dts), length(methods))*NaN;

for i=1:length(dts(1:10))
    dt = dts(i);
    for j=1:length(methods)
        method = methods{j};
        in = load(sprintf(base, method, dt));
        
        errors_H_VaLe(i,j) = norm(in.energies(1:L) - H_ref, 2);
        errors_r_VaLe(i,j) = sqrt(sum(sum(sum((in.vortices - r_ref).^2))));
    end
end

for j=1:size(errors_H,2)
    [~,best] = min(errors_H(:,j));
    errors_H(best+1:end, j) = NaN;
end

for j=1:size(errors_H_VaLe,2)
    [~,best] = min(errors_H_VaLe(:,j));
    errors_H_VaLe(best+1:end, j) = NaN;
end

workloads = bsxfun(@times, [1 2 [3 5 7 9]*2]*1000, 1./dts');
workloads_VaLe = bsxfun(@times, [2 4]*1000, 1./dts');

figure(1)
for j=1:length(compos)
    loglog(workloads(:,j),errors_H(:,j),'-','Color',map(j,:), 'Marker', markers{j})
    hold on
end
for j=1:length(methods)
    loglog(workloads_VaLe(:,j),errors_H_VaLe(:,j),'-.','Color',map(j,:), 'Marker', markers{j})
    hold on
end

figure(2)
for j=1:length(compos)
    loglog(workloads(:,j),errors_r(:,j),'-','Color',map(j,:), 'Marker', markers{j})
    hold on
end
for j=1:length(methods)
    loglog(workloads_VaLe(:,j),errors_r_VaLe(:,j),'-.','Color',map(j,:), 'Marker', markers{j})
    hold on
end

%%

base = 'convergence/street/Tend1.00e+02/C1_compo%s_synch1_N12_orderfractal_dt%.2e';

compos = {'LT','S','Y4','M4','Y6','M6'};
dts = kron(10.^[-1 -2 -3 -4], [5 2 1]);

errors_H = zeros(length(dts), length(compos));
errors_r = zeros(length(dts), length(compos));

path = sprintf(base, 'M6', .0001);
list = dir(path);
in_ref = load(sprintf('%s/%s', path, list(end).name));
H_ref = in_ref.H(1);
x_ref = in_ref.x;
y_ref = in_ref.y;
z_ref = in_ref.z;

L = 10;

for i=1:length(dts)
    dt = dts(i);
    for j=1:length(compos)
        compo = compos{j};
        path = sprintf(base, compo, dt);
        list = dir(path);
        in = load(sprintf('%s/%s', path, list(end).name));
        
        errors_H(i,j) = norm(in.H(1:L) - H_ref, 2);
        errors_r(i,j) = sqrt(sum(sum((in.x-x_ref).^2 + (in.y-y_ref).^2 + (in.z-z_ref).^2)));
    end
end

base = '../../VaLe/worksheets/vortex_street/data/svs5_%s_h%+.2e.mat';
methods = {'rk2','rk4'};

in_ref = load(sprintf(base, 'rk4', .0001));
H_ref = in_ref.energies(1);
r_ref = in_ref.vortices;

errors_H_VaLe = zeros(length(dts), length(methods))*NaN;
errors_r_VaLe = zeros(length(dts), length(methods))*NaN;

for i=1:length(dts(1:10))
    dt = dts(i);
    for j=1:length(methods)
        method = methods{j};
        in = load(sprintf(base, method, dt));
        
        errors_H_VaLe(i,j) = norm(in.energies(1:L) - H_ref, 2);
        errors_r_VaLe(i,j) = sqrt(sum(sum(sum((in.vortices - r_ref).^2))));
    end
end

for j=1:size(errors_H,2)
    [~,best] = min(errors_H(:,j));
    errors_H(best+1:end, j) = NaN;
end

for j=1:size(errors_H_VaLe,2)
    [~,best] = min(errors_H_VaLe(:,j));
    errors_H_VaLe(best+1:end, j) = NaN;
end

figure(3)
for j=1:length(compos)
    loglog(workloads(:,j),errors_H(:,j),'-','Color',map(j,:), 'Marker', markers{j})
    hold on
end
for j=1:length(methods)
    loglog(workloads_VaLe(:,j),errors_H_VaLe(:,j),'-.','Color',map(j,:), 'Marker', markers{j})
    hold on
end

figure(4)
for j=1:length(compos)
    loglog(workloads(:,j),errors_r(:,j),'-','Color',map(j,:), 'Marker', markers{j})
    hold on
end
for j=1:length(methods)
    loglog(workloads_VaLe(:,j),errors_r_VaLe(:,j),'-.','Color',map(j,:), 'Marker', markers{j})
    hold on
end

%%

base = 'convergence/collapse/Tend5.00e+01/C1_compo%s_synch1_N3_orderfractal_dt%.2e';

compos = {'LT','S','Y4','M4','Y6','M6'};
dts = kron(10.^[-1 -2 -3 -4], [5 2 1]);

errors_H = zeros(length(dts), length(compos))*NaN;
errors_r = zeros(length(dts), length(compos))*NaN;

path = sprintf(base, 'M6', .0001);
list = dir(path);
in_ref = load(sprintf('%s/%s', path, list(end).name));
H_ref = in_ref.H(1);
x_ref = in_ref.x;
y_ref = in_ref.y;
z_ref = in_ref.z;

LH = 10;
L = 10;

for i=4:length(dts)
    dt = dts(i);
    for j=1:length(compos)
        compo = compos{j};
        path = sprintf(base, compo, dt);
        list = dir(path);
        in = load(sprintf('%s/%s', path, list(end).name));
        
        errors_H(i,j) = norm(in.H(1:LH) - H_ref, 2);
        errors_r(i,j) = sqrt(sum(sum((in.x(:,1:L)-x_ref(:,1:L)).^2 + (in.y(:,1:L)-y_ref(:,1:L)).^2 + (in.z(:,1:L)-z_ref(:,1:L)).^2)));
    end
    figure(7)
    plot3(in.x', in.y', in.z', '-', 'Color', map(i,:))
    hold on
end

base = '../../VaLe/worksheets/collapse3/data/collapse3_%s_h%+.2e.mat';
methods = {'rk2','rk4'};

in_ref = load(sprintf(base, 'rk4', .0001));
H_ref = in_ref.energies(1);
r_ref = in_ref.vortices;

errors_H_VaLe = zeros(length(dts), length(methods))*NaN;
errors_r_VaLe = zeros(length(dts), length(methods))*NaN;

for i=4:length(dts(1:10))
    dt = dts(i);
    for j=1:length(methods)
        method = methods{j};
        in = load(sprintf(base, method, dt));
        
        errors_H_VaLe(i,j) = norm(in.energies(1:LH) - H_ref, 2);
        errors_r_VaLe(i,j) = sqrt(sum(sum(sum((in.vortices(1:L,:,:) - r_ref(1:L,:,:)).^2))));
    end
    figure(7)
    plot3(in.vortices(:,:,1), in.vortices(:,:,2), in.vortices(:,:,3), '-.', 'Color', map(i,:))
    hold on
end

for j=1:size(errors_H,2)
    [~,best] = min(errors_H(:,j));
    errors_H(best+1:end, j) = NaN;
end

for j=1:size(errors_H_VaLe,2)
    [~,best] = min(errors_H_VaLe(:,j));
    errors_H_VaLe(best+1:end, j) = NaN;
end

figure(5)
for j=1:length(compos)
    loglog(workloads(:,j),errors_H(:,j),'-','Color',map(j,:), 'Marker', markers{j})
    hold on
end
for j=1:length(methods)
    loglog(workloads_VaLe(:,j),errors_H_VaLe(:,j),'-.','Color',map(j,:), 'Marker', markers{j})
    hold on
end

figure(6)
for j=1:length(compos)
    loglog(workloads(:,j),errors_r(:,j),'-','Color',map(j,:), 'Marker', markers{j})
    hold on
end
for j=1:length(methods)
    loglog(workloads_VaLe(:,j),errors_r_VaLe(:,j),'-.','Color',map(j,:), 'Marker', markers{j})
    hold on
end

%%
for i=1:6
    fig=figure(i);
    if i > 4
        mino = 3;
    else
        mino = 3;
    end
    xlim(10.^[mino, 7])
    xlabel('workload','Interpreter','LaTeX')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'Xtick', 10.^(mino:7))
    set(findall(fig,'-property','FontSize'),'FontSize',14)
    if i > 4
        legend({compos{:}, 'RK2', 'RK4'},'Interpreter','LaTeX', 'Location', 'West','FontSize',18)
    end
end

energies = {'rings','street','collapse'}
for i=1:3
    i1 = 2*i-1;
    i2 = 2*i;
    energy = energies{i};
    
    figure(i1)
    title(names{i}, 'Interpreter', 'LaTeX')
    ylabel('energy error','Interpreter','LaTeX')
    saveas(i1,sprintf('../figures/errors_H_%s.svg', energy),'svg')
    
    figure(i2)
    title(names{i}, 'Interpreter', 'LaTeX')
    ylabel('position error','Interpreter','LaTeX')
    saveas(i2,sprintf('../figures/errors_r_%s.svg', energy),'svg')
%     legend({compos{:}, 'RK2', 'RK4'},'Interpreter','LaTeX', 'Location', 'EastOutside')
end