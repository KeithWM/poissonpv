clear all
close all

Ns = 50;
[xs,ys,zs] = sphere(Ns);
scale = .99
% xs = scale*xs;
% ys = scale*ys;
% zs = scale*zs;

load goldenmap.mat
icscale = .1

epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';

%% KH (large ring)
Tend = 1000
dt = .1
path = sprintf('C/VaLe_figs/KH/Tend%.2e/C1_compoS_synch1_N40_orderfractal_dt%.2e', Tend, dt);
list = dir(path);
path = [path, '/', list(end).name];

in = load(path);

nz = 201
figure(11)
set(11,'PaperUnits','centimeters','PaperPosition',[0 0 10 6]);
set(11,'Position',[400 400 500 300]);
% set(gca, 'ColorOrder', map)
plot(in.t(1:nz),in.z(:,1:nz)')
xlim([0, in.t(nz)])
xlabel('$t$','Interpreter','LaTeX')
ylabel('$z$','Interpreter','LaTeX')
% hgexport(11,'figures/KHz.eps',epsfig,'Format','eps')
saveas(11,'figures/KHz.svg','svg')

plotics(in, 'KH', .03);

%% collapse
Tend = 200
dt = .001
path = sprintf('C/VaLe_figs/collapse/Tend%.2e/C1_compoS_synch1_N3_orderfractal_dt%.2e', Tend, dt);
list = dir(path);
path = [path, '/', list(end).name];

in = load(path);
in_null = in;
in_null.G = zeros(0,1);
in_null.x = zeros(0,1);
in_null.y = zeros(0,1);
in_null.z = zeros(0,1);

handle = plotics(in_null, 'temp', icscale);
n3d = 8454;
figure(handle)
% hold off
% set(12,'PaperUnits','centimeters','PaperPosition',[0 0 8 8]);
% set(12,'Position',[400 400 400 400]);
% for r=1:-.01:.95
%     surf(r*xs,r*ys,r*zs,'FaceColor',[1 1 1]*.85,'FaceAlpha',.2,'EdgeColor','none') % sphere
%     hold all
% end
% set(gca, 'ColorOrder', map)
% plot3(in.x(:,1:n3d)', in.y(:,1:n3d)', in.z(:,1:n3d)')
linemap = [[.3,0,0];[.9,1,.9];[0,0,1]];
for n=1:size(in.x,1)
    plot3(in.x(n,1:n3d)', in.y(n,1:n3d)', in.z(n,1:n3d)', 'Color', linemap(n,:))
end
view([-60 60])
axis('equal')
axis('off')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hgexport(12,'figures/collapse3d.eps',epsfig,'Format','eps')
% set(12, 'InvertHardcopy', 'off')
saveas(handle,'figures/collapse3d.svg','svg')

handle = plotics(in_null, 'temp', icscale);
n3d2 = 64499;
figure(handle)
% hold off
% set(112,'PaperUnits','centimeters','PaperPosition',[0 0 8 8]);
% set(112,'Position',[400 400 400 400]);
% for r=1:-.01:.95
%     surf(r*xs,r*ys,r*zs,'FaceColor',[1 1 1]*.85,'FaceAlpha',.2,'EdgeColor','none') % sphere
%     hold all
% end
% set(gca, 'ColorOrder', map)
% plot3(in.x(:,n3d:n3d2)', in.y(:,n3d:n3d2)', in.z(:,n3d:n3d2)')
for n=1:size(in.x,1)
    plot3(in.x(n,n3d:n3d2)', in.y(n,n3d:n3d2)', in.z(n,n3d:n3d2)', 'Color', linemap(n,:))
end
view([-60 60])
axis('equal')
axis('off')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hgexport(112,'figures/collapse3d_post.eps',epsfig,'Format','eps')
% set(112, 'InvertHardcopy', 'off')
saveas(handle,'figures/collapse3d_post.svg','svg')

plotics(in, 'collapse', icscale);

%% ring
Tend = 1000
dt = .05
path = sprintf('C/VaLe_figs/rings/Tend%.2e/C1_compoS_synch1_N6_orderfractal_dt%.2e', Tend, dt);
list = dir(path);
path = [path, '/', list(end).name];

in = load(path);

nz = 10001
figure(13)
set(13,'PaperUnits','centimeters','PaperPosition',[0 0 10 6]);
set(13,'Position',[400 400 500 300]);
% set(gca, 'ColorOrder', map)
plot(in.t(1:nz),in.z(:,1:nz)'-.92)
xlim([0, in.t(nz)])
xlabel('$t$','Interpreter','LaTeX')
ylabel('$z - 0.92$','Interpreter','LaTeX')
% hgexport(13,'figures/ringz.eps',epsfig,'Format','eps')
saveas(13,'figures/ringz.svg','svg')

plotics(in, 'ring', icscale);

%% street
Tend = 10000
dt = .5
path = sprintf('C/VaLe_figs/street/Tend%.2e/C1_compoS_synch1_N12_orderfractal_dt%.2e', Tend, dt);
list = dir(path);
path = [path, '/', list(end).name];

in = load(path);

nz = 101
figure(14)
set(14,'PaperUnits','centimeters','PaperPosition',[0 0 10 6]);
% set(14,'Position',[400 400 500 300]);
set(gca, 'ColorOrder', map)
plot(in.t(1:nz),in.z(:,1:nz)')
xlim([0, in.t(nz)])
xlabel('$t$','Interpreter','LaTeX')
ylabel('$z$','Interpreter','LaTeX')
% hgexport(14,'figures/streetz.eps',epsfig,'Format','eps')
saveas(14,'figures/streetz.svg','svg')

plotics(in, 'street', icscale);

%% generic
energies = {'lowest', 'highest'};

for energy = energies
    in = load(sprintf('ic/ics_%s_M8N48.mat', energy{1}))
    mod.G = in.Gamma;
    mod.x = in.X0(1,:)';
    mod.y = in.X0(2,:)';
    mod.z = in.X0(3,:)';
    plotics(mod, sprintf('%s_M8N48', energy{1}), .1)
end

