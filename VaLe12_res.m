close all
clear all

ref = true;
refs = {'sphere','rk4'};
ref_names = {'Hopf','RK4'}
leg_names = {'Strang splitting', ref_names{:}};
leg_namesM6 = {'Strang splitting', ref_names{:}, 'McLachlan 6th order'};

styles = {'-','-.','--',':'};
markers = {'s','o','^','v','+','x'};
figures = 1:3; % stable vortex ring
% figures = [6:7, 16:17]; % collapse
% figures = [4:5, 14]; % Von Karman streets
% figures = 8:9; % unstable vortex ring
% figures = 10:15; % generic ics
% figures = 1:17;
for i=figures
    figure(i)
%     set(gca,'ColorOrder',[.5 0 0; 0 1 1]);
    set(gca,'ColorOrder',[0 0 0]);
    set(gca,'LineStyleOrder',styles);
    set(i,'PaperUnits','centimeters','PaperPosition',[0 0 10 6]);
    set(i,'Position',[400 400 500 300]);
    hold on
end
load goldenmap.mat
map = map(1:length(styles),:);

%% vortex rings
if nnz(figures == 1)
    inputdir = 'C/VaLe12/rings/Tend1.00e+03/C1_compoS_synch1_N6_orderfractal_dt1.00e-01';
    list = dir(sprintf('%s/expid*', inputdir));
    in = load(sprintf('%s/%s', inputdir, list(end).name));

    figure(1)
    set(gca,'YScale','log')
    line = plot(in.t, abs(in.H-in.H(1)));
    set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
    hold all

    figure(2)
    set(gca,'YScale','log')
    line = plot(in.t, sqrt(sum(bsxfun(@minus, in.J(1:3,:), in.J(1:3,1)).^2, 1)));
    set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
    hold all
    
    inputdir = 'C/VaLe12/rings/Tend1.00e+03/C1_compoM6_synch1_N6_ordersmart_dt1.00e-01';
    list = dir(sprintf('%s/expid*', inputdir));
    in = load(sprintf('%s/%s', inputdir, list(end).name));

    if ref
        for r = 1:length(refs)
            inref = load(sprintf('../VaLe/worksheets/polvani_dritschel/data/pd_%s_long.mat', refs{r}));

            figure(1)
            line = plot(inref.times, abs(inref.energies-inref.energies(1)));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
            hold all

            figure(2)
            line = plot(inref.times, sqrt(sum(bsxfun(@minus, inref.moments, inref.moments(1,:)).^2, 2)));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
            hold all
        end
    end

    figure(1)
    set(gca,'YScale','log')
    line = plot(in.t, abs(in.H-in.H(1)));
    set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));

    figure(2)
    set(gca,'YScale','log')
    line = plot(in.t, sqrt(sum(bsxfun(@minus, in.J(1:3,:), in.J(1:3,1)).^2, 1)));
    set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));

    figure(1)
    xlim([0 max(in.t)])
    xlabel('$t$','Interp','LaTeX');
    ylabel('Energy error','Interp','LaTeX');
    leg = legend(leg_namesM6);
    set(leg, 'Interp','LaTeX');
    saveas(1, 'figures/ringsH.eps', 'epsc');

    figure(2)
    xlim([0 max(in.t)])
    xlabel('$t$','Interp','LaTeX');
    ylabel('Momentum error','Interp','LaTeX');
    leg = legend(leg_namesM6);
    set(leg, 'Interp','LaTeX');
    saveas(2, 'figures/ringsJ.eps', 'epsc');

    ids = in.t < 100;

    figure(3)
    line = plot(in.t(ids), in.x(1,ids),'-','Color',map(1,:));
    hold all
    line = plot(in.t(ids), in.y(1,ids),'-','Color',map(2,:));
    line = plot(in.t(ids), in.z(1,ids),'-','Color',map(3,:));
    leg = legend('x','y','z');
    set(leg,'Interp','LaTeX');
    xlabel('$t$','Interp','LaTeX');
    ylabel('Components of the first vortex','Interp','LaTeX');
    saveas(3, 'figures/ringsx.eps', 'epsc');
end

%% spherical Von Karman vortex sheet
if nnz(figures == 4)
    inputdir = 'C/VaLe12/street/Tend1.00e+04/C1_compoS_synch1_N12_orderfractal_dt5.00e-01';
%     inputdir = 'C/VaLe12/sheet/Tend1.00e+04/C1_compoM6_synch1_N12_ordersmart_dt5.00e-01';
    list = dir(sprintf('%s/expid*', inputdir));
    in = load(sprintf('%s/%s', inputdir, list(end).name));

    figure(4)
    % line = plot(in.t, abs(in.H-in.H(1)))
    line = plot(in.t, in.H-in.H(1));
    set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
    hold all

    figure(5)
    set(gca,'YScale','log')
    line = plot(in.t, abs(in.J(4,:)-in.J(4,1)));
    set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
    hold all

    figure(14)
    line = plot(in.t, in.z(11,:));
    set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
    hold all

    if ref
        for r = 1:length(refs)
            inref = load(sprintf('../VaLe/worksheets/vortex_street/data/svs5_poles_%s_long.mat', refs{r}));

            figure(4)
            line = plot(inref.times, inref.energies-inref.energies(1));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));

            figure(5)
            line = plot(inref.times, sqrt(sum(bsxfun(@minus, inref.moments, inref.moments(1,:)).^2, 2)));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));

            inref = load(sprintf('../VaLe/worksheets/vortex_street/data/svs5_poles_%s.mat', refs{r}));
            figure(14)
            line = plot(inref.times, inref.vortices(:,11,3));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
        end
    end
    
%     figure(14)
%     hold on
%     line = plot(in.t, in.z(12,:));
%     set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
%     hold all
%     if ref
%         for r = 1:length(refs)
%             inref = load(sprintf('../VaLe/worksheets/vortex_street/data/svs5_poles_%s.mat', refs{r}));
%             figure(14)
%             line = plot(inref.times, inref.vortices(:,12,3));
%             set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
%         end
%     end

    figure(4)
    xlim([0 max(in.t)])
    xlabel('$t$','Interp','LaTeX');
    ylabel('Energy error','Interp','LaTeX');
    leg = legend(leg_names);
    set(leg, 'Interp','LaTeX');
    saveas(4, 'figures/sheetH.eps', 'epsc');

    figure(5)
    xlim([0 max(in.t)])
    xlabel('$t$','Interp','LaTeX');
    ylabel('Momentum error','Interp','LaTeX');
    leg = legend(leg_names);
    set(leg, 'Interp','LaTeX');
    saveas(5, 'figures/sheetJ.eps', 'epsc');

    figure(14)
%     xlim([0 50])
    leg = legend(leg_names);
    set(leg,'Interp','LaTeX','Location','West');
    xlabel('$t$','Interp','LaTeX');
    ylabel('$z_{11}$, $z_{12}$','Interp','LaTeX');
    saveas(14, 'figures/sheetz1112.eps', 'epsc');
end

%% collapse
if nnz(figures == 6)
    for pt = 1:4
        inputdir = sprintf('C/VaLe12/collapse/Tend5.00e+02/C1_compoS_synch1_N3_orderfractal_dt1.00e-0%d', pt);
        list = dir(sprintf('%s/expid*', inputdir));
        in = load(sprintf('%s/%s', inputdir, list(end).name));
        ids = in.t < 500;

        figure(6)
        line = plot(in.t(ids), abs(in.H(ids)-in.H(1)));
        set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
        hold all

        figure(7)
        set(gca,'YScale','log')
        line = plot(in.t(ids), sqrt(sum(bsxfun(@minus, in.J(1:3,ids), in.J(1:3,1)).^2)));
        set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
        hold all

        ids = in.t < 15;
        figure(16)
        line = plot(in.t(ids), in.H(ids)-in.H(1));
        set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
        hold all

        figure(17)
        set(gca,'YScale','log')
        line = plot(in.t(ids), sqrt(sum(bsxfun(@minus, in.J(1:3,ids), in.J(1:3,1)).^2)));
        set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
        hold all
    end

%     inputdir = sprintf('C/VaLe12/collapse/Tend5.00e+02/C1_compoM6_synch1_N3_ordersmart_dt1.00e-02');
%     list = dir(sprintf('%s/expid*', inputdir));
%     in = load(sprintf('%s/%s', inputdir, list(end).name));
%     ids = in.t < 15;
%     
%     figure(16)
%     line = plot(in.t(ids), in.H(ids)-in.H(1),'k');
% %     set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
%     hold all
%     
%     figure(17)
%     line = plot(in.t(ids), sqrt(sum(bsxfun(@minus, in.J(1:3,ids), in.J(1:3,1)).^2)),'k');
% %     set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
%     hold all

    if 0
        for r = 1:length(refs)
            inref = load(sprintf('../VaLe/worksheets/collapse3/data/collapse3_%s_sigma10_sim_long.mat', refs{r}));

            figure(6)
            line = plot(inref.times, abs(inref.energies-inref.energies(1)));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
            hold all

            figure(7)
            line = plot(inref.times, sqrt(sum(bsxfun(@minus, inref.moments, inref.moments(1,:)).^2, 2)));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
            hold all

            inref = load(sprintf('../VaLe/worksheets/collapse3/data/collapse3_%s_sigma10_sim.mat', refs{r}));
            figure(16)
            line = plot(inref.times, inref.energies-inref.energies(1));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
            hold all

            figure(17)
            line = plot(inref.times, sqrt(sum(bsxfun(@minus, inref.moments, inref.moments(1,:)).^2, 2)))
            hold all
        end
    end

    figure(16)
    T = 4*pi*(sqrt(23)-sqrt(17));
    line = plot([T T], [0 5e-3],'k:');
    hold on

    ids = in.t < 500;
    figure(6)
    set(gca, 'YScale', 'Log');
    xlim([0 max(in.t(ids))])
    xlabel('$t$','Interp','LaTeX');
    ylabel('Energy error','Interp','LaTeX');
    leg = legend('$\tau = 10^{-1}$','$\tau = 10^{-2}$','$\tau = 10^{-3}$','$\tau = 10^{-4}$');
    set(leg,'Interp','LaTeX','Location','North');
    saveas(6, 'figures/collapseH.eps', 'epsc');

    figure(7)
    xlim([0 max(in.t(ids))])
    xlabel('$t$','Interp','LaTeX');
    ylabel('Momentum error','Interp','LaTeX');
    leg = legend('$\tau = 10^{-1}$','$\tau = 10^{-2}$','$\tau = 10^{-3}$','$\tau = 10^{-4}$');
    set(leg,'Interp','LaTeX','Location','NorthWest');
    saveas(7, 'figures/collapseJ.eps', 'epsc');

    ids = in.t < 15;
    figure(16)
    xlim([0 max(in.t(ids))])
    xlabel('$t$','Interp','LaTeX');
    ylabel('Energy error','Interp','LaTeX');
    leg = legend('$\tau = 10^{-1}$','$\tau = 10^{-2}$','$\tau = 10^{-3}$','$\tau = 10^{-4}$');
    set(leg,'Interp','LaTeX','Location','NorthWest');
    saveas(16, 'figures/collapse_shortH.eps', 'epsc');

    figure(17)
    xlim([0 max(in.t(ids))])
    xlabel('$t$','Interp','LaTeX');
    ylabel('Momentum error','Interp','LaTeX');
    leg = legend('$\tau = 10^{-1}$','$\tau = 10^{-2}$','$\tau = 10^{-3}$','$\tau = 10^{-4}$');
    set(leg,'Interp','LaTeX','Location','NorthWest');
    saveas(17, 'figures/collapse_shortJ.eps', 'epsc');
end



%% vortex rings "Large" ensemble Kelvin-Helmholts inst sheet
if nnz(figures == 8)
%     inputdir = 'C/VaLe12/KH/Tend1.00e+03/C1_compoS_synch1_N40_orderfractal_dt1.00e-01';
    inputdir = 'C/VaLe12/KH/Tend1.00e+03/C1_compoM6_synch1_N40_ordersmart_dt1.00e-01';
    list = dir(sprintf('%s/expid*', inputdir));
    in = load(sprintf('%s/%s', inputdir, list(end).name));

    figure(8)
    line = plot(in.t, in.H-in.H(1));
    set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
    hold all

    figure(9)
    line = plot(in.t, abs(in.J(4,:)-in.J(4,1)));
    set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
    hold all

    if ref
        for r = 1:length(refs)
            inref = load(sprintf('../VaLe/worksheets/polvani_dritschel/data/equatorial_%s_long.mat', refs{r}));

            figure(8)
            line = plot(inref.times, inref.energies-inref.energies(1));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));

            figure(9)
            line = plot(inref.times, sqrt(sum(bsxfun(@minus, inref.moments, inref.moments(1,:)).^2, 2)));
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));

%             figure(18)
%             line = plot(inref.times, inref);
%             set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
        end
    end

    figure(8)
    xlim([0 300])
    leg = legend(leg_names);
    set(leg,'Interp','LaTeX','Location','NorthEast');
    xlabel('$t$','Interp','LaTeX');
    ylabel('Energy error','Interp','LaTeX');
    saveas(8, 'figures/KHH.eps', 'epsc');

    figure(9)
    set(gca,'YScale','log')
    xlim([0 300])
    leg = legend(leg_names);
    set(leg,'Interp','LaTeX','Location','East');
    xlabel('$t$','Interp','LaTeX');
    ylabel('Momentum error','Interp','LaTeX');
    saveas(9, 'figures/KHJ.eps', 'epsc');
end

%% many vortices
if nnz(figures == 10)
%     long = '_long';
%     Tend = 100;
%     dt = 1e-1;
%     refrs = [2];
    long = '';
    Tend = 1;
    dt = 1e-4;
    refrs = [1:3];

    refs = {'sphere','rk4','mp','lp'};
    ref_names = {'Hopf','RK4','Midpoint','Lie-Poisson'}
    leg_names = {'splitting', ref_names{refrs}};
    styles = {styles{1}, styles{refrs+1}}
    map = map([1 refrs+1],:)
    energies = {'lowest','neutral','highest'};


    for t=[1 3]
        energy = energies{t};
        inputdir = sprintf('C/many/%s/Tend%.2e/C3_compoS_synch1_N48_ordersmart_dt%.2e', energy, Tend, dt);
        list = dir(sprintf('%s/expid*', inputdir));
        in = load(sprintf('%s/%s', inputdir, list(end).name));

        figure(9+t);
        set(gca,'LineStyleOrder',styles);
        line = plot(in.t, abs(in.H-in.H(1)));
%         line = plot(in.t, in.H);
        set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
        hold all
        for r = refrs
            inref = load(sprintf('../VaLe/worksheets/many/data/%s_M8N48%s_%s.mat', energy, long, refs{r}));
            line = plot(inref.times, abs(inref.energies-inref.energies(1)));
%             line = plot(inref.times, inref.energies);
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
            hold all
        end
        set(gca, 'YScale', 'log');
        leg = legend(leg_names);
        set(leg,'Interp','LaTeX','Location','SouthEast');
        xlim([0 max(inref.times)])
%         ylim(10.^[-2 1]);
        xlabel('$t$','Interp','LaTeX');
        ylabel('Energy error','Interp','LaTeX');
        saveas(9+t, sprintf('figures/%s_M8N48%s.eps', energy, long), 'epsc');

        figure(12+t);
        set(gca,'LineStyleOrder',styles);
        line = plot(in.t, sqrt(sum(bsxfun(@minus, in.J(1:3,:), in.J(1:3,1)).^2)));
%         line = plot(in.t, in.H);
        set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
        hold all
        for r = refrs
            inref = load(sprintf('../VaLe/worksheets/many/data/%s_M8N48%s_%s.mat', energy, long, refs{r}));
            line = plot(inref.times, sqrt(sum(bsxfun(@minus, inref.moments, inref.moments(1,:)).^2, 2)))
%             line = plot(inref.times, inref.energies);
            set(line, 'Color', map(strcmp(styles, get(line, 'LineStyle')),:));
            hold all
        end
        set(gca, 'YScale', 'log');
        leg = legend(leg_names);
        set(leg,'Interp','LaTeX','Location','NorthWest');
        xlim([0 max(inref.times)])
%         ylim(10.^[-2 1]);
        xlabel('$t$','Interp','LaTeX');
        ylabel('Momentum error','Interp','LaTeX');
%         saveas(12+t, sprintf('figures/momentum_%s_M8N48%s.eps', energy, long), 'epsc');
    end
end