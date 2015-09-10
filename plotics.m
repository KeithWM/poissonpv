function handle = plotics(in, name, icscale)
    epsfig = hgexport('factorystyle');
    epsfig.Format = 'eps';
    disp(in.G)
    
    Ns = 50;
    [xs,ys,zs] = sphere(Ns);
    xsic = icscale*xs;
    ysic = icscale*ys;
    zsic = icscale*zs;
    csic = 0*zs; % a zero-sphere for use in colouring the balls
    scale = .99;
%     xs = scale*xs;
%     ys = scale*ys;
%     zs = scale*zs;
    
    cscale = .99
    tc = (0:2*Ns)/Ns*pi;
    xc = cscale*cos(tc);
    yc = cscale*sin(tc);
    zc = cscale*ones(size(tc));
    
    handle = figure('Name',sprintf('ic_%s', name));
    set(handle,'PaperUnits','centimeters','PaperPosition',[0 0 8 8]);
    set(handle,'Position',[400 400 400 400]);
    colormap('hot')
%     colormap('parula')
%     set(gca, 'ColorOrder', map)
    for r=1:-.01:.95
        surf(r*xs,r*ys,r*zs,'FaceColor',[1 1 1]*.85,'FaceAlpha',.2,'EdgeColor','none') % sphere
        hold all
    end
%     surf(xc,yc,zc,'FaceColor','k','FaceAlpha',.5,'EdgeColor','none') % circle
    for n=1:size(in.x,1)
        surf(in.x(n,1)+xsic,in.y(n,1)+ysic,in.z(n,1)+zsic,in.G(n)+csic,'FaceAlpha',1,'EdgeColor','none')
    end
    for lat = -60:30:60
        z = sin(lat/180*pi);
        r = sqrt(1-z^2);
        plot3(r*xc, r*yc, z*zc, '-','Color',[1,1,1]*.75)
    end
    for lon = -180:45:179
        phi = lon/180*pi;
        plot3(xc*cos(phi), xc*sin(phi), yc, '-','Color',[1,1,1]*.75)
    end
    view([-60, 60])
    if strcmp(name, 'street')
        view([-60 30])
    end
    set(gca, 'clim', [-1, 1]*1.25)
    axis('equal')
    axis('off')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     hgexport(handle,sprintf('figures/ic_%s.eps', name),epsfig,'Format','eps')
%     saveas(handle,sprintf('figures/ic_%s_saveas.eps', name),'epsc')
    saveas(handle,sprintf('figures/ic_%s.svg', name),'svg')
end