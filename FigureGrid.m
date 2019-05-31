% figure showing ways the grid can be conceptualized
% eric zilli - 20111106 - v1.0

%% Make some colormaps
% red to white to blue
RWB = ones(64,3); RWB(1:32,1) = linspace(0,1,32); RWB(1:32,2) = RWB(1:32,1); RWB(33:64,2) = flipud(RWB(1:32,1)); RWB(33:64,3) = RWB(33:64,2);
startColor = [0.07 0.23 0.28]/3;
midColor = [1 1 1];
endColor = [0.55 0.48 0.05]/1;

% less bright red to white to less bright blue
startColor = [0 0 0.9];
midColor = [1 1 1];
endColor = [0.9 0 0];
mRWB = [linspace(startColor(1),midColor(1),32) linspace(midColor(1),endColor(1),32); linspace(startColor(2),midColor(2),32) linspace(midColor(2),endColor(2),32); linspace(startColor(3),midColor(3),32) linspace(midColor(3),endColor(3),32)]';

% gray to white to gold
startColor = [0.07 0.23 0.28]/3;
midColor = [1 1 1];
endColor = [0.55 0.48 0.05]/1;
GWG = [linspace(startColor(1),midColor(1),32) linspace(midColor(1),endColor(1),32); linspace(startColor(2),midColor(2),32) linspace(midColor(2),endColor(2),32); linspace(startColor(3),midColor(3),32) linspace(midColor(3),endColor(3),32)]';

useColorMap = [flipud(GWG); RWB; mRWB];

%% Load and generate data to plot
% loads variable M
load data/BlairEtAl2007_Readout.mat

% loads variable thrt2
load data/Spatial_interference.mat

% Plot a matrix scaled and shifted so that its zero point corresponds
% to the halfway point in a 64-item colormap.
imagepn = @(v)image(64*(v./2/abs(eps+max(abs([max(max(v)) min(min(v))])))+0.5));
imagepnc = @(v)(64*(v./2/abs(eps+max(abs([max(max(v)) min(min(v))])))+0.5));

% make mesh to evaluate grid on
x = linspace(0,2,625);
y = linspace(0,2,625);
[X Y] = meshgrid(x,y);
Xp = reshape(X,1,[]);
Yp = reshape(Y,1,[]);

% directions of grid basis vectors
theta = [0 pi/3 2*pi/3]+pi/2;

% Project the mesh onto the preferred directions
H = [cos(theta') sin(theta')];
projMesh = H*[Xp; Yp];

% Evaluate the three gratings and multiply them into the grid
grating1 = .5+.5*cos(projMesh(1,:)*2*pi*5);
grating2 = .5+.5*cos(projMesh(2,:)*2*pi*5);
grating3 = .5+.5*cos(projMesh(3,:)*2*pi*5);
gridPattern = reshape(grating1.*grating2.*grating3,length(x),[]);

%% Figure set up
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')

% size on paper:
widthOnPaper = 18; % cm
heightOnPaper = 7; % cm

figure('units','centimeters','position',[1 1 widthOnPaper heightOnPaper],'color','w');
set(gcf, 'renderer', 'painter')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [widthOnPaper heightOnPaper]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 widthOnPaper heightOnPaper]);
colormap(useColorMap)

leftMargin = 0.03;
bottomMargin = 0.02;
nCols = 3;
nRows = 2;
lefts = leftMargin + 0.95*(0:nCols-1)/nCols;
bottoms = bottomMargin + [0 1/2];
widths = 0.9*[1/2 2.25/4 1/4];
heights = 0.9*[1.3 2.25/2];

shapeLineWidth = 2;

% the 2D map goes from 0,0 to 1,1 so we need some scaling factors to let us
% identify field center
yspacing = 1*0.1; % i.e. 10 fields along the y directionn
xspacing = yspacing*2/sqrt(3);


% 2D Plot
%% 2D arena with grid and spatial phases
axes('position',[lefts(1) bottoms(1)-0.16 widths(1) heights(1)])

% Due to a limitation in the eps file format we have to break this into
% smaller chunks.
nchunks = 5;
for rchunk=1:nchunks
  for cchunk=1:nchunks
    % subscripts into the grid pattern:
    rsubs = ((rchunk-1)*size(gridPattern,1)/nchunks+1):(rchunk*size(gridPattern,1)/nchunks);
    csubs = ((cchunk-1)*size(gridPattern,2)/nchunks+1):(cchunk*size(gridPattern,2)/nchunks);
    % locations where we'll draw that bit:
    xd = 1*[(cchunk-1)/nchunks cchunk/nchunks];
    yd = 1*[(nchunks-rchunk)/nchunks (nchunks-rchunk+1)/nchunks];
    surface('xdata',xd,'ydata',yd,'zdata',-1e-9*ones(2),'cdata',128+round(imagepnc(flipud(gridPattern(rsubs,csubs)))),'facecolor','texturemap','edgecolor','none','cdatamapping','direct');
  end
end

% Draw the "walls" of the spatial environment
surface('xdata',[0 0],'ydata',[0 1],'zdata',[0 4; 0 4],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
surface('xdata',[1 1],'ydata',[0 1],'zdata',[0 4; 0 4],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
surface('xdata',[0 1],'ydata',[0 0],'zdata',[4 4; 0 0],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
surface('xdata',[0 1],'ydata',[1 1],'zdata',[0 0; 4 4],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')

hold on;

set(gca,'cameraposition',[0 -1 36])
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
axis off

%% triangle
xs = xspacing/2+[xspacing 2*xspacing 1.5*xspacing xspacing 2*xspacing];
ys = [9 9 8 9 9]*yspacing;
p = patch(xs,ys,'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')

text(0.13*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     0.93*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(A)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')


%% hexagon
left = 3*xspacing;
xs = left+[0 xspacing/2 3*xspacing/2 2*xspacing 3*xspacing/2 xspacing/2 0 xspacing/2];
ys = [6+[2 3 3 2 1 1 2 3]]*yspacing;
p = patch(xs,ys,'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')

text(0.36*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     0.93*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(B)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')


% packed discs
circleAngles = linspace(0,2*pi+0.1,36);
dleft = 6*xspacing;
dbottom = 9.5*yspacing;
left = dleft + 0.5*xspacing;
bottom = dbottom-yspacing/2;
p = patch(left+xspacing/2*cos(circleAngles),bottom+2/sqrt(3)*yspacing/2*sin(circleAngles),'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')
left = dleft + 1.5*xspacing;
bottom = dbottom-xspacing/2+0.007;
p = patch(left+xspacing/2*cos(circleAngles),bottom+2/sqrt(3)*yspacing/2*sin(circleAngles),'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')
left = dleft + 0*xspacing;
bottom = dbottom-yspacing-xspacing/2+0.007;
p = patch(left+xspacing/2*cos(circleAngles),bottom+2/sqrt(3)*yspacing/2*sin(circleAngles),'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')
left = dleft + 1*xspacing;
bottom = dbottom-yspacing-xspacing/2+0.007;
p = patch(left+xspacing/2*cos(circleAngles),bottom+2/sqrt(3)*yspacing/2*sin(circleAngles),'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')
left = dleft + 2*xspacing;
bottom = dbottom-yspacing-xspacing/2+0.007;
p = patch(left+xspacing/2*cos(circleAngles),bottom+2/sqrt(3)*yspacing/2*sin(circleAngles),'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')
left = dleft + 0.5*xspacing;
bottom = dbottom-2*yspacing-xspacing/2+0.007;
p = patch(left+xspacing/2*cos(circleAngles),bottom+2/sqrt(3)*yspacing/2*sin(circleAngles),'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')
left = dleft + 1.5*xspacing;
bottom = dbottom-2*yspacing-xspacing/2+0.007;
p = patch(left+xspacing/2*cos(circleAngles),bottom+2/sqrt(3)*yspacing/2*sin(circleAngles),'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')

text(0.60*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     0.93*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(C)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')


%% rhombus
left = 3.5*xspacing;
xs = left-[xspacing 0 xspacing/2 3*xspacing/2 xspacing 0];
ys = 3*yspacing + [2 2 1 1 2 2]*yspacing;
p = patch(xs,ys,'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')

text(0.29*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     0.57*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(D)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')


%% wrapped torus
left = 4.5*xspacing;
bottom = 4*yspacing;
xs = left+[0 xspacing xspacing 0 0 xspacing];
ys = bottom+[yspacing yspacing 0 0 yspacing yspacing];
p = patch(xs,ys,'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')
xs = left+[0 xspacing xspacing 0 0 xspacing]-xspacing/2;
ys = bottom+[yspacing yspacing 0 0 yspacing yspacing]-yspacing;
p = patch(xs,ys,'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none','linestyle',':')

text(0.52*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     0.57*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(E)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')


%% two bumps on torus
left = 7*xspacing;
bottom = 4*yspacing;
xs = left+[0 xspacing xspacing 0 0 xspacing]-xspacing/2;
ys = bottom+[2*yspacing 2*yspacing 0 0 2*yspacing 2*yspacing]-yspacing;
p = patch(xs,ys,'k','linewidth',shapeLineWidth);
set(p,'edgecolor','k','facecolor','none')

text(0.75*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     0.57*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(F)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')

%% Blair et al 2007 Moire model
left = lefts(3)-0.15; bottom = bottoms(2)-0.28;
axes('position',[left bottom widths(2) heights(2)]);
set(gca,'cameraposition',[-2.8 -6. 11])
axis equal
set(gca,'xlim',[0 1],'ylim',[0 1])
axis off

[rs cs] = size(M);
winlen = 128;
rs = ceil(rs/winlen);
cs = ceil(cs/winlen);
% If making pdf with export_fig, can't use surface to plot this matrix because
% it is too big and breaks a string limit in the eps file format. Instead
% we break it into windows of data <2^16 bytes long.
for rind=1:rs
  for cind=1:cs
    if rind<max(rs) && cind<max(cs)
      dataWin = M((rind-1)*winlen+(1:winlen),(cind-1)*winlen+(1:winlen));
    elseif rind<max(rs) && cind==max(cs)
      dataWin = M((rind-1)*winlen+(1:winlen),(cind-1)*winlen+1:end);
    elseif rind==max(rs) && cind==max(cs)
      dataWin = M((rind-1)*winlen+1:end,(cind-1)*winlen+1:end);
    else
      dataWin = M((rind-1)*winlen+1:end,(cind-1)*winlen+(1:winlen));
    end
    surface('xdata',[(cind-1)/cs cind/cs],'ydata',[(rind-1)/rs rind/rs],'zdata',zeros(2),'cdata',dataWin,'facecolor','texturemap','edgecolor','none','cdatamapping','direct');
  end
end

text(.3*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     1.25*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(G)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')
   
   
%% Spatial interference of bands
left = lefts(3)-0.15; bottom = bottoms(1)-0.25;
axes('position',[left bottom widths(2) heights(2)])

surface('xdata',[0 1],'ydata',[0 1],'zdata',zeros(2),'cdata',thrt2,'facecolor','texturemap','edgecolor','k','cdatamapping','direct');

set(gca,'cameraposition',[-2.8 -6. 11])
axis equal
axis off

text(.3*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     1.25*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(H)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')
