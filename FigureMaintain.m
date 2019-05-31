% figure showing ways various models encode positions
% eric zilli - 20111005 - v1.0

% this figure is drawn very piecemeal. after the opening bit where
% the colormaps are made and all the data is loaded, the figure is
% made and then a separate section draws each plot. each plot is a
% separate cell so use cell titles or ctrl+up/down to navigate this file


%% Make colormaps
% Playing around with different colormaps:
% Colormap from blue to white to red
startColor = [0 0 1];
midColor = [1 1 1];
endColor = [1 0 0];
RWB = [linspace(startColor(1),midColor(1),32) linspace(midColor(1),endColor(1),32); linspace(startColor(2),midColor(2),32) linspace(midColor(2),endColor(2),32); linspace(startColor(3),midColor(3),32) linspace(midColor(3),endColor(3),32)]';

% Darker blue to white to darker red
startColor = [0 0 0.9];
midColor = [1 1 1];
endColor = [0.9 0 0];
mRWB = [linspace(startColor(1),midColor(1),32) linspace(midColor(1),endColor(1),32); linspace(startColor(2),midColor(2),32) linspace(midColor(2),endColor(2),32); linspace(startColor(3),midColor(3),32) linspace(midColor(3),endColor(3),32)]';

% Gray to white to gold
startColor = [0.07 0.23 0.28]/3;
midColor = [1 1 1];
endColor = [0.55 0.48 0.05]/1;
GWG = [linspace(startColor(1),midColor(1),32) linspace(midColor(1),endColor(1),32); linspace(startColor(2),midColor(2),32) linspace(midColor(2),endColor(2),32); linspace(startColor(3),midColor(3),32) linspace(midColor(3),endColor(3),32)]';

useColorMap = [flipud(GWG); RWB; mRWB];

% Plot a matrix scaled and shifted so that its zero point corresponds
% to the halfway point in a 64-item colormap. The first function draws
% the image, the second just returns the matrix.
imagepn = @(v)image(64*(v./2/abs(eps+max(abs([max(max(v)) min(min(v))])))+0.5));
imagepnc = @(v)(64*(v./2/abs(eps+max(abs([max(max(v)) min(min(v))])))+0.5));

%% Generate and load data for figure
% Approximate a wrapped-normal distribution to shape a bump for ring
% attractor plots (only accurate for small sigma and mu, X in [0,2*pi])
wnormpdf = @(X,mu,sigma)normpdf(X-4*pi,mu,sigma)+normpdf(X-2*pi,mu,sigma)+normpdf(X,mu,sigma)+normpdf(X+2*pi,mu,sigma)+normpdf(X+4*pi,mu,sigma);

% Position to illustrate in figure
% NB only about 75% of the elements in the figure will dynamically
% change if you change this. there's a bit of hard-coding here and there
currentPosition = [pi/2 2*pi]; % spatial phase, rad

nPos = 400;
pos1D = linspace(0,4*pi,nPos);

% Encode the spatial phase currentPosition
% in terms of the oscillations of each model:
Bu07baseline = cos(pos1D);
Bu07active1 = cos(pos1D-currentPosition(1));
Bu07active2 = cos(pos1D-currentPosition(2));

Bu08baseline = (0.5+0.5*cos(pos1D)).^50;
Bu08active1 = (0.5+0.5*cos(pos1D-currentPosition(1))).^50;
Bu08active2 = (0.5+0.5*cos(pos1D-currentPosition(2))).^50;

% This one has no baseline, but we need the third oscillator and it must
% be consistent in phase with the other two. Assuming the first two were
% at directions 0 and 60 and the third is at 120 degrees:
thirdPhase = [cos(2*pi/3) sin(2*pi/3)]*inv([cos(0) sin(0); cos(pi/3) sin(pi/3)])*currentPosition';
Ha08baseline = heaviside(cos(pos1D));
Ha08active1 = [0 heaviside(cos(pos1D-currentPosition(1)))];
Ha08active2 = [0 heaviside(cos(pos1D-currentPosition(2)))];
Ha08active3 = [heaviside(cos(pos1D-thirdPhase))];

% load variable vtraces
load data/ZilliHasselmo2010_voltage_traces.mat
Zi10windowlen = 27e3;
Zi10baseline = vtraces(:,1:Zi10windowlen)';
Zi10active1 = vtraces(:,round((2*pi-currentPosition(1))*Zi10windowlen/2/2/pi)+(1:Zi10windowlen))';
Zi10active2 = vtraces(:,round((2*pi-currentPosition(2))*Zi10windowlen/2/2/pi)+(1:Zi10windowlen))';

nDiscretizationCells = 16;
nModuloCells = nDiscretizationCells/2;
Ga07place = zeros(1,nDiscretizationCells);
Ga07place(nModuloCells + ceil(currentPosition(1)/2*pi)) = 1;
Ga07grid = Ga07place(1:nModuloCells)+Ga07place(nModuloCells+(1:nModuloCells));
Ga07place2 = zeros(2,nDiscretizationCells);
Ga07place2(1,nModuloCells + ceil(mod(currentPosition(1),2*pi)/2*pi)) = 1;
Ga07place2(2,nModuloCells + ceil(mod(currentPosition(2),2*pi)/2*pi)) = 1;
Ga07grid2 = Ga07place2(:,1:nModuloCells)+Ga07place2(:,nModuloCells+(1:nModuloCells));

load data/Fu06_WeightFigure_vars.mat;
f = flipud(Fu06_full_act);

load data/Gu07_WeightFigure_vars.mat;
A = flipud(Gu07_full_act);

load data/Bu09_WeightFigure_vars.mat;
s = flipud(Bu09_full_act);

load data/generalGridPattern.mat

%% Figure properties
colorFigures = 1;

leftMargin = 0.042;
bottomMargin1 = 0.03;
bottomMargin2 = 0.00;
nCols = 2;
nRows = 8;
lefts = leftMargin + [0 0.5 0.7125]; %0.95*(0:nCols-1)/nCols;
bottoms = bottomMargin1 + .91*(nRows-1:-1:0)/nRows;
bottoms(1) = 0.875;
width = 0.78/nCols;
height = .63/nRows;
ringHeight = 0.75/nRows;

nRowsCol2 = 5;
bottomsCol2 = bottomMargin2 + (nRowsCol2-1:-1:0)/nRowsCol2;
heightCol2 = 1.25*.65/nRowsCol2;

% shifts the phase-difference plot bottoms in the left column
leftColPhaseDiffBottomShift = 0.04;

% figure options:
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')

% size on paper:
widthOnPaper = 8.5; % cm
heightOnPaper = 16; % cm

figure('units','centimeters','position',[2 2 widthOnPaper heightOnPaper],'color','w');
set(gcf, 'renderer', 'painter')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [widthOnPaper heightOnPaper]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 widthOnPaper heightOnPaper]);

if colorFigures
  colors = 1*[useColorMap(1,:); useColorMap(64,:); useColorMap(64,:)*4];
else
  colors = [0 0 0; 0 0 0; 0 0 0];
end

ABlabelsX = -0.18;
xlabelY = -0.26;

%% 1D track with grid and spatial phases
left = lefts(1); bottom = bottoms(1)+.045;
axes('position',[left-width*.09 bottom width*1.1 height/1.5]);
surface('xdata',[0 1],'ydata',[0 1],'zdata',zeros(2),'cdata',128+round(imagepnc(gridPattern(1:73,58:69)')),'facecolor','texturemap','edgecolor','none','cdatamapping','direct');
surface('xdata',[0 0],'ydata',[0 1],'zdata',[0 4; 0 4],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
surface('xdata',[1 1],'ydata',[0 1],'zdata',[0 4; 0 4],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
surface('xdata',[0 1],'ydata',[0 0],'zdata',[4 4; 0 0],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
surface('xdata',[0 1],'ydata',[1 1],'zdata',[0 0; 4 4],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
line([0 1],[.5 .5],[0 0.00001],'color','k') 
set(gca,'cameraposition',[0 -5 75])
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
axis off

% "Linear coding"
text(.7*diff(get(gca,'xlim'))+min(get(gca,'xlim')),1.8*diff(get(gca,'ylim'))+min(get(gca,'ylim')),'Linear coding','horizontalalignment','center','fontsize',11)

% "Spatial phase"
text(0.4*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-1*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Spatial phase','horizontalalignment','center','fontsize',9)

% Panel label
text(0.14*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     1.91*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(A)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')

% Mark-up linear track with lines and dots to indicate how spatial phases
% map to positions
hold on;
plot3(0.5,0.5,0.001,'k.','markersize',8);
plot3(0.625,0.5,0.001,'ko','markersize',5);
plot3(1,0.5,0.001,'k.','markersize',8);
text(0*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.4*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'-360^\circ','horizontalalignment','center','fontsize',9)
text(0.5*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.4*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'0^\circ','horizontalalignment','center','fontsize',9)
text(5/8*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.4*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'90^\circ','horizontalalignment','center','fontsize',9)
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.4*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'360^\circ','horizontalalignment','center','fontsize',9)
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'box','off')
axis off

%% Burgess et al 2007 1D
left = lefts(1); bottom = bottoms(2)+leftColPhaseDiffBottomShift;
axes('position',[left bottom+height/2 width height/2]);
% Draw baseline oscillation
plot(pos1D,Bu07baseline,'color',colors(2,:));
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'ylim',[-1.1 1.1])
set(gca,'xlim',[0 max(pos1D)])
set(gca,'box','off')

% Category label - "Phase difference"
text(0.5*diff(get(gca,'xlim'))+min(get(gca,'xlim')),1.3*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Phase Difference','horizontalalignment','center')

% Draw active oscillation
axes('position',[left bottom width height/2]);
plot(pos1D,Bu07active1,'color',colors(1,:));
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'ylim',[-1.1 1.1])
set(gca,'xlim',[0 max(pos1D)])
set(gca,'box','off')

% Author label
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.35*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Burgess et al. 2007','horizontalalignment','right','fontsize',8)
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.8*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Giocomo et al. 2007','horizontalalignment','right','fontsize',8)


%% Burgess 2008 1D
left = lefts(1); bottom = bottoms(3)-0.005+leftColPhaseDiffBottomShift;
axes('position',[left bottom+height/2 width height/2]);
% Draw baseline oscillation
plot(pos1D,Bu08baseline,'color',colors(2,:));
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'ylim',[-0.1 1.1])
set(gca,'xlim',[0 max(pos1D)])
set(gca,'box','off')

% Draw active oscillation
axes('position',[left bottom width height/2]);
plot(pos1D,Bu08active1,'color',colors(1,:));
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'ylim',[-0.1 1.1])
set(gca,'xlim',[0 max(pos1D)])
set(gca,'box','off')

% Author label
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.35*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Burgess 2008','horizontalalignment','right','fontsize',8)

%% Hasselmo 2008
left = lefts(1); bottom = bottoms(3)-height/2;
axes('position',[left bottom width height/2]);
% plot(pos1D,Ha08active3,'color',[.55 .05 .48]);
plot(pos1D,Ha08baseline,'color',colors(2,:));
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'ylim',[-0.1 1.1])
set(gca,'xlim',[0 max(pos1D)])
set(gca,'box','off')

axes('position',[left bottom-height/2 width height/2]);
plot([0 pos1D],Ha08active1,'color',colors(1,:));
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'ylim',[-0.1 1.1])
set(gca,'xlim',[0 max(pos1D)])
set(gca,'box','off')

% axes('position',[left bottom width height/2]);
% plot([0 pos1D],Ha08active2,'color',colors(3,:));
% set(gca,'ydir','normal')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% set(gca,'xticklabel',[])
% set(gca,'ylim',[-0.1 1.1])
% set(gca,'xlim',[0 max(pos1D)])
% set(gca,'box','off')

% author label
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.35*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Hasselmo 2008','horizontalalignment','right','fontsize',8)


%% Zilli and Hasselmo 2010 1D
left = lefts(1); bottom = bottoms(5)+leftColPhaseDiffBottomShift;
axes('position',[left bottom+height/2 width height/2]);
% Draw baseline oscillation
plot(Zi10baseline,'color',colors(2,:));
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ylim',[-70 50])
set(gca,'xticklabel',[])
set(gca,'xlim',[0 length(Zi10baseline)])
set(gca,'box','off')

% Draw active oscillation
axes('position',[left bottom width height/2]);
plot(Zi10active1,'color',colors(1,:));
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ylim',[-70 50])
set(gca,'xticklabel',[])
set(gca,'xlim',[0 length(Zi10baseline)])
set(gca,'box','off')

% Author label
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.35*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Zilli and Hasselmo 2010','horizontalalignment','right','fontsize',8)

%% Blair et al. 2008 1D
left = lefts(1); bottom = bottoms(6)+0.01+leftColPhaseDiffBottomShift;
nRings = 7;
nCellsPerRing = 8;
ringFillColors = 1-27*(normpdf(1:nCellsPerRing,1,1)'*[.085 .085 .08] + normpdf(1:nCellsPerRing,nCellsPerRing+1,1)'*[.085 .085 .08]);
highlightCells = mod(ceil(linspace(nCellsPerRing/4,9*nCellsPerRing/4,nRings)),nCellsPerRing)+1;
% Baseline oscillator:
for ring=1:nRings
  axes('position',[left+(ring-1)*width/nRings bottom+height/2 width*1/nRings ringHeight/2]);
  drawSimpleRing(nCellsPerRing,1,0.7,0.5,highlightCells(ring),ringFillColors,ring==1);
  axis equal
  xlim([-1 1])
  ylim([-1 2])
  axis off
end
% Active oscillator:
ringFillColors = 1-9*(normpdf(1:nCellsPerRing,1,1)'*[.08 .08 .25] + normpdf(1:nCellsPerRing,nCellsPerRing+1,1)'*[.08 .08 .25]);
ringFillColors(ringFillColors<0) = 0;
highlightCells = mod(ceil(nCellsPerRing*currentPosition(1)/2/pi+linspace(nCellsPerRing/4,9*nCellsPerRing/4,nRings)),nCellsPerRing)+1;
for ring=1:nRings
  axes('position',[left+(ring-1)*width/nRings-width/2/nRings bottom-0.01 width*1.5/nRings ringHeight/2]);
  drawSimpleRing(nCellsPerRing,1,0.7,0.5,highlightCells(ring),ringFillColors,ring==1);
  axis equal
  xlim([-2 1])
  ylim([-2 1])
  axis off
end

% Author label
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),0.05*diff(get(gca,'ylim'))+min(get(gca,'ylim')),'Blair et al. 2008','horizontalalignment','right','fontsize',8)
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.33*diff(get(gca,'ylim'))+min(get(gca,'ylim')),'Welday et al. 2011','horizontalalignment','right','fontsize',8)

%% Gaussier et al. 2007
left = lefts(1); bottom = bottoms(7)-0.025;
% Draw coordinate cell
axes('position',[left+width/2+9*width/nDiscretizationCells/2 bottom+5*height/4 width height/3]);
axis off
colorInds = round(imagepnc(Ga07place));
for pind=1:1
  rectWidth = width/nDiscretizationCells;
  rectangle('position',[(pind-1)*rectWidth 0 rectWidth height/4],'edgecolor','k','facecolor',GWG(65-colorInds(11),:));
end
xlim([0 width])
% "Coordinate cell"
text(-0.6*diff(get(gca,'xlim'))+min(get(gca,'xlim')),.3*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Coordinate cell','horizontalalignment','left','fontsize',8)
set(gca,'xlim',[0 rectWidth*16],'ylim',[0 height/2])
colormap(useColorMap);
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')

% category label:
text(-0.25*diff(get(gca,'xlim'))+min(get(gca,'xlim')),1.16*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Spatial interference','horizontalalignment','center')

% Draw arrows from coordinate cell to discretization cells
axes('position',[0 0 1 1]);
axis off
for ar=1:nDiscretizationCells
  if ar==11
    annotation('arrow',[lefts(1)+width/2+11*width/nDiscretizationCells/2 left+2*ar*width/nDiscretizationCells/2]-width/nDiscretizationCells/2,[bottom+0.95*5*height/4 bottom+3.6*height/4],'linewidth',0.25,'headwidth',4,'headlength',4,'color','k')
  elseif ar>nDiscretizationCells/2
    annotation('arrow',[lefts(1)+width/2+11*width/nDiscretizationCells/2 left+2*ar*width/nDiscretizationCells/2]-width/nDiscretizationCells/2,[bottom+0.95*5*height/4 bottom+3.6*height/4],'linewidth',0.25,'headwidth',4,'headlength',4,'color',0.8*[1 1 1])
  else
    % annotations are always over text, so we use lines here so text can be on top
    line([lefts(1)+width/2+11*width/nDiscretizationCells/2 lefts(1)+2*ar*width/nDiscretizationCells/2]-width/nDiscretizationCells/2,[bottom+0.95*5*height/4 bottom+3.6*height/4],'linewidth',0.25,'color',0.8*[1 1 1])
  end
end
set(gca,'xlim',[0 1],'ylim',[0 1])


% Draw modulo/stripe cells
axes('position',[left+width/4 bottom-0.005+0.01 width/2 height/5]);
colormap(useColorMap);
colorInds = round(imagepnc(Ga07grid));
for pind=1:nModuloCells
  rectWidth = width/nModuloCells;
  rectangle('position',[(pind-1)*rectWidth 0 rectWidth height/4],'edgecolor','k','facecolor',GWG(65-colorInds(pind),:));
end
xlim([0 width])
axis off
set(gca,'ytick',[])
set(gca,'ydir','normal')
set(gca,'xticklabel',[])
set(gca,'box','off')

% Arrows from discretization cells to modulo cells
for ar=1:nModuloCells
  annotation('arrow',[lefts(1)+2*ar*width/nDiscretizationCells/2 left+width/4+2*ar*width/nDiscretizationCells/2]-width/nDiscretizationCells/2,[bottom+.8*3*height/4 bottom+1.15*height/4],'linewidth',0.25,'headwidth',4,'headlength',4,'color',0.8*[1 1 1])
end
for ar=(nModuloCells+1):2*nModuloCells
  if ar==11
    annotation('arrow',[lefts(1)+2*ar*width/nDiscretizationCells/2 left+width/4+2*(ar-nModuloCells)*width/nDiscretizationCells/2]-width/nDiscretizationCells/2,[bottom+.8*3*height/4 bottom+1.15*height/4],'linewidth',0.25,'headwidth',4,'headlength',4,'color','k')
  else
    annotation('arrow',[lefts(1)+2*ar*width/nDiscretizationCells/2 left+width/4+2*(ar-nModuloCells)*width/nDiscretizationCells/2]-width/nDiscretizationCells/2,[bottom+.8*3*height/4 bottom+1.15*height/4],'linewidth',0.25,'headwidth',4,'headlength',4,'color',0.8*[1 1 1])
  end
end

% Draw discretization cells
axes('position',[left bottom+2.7*height/4 width height/5]);
axis off
colorInds = round(imagepnc(Ga07place));
for pind=1:nDiscretizationCells
  rectWidth = width/nDiscretizationCells;
  rectangle('position',[(pind-1)*rectWidth 0 rectWidth height/4],'edgecolor','k','facecolor',GWG(65-colorInds(pind),:));
end
xlim([0 width])

% "Discretization cells"
text(-0.04*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.7*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),100,'cells','horizontalalignment','left','fontsize',8)
text(-0.04*diff(get(gca,'xlim'))+min(get(gca,'xlim')),1.5*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),100,'Discretization','horizontalalignment','left','fontsize',8)

% "Modulo cells"
text(-0.04*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-2.4*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Modulo','horizontalalignment','left','fontsize',8)
text(-0.04*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-3.4*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'cells','horizontalalignment','left','fontsize',8)

% author label:
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-4*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Gaussier et al. 2007','horizontalalignment','right','fontsize',8)
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-5.1*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Hasselmo and Brandon 2008','horizontalalignment','right','fontsize',8)


colormap(useColorMap);
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')

%% Mhatre et al. 2010 1D
left = lefts(1); bottom = bottoms(8)-0.034;
nRings = 7;
nCellsPerRing = 8;
ringFillColors = 1-27*(normpdf(1:nCellsPerRing,1,1)'*[.085 .085 .08] + normpdf(1:nCellsPerRing,nCellsPerRing+1,1)'*[.085 .085 .08]);
for ring=1:nRings
  axes('position',[left+(ring-1)*width/nRings bottom+height/2 width*1/nRings ringHeight/3]);
  drawSimpleRing(nCellsPerRing,1,0.7,0.5,1+ceil(nCellsPerRing*(0.25+currentPosition(1)/2/pi)),ringFillColors,0);
  axis off
end

% author labels:
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.35*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Mhatre et al. 2010','horizontalalignment','right','fontsize',8)
text(1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.85*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Navratilova et al. 2011','horizontalalignment','right','fontsize',8)

% %% Divider between Col 1 and 2
% annotation('line',[lefts(2) lefts(2)]-0.024,[bottoms(1)-0.04 bottoms(7)]);

%% Column 2
% 2D Plot
%% 2D arena with grid and spatial phases
left = lefts(2)-0.1; bottom = bottoms(1)-0.08;
axes('position',[left bottom-heightCol2/2 width*1.4 heightCol2*1.4]);
surface('xdata',[0 1],'ydata',[0 1],'zdata',zeros(2),'cdata',128+round(imagepnc(flipud(gridPattern))),'facecolor','texturemap','edgecolor','none','cdatamapping','direct');
surface('xdata',[0 0],'ydata',[0 1],'zdata',[0 4; 0 4],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
surface('xdata',[1 1],'ydata',[0 1],'zdata',[0 4; 0 4],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
surface('xdata',[0 1],'ydata',[0 0],'zdata',[4 4; 0 0],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
surface('xdata',[0 1],'ydata',[1 1],'zdata',[0 0; 4 4],'facecolor','w','edgecolor','k','DiffuseStrength',1,'edgelighting','phong','facelighting','phong')
% Lines indicating spatial phases:
phaseLinesShift = [.3 0.099];
zeroField = [.31 .9]+phaseLinesShift;
field1 = [.012 .35]+phaseLinesShift;
field2 = [.64 .36]+phaseLinesShift;
direction1line = [zeroField; field1]';
slope1 = diff(direction1line');
halfpoint1 = direction1line(:,1) + 1/2*slope1';
position1 = direction1line(:,1) + 5/8*slope1';
endpoint1 = position1 + -.4*[slope1(2); -slope1(1)];
direction2line = [zeroField; field2]';
slope2 = diff(direction2line');
halfpoint2 = direction2line(:,1) + 1/2*slope2';
position2 = direction2line(:,1) + 4/8*slope2';
endpoint2 = position2 + -.5*[-slope2(2); slope2(1)];
line(direction1line(1,:),direction1line(2,:),[0 0.00001],'color','k')
line(direction2line(1,:),direction2line(2,:),[0 0.00001],'color','k')
line([position1(1) endpoint1(1)],[position1(2) endpoint1(2)],[0 0.00001],'color','k','linestyle',':','linewidth',1.2)
line([position2(1) endpoint2(1)],[position2(2) endpoint2(2)],[0 0.00001],'color','k','linestyle',':','linewidth',1.2)

hold on;
plot3(direction1line(1,1),direction1line(2,1),0.001,'k.','markersize',8);
plot3(direction1line(1,2),direction1line(2,2),0.001,'k.','markersize',8);
plot3(direction2line(1,2),direction2line(2,2),0.001,'k.','markersize',8);
plot3(position1(1),position1(2),0.001,'ko','markersize',5);
plot3(position2(1),position2(2),0.001,'ko','markersize',5);
plot3(halfpoint1(1),halfpoint1(2),0.001,'k.','markersize',8);
plot3(halfpoint2(1),halfpoint2(2),0.001,'k.','markersize',8);

text(0.32+phaseLinesShift(1), 0.83+phaseLinesShift(2), 4.1, '0^\circ')
text(0.46+phaseLinesShift(1), 0.79+phaseLinesShift(2), .1, '360^\circ')
text(0.02+phaseLinesShift(1), 0.75+phaseLinesShift(2), .1, '360^\circ')
set(gca,'cameraposition',[0 -1 36])
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
axis off

% "Planar coding"
text(0.88*diff(get(gca,'xlim'))+min(get(gca,'xlim')),1.56*diff(get(gca,'ylim'))+min(get(gca,'ylim')),'Planar coding','horizontalalignment','center','fontsize',11)

% panel label
text(0.33*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
     1.76*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
     '(B)',...
     'FontSize',10,...
     'FontWeight','bold',...
     'HorizontalAlignment','center')

%% Place-driven model
left = lefts(2)-0.015; bottom = bottoms(4)-0.01; %bottomsCol2(2)-0.15;%0.09;
axes('position',[left-0.02 bottom 1.2*[width heightCol2]]);

% gridEllipse = [left-0.01 bottom+0.005 .02 .01];
gridEllipse = [left-0.05 bottom+0.052 .02 .01];
gridPos = [gridEllipse(1)+gridEllipse(3)/2 gridEllipse(2)+gridEllipse(4)/2];

% Top grid field arrow/lines
lineStart = [left-0.0185 bottom+0.127];
lineLen = 0.9;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('arrow',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linewidth',0.5,'headwidth',6,'headlength',6)

lineStart = [left-0.006 bottom+0.118];
lineLen = 0.5;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('line',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linestyle',':')

% Middle grid field arrow/lines
lineStart = [left+0.03 bottom+0.0905];
lineLen = 0.85;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('arrow',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linewidth',0.5,'headwidth',6,'headlength',6)

lineStart = [left+0.02 bottom+0.0985];
lineLen = 0.5;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('line',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linestyle',':')


lineStart = [left+0.043 bottom+0.081];
lineLen = 0.46;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('line',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linestyle',':')

% lineStart = [left+0.032 bottom+0.051];
% lineLen = 0.5;
% lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
% annotation('line',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linestyle',':')

% Bottom grid field arrow/lines
lineStart = [left+0.082 bottom+0.052];
lineLen = 0.88;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('arrow',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linewidth',0.5,'headwidth',6,'headlength',6)

lineStart = [left+0.07 bottom+.061];
lineLen = 0.5;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('line',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linestyle',':')

lineStart = [left+0.093 bottom+.042];
lineLen = 0.5;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('line',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linestyle',':')


% Grayed arrow/lines
lineStart = [left+0.076 bottom+0.124];
lineLen = 0.29;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('line',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linestyle',':','color',[0.6 0.6 0.6])

lineStart = [left+0.138 bottom+0.115];
lineLen = 0.33;
lineEnd = [lineStart(1)+lineLen*(gridPos(1)-lineStart(1)) lineStart(2)+lineLen*(gridPos(2)-lineStart(2))];
annotation('line',[lineStart(1) lineEnd(1)],[lineStart(2) lineEnd(2)],'linestyle',':','color',[0.6 0.6 0.6])

annotation('ellipse',gridEllipse,'facecolor',useColorMap(127,:))
surface('xdata',[0 1],'ydata',[0 1],'zdata',zeros(2),'cdata',round(imagepnc(flipud(gridPattern)>1)),'facecolor','texturemap','edgecolor','k','cdatamapping','direct');
set(gca,'cameraposition',[-0.77 -2.5 16])
axis off
% citations:
text(0.82*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.50*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Kropff and Treves 2008','horizontalalignment','right','fontsize',8)

% "Place cells"
text(0.5*diff(get(gca,'xlim'))+min(get(gca,'xlim')),1.11*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Place cells','horizontalalignment','right','fontsize',8,'rotation',14)

% "Grid cells"
text(-0.12*diff(get(gca,'xlim'))+min(get(gca,'xlim')),0.2*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Grid cell','horizontalalignment','right','fontsize',8)

% category label:
text(0.8*diff(get(gca,'xlim'))+min(get(gca,'xlim')),1.3*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Place-driven','horizontalalignment','center')


col4SizeScale = 1.7/2;
col4Lefts = -0.02;
%% Toroidal network/rectangular grid
% we'll generate a simple rectangular grid pattern:
left = lefts(2)-0.1+col4Lefts; bottom = bottoms(6);%+0*height/2;
axes('position',[left bottom col4SizeScale*[width heightCol2]]);
n = 11;
[X,Y] = ndgrid(linspace(-1,1,n),linspace(-1,1,n));
R = @(x)(x.*(x>0.5));
rect = fftshift(R(exp(-sqrt(X.^2+Y.^2)/1)));
surface('xdata',[0 1],'ydata',[0 1],'zdata',zeros(2),'cdata',0+round(imagepnc(rect)),'facecolor','texturemap','edgecolor','k','cdatamapping','direct');
axis equal
axis off
% author label:
text(0.9*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.2*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),{'Toroidal','network'},'horizontalalignment','right','fontsize',8)

%% Fuhs and Touretzky model
left = lefts(2)+0.18+col4Lefts; bottom = bottoms(6);%+height/2;%+1.5*height/2;
axes('position',[left bottom col4SizeScale*[width heightCol2]]);
surface('xdata',[0 1],'ydata',[0 1],'zdata',zeros(2),'cdata',0+round(imagepnc(f)),'facecolor','texturemap','edgecolor','k','cdatamapping','direct');
line([0 1],1+[1e-3 1e-3],'color','k')
rectWidth = 1/62;
rectHeight = 1/62;
rectangle('position',[28*rectWidth 35*rectHeight rectWidth rectHeight],'edgecolor','r')
axis equal
axis off
% author label:
text(0.9*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.2*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),{'Fuhs and','Touretzky 2006'},'horizontalalignment','right','fontsize',8)
% text(0.8*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.07*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Fuhs and Touretzky 2006','horizontalalignment','right','fontsize',8)

% category label:
text(0*diff(get(gca,'xlim'))+min(get(gca,'xlim')),1.15*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),'Population Activity','horizontalalignment','center','fontsize',10)

%% Guanella model
left = lefts(2)-0.1+col4Lefts; bottom = bottoms(8)+height/3;
axes('position',[left bottom col4SizeScale*[width heightCol2]]);
surface('xdata',[0 1],'ydata',[0 1],'zdata',zeros(2),'cdata',0+round(imagepnc(A)),'facecolor','texturemap','edgecolor','k','cdatamapping','direct');
rectWidth = 1/10;
rectHeight = 1/9;
rectangle('position',[3*rectWidth 5*rectHeight rectWidth rectHeight],'edgecolor','r')
axis equal
axis off
% author label:
text(0.9*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.2*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),{'Guanella et','al. 2007'},'horizontalalignment','right','fontsize',8)

%% Burak and Fiete model
left = lefts(2)+0.18+col4Lefts; bottom = bottoms(8)+height/3;
axes('position',[left bottom col4SizeScale*[width heightCol2]]);
surface('xdata',[0 1],'ydata',[0 1],'zdata',zeros(2),'cdata',0+round(imagepnc(s)),'facecolor','texturemap','edgecolor','k','cdatamapping','direct');
rectWidth = 1/128;
rectHeight = 1/128;
rectangle('position',[58*rectWidth 60*rectHeight rectWidth rectHeight],'edgecolor','r')
axis equal
axis off
% author label:
text(0.9*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-0.2*(max(get(gca,'ylim'))-min(get(gca,'ylim')))+min(get(gca,'ylim')),{'Burak and','Fiete 2009'},'horizontalalignment','right','fontsize',8)

