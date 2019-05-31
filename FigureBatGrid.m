% bat-grid--like grid cell activity from oscillatory interference models
% eric zilli - 20111108 - v1.0

% Columns:        Burgess et al. 2007   Burgess 2008    Hasselmo 2008
% Rows:
% Baseline = 5    [rate map]   [ISIs]
% Baseline = 0.5
% Baseline = 0

caption = {'With a 5 Hz baseline, most',...
           '(but not all) temporal',...
           'interference models show grid',...
           'cell firing rate modulation at 5 Hz.',...
           'With the precession mechanisms',...
           'described in these papers, the',...
           'baseline can go low enough that',...
           'theta modulation disappears while',...
           'baseline modulation remains.',...
           'These even work with a baseline',...
           'of 0 Hz, which turns them into',...
           'spatial interference models.'};

         
methods = {'200 s simulations. Hafting et al. 2005 trajectory. Unbiased autocorrelations. Gaussian smoothed rate maps,  \sigma = 3.3 cm. For code/comments/questions: zilli@bu.edu. Eric Zilli.'};

% Create the figure we're making
widthOnPaper = 9;
heightOnPaper = 5;
% figure options:
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
figure('units','inches','position',[1 1 widthOnPaper heightOnPaper],'color','w');
set(gcf, 'renderer', 'painter')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [widthOnPaper heightOnPaper]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 widthOnPaper heightOnPaper]);

% Variables for positioning plots
leftMargin = 0.25;
bottomMargin = 0.1;
nCols = 6;
nRows = 3;
lefts = leftMargin + 0.77*(0:nCols-1)/nCols;
bottoms = bottomMargin + .8*(0:nRows-1)/nRows;
bottoms2 = 0.01 + bottomMargin + .8*(0:nRows-1)/nRows;
width = 0.55/nCols;
height = .7/nRows;
height2 = 0.2/nRows;

% figure options:
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')

% kernel to smooth rate maps:
gaussian = fspecial('gaussian',[5 5],1);

% length simulations will run (set separately inside each)
simdur = 200; % s

% X axis width for the autocorrelations
xcorrWidth = 0.5; % s
xcorrWidth2 = 5; % s

% Position (m) corresponding to each spatial bin of the rate map
posBins = linspace(-1,1,60);

% line spacing of caption text
linespacing = .33;

% Relative (to rows) y position of caption text
captionY = 0.7;
captionY2 = 0.8;
captionY3 = 0.65;

% Relative (to rate maps) position of reference names over top row of plots
referenceTextX = 1.2;
referenceTextY = -0.7;

%% Run Burgess et al. 2007 script, not bat mode, hide figures
[dt spikes occupancy spikeTimes] = BurgessEtAl2007_bat(0,0);

% Calculate and plot smoothed rate map
axes('position',[lefts(1) bottoms(3) width height]);
imagesc(posBins,posBins,conv2(gaussian,spikes./(occupancy+eps)));
axis square;
title('Rate map')

% Add labels
text(-0.85*diff(get(gca,'xlim'))+min(get(gca,'xlim')),0.5*diff(get(gca,'ylim'))+min(get(gca,'ylim')),{'Rodent','f = 5 Hz'},'fontsize',12,'horizontalalignment','right')
text(referenceTextX*diff(get(gca,'xlim'))+min(get(gca,'xlim')),referenceTextY*diff(get(gca,'ylim'))+min(get(gca,'ylim')),{'Burgess et al. 2007'},'fontsize',11,'horizontalalignment','center','verticalalignment','middle')

% Add caption
for ind=1:4
  text(-0.25*diff(get(gca,'xlim'))+min(get(gca,'xlim')),linespacing*ind-captionY*diff(get(gca,'ylim'))+min(get(gca,'ylim')),caption{ind},'fontsize',9,'horizontalalignment','right')
end

% Calculate and plot first autocorrelation
axes('position',[lefts(2) bottoms2(3)+height2*2 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth:dt:xcorrWidth,xcorr(spikes,xcorrWidth/dt,'unbiased'))
title({'Grid cell spike','train autocorrelation'})

% Calculate and plot second autocorrelation
axes('position',[lefts(2) bottoms2(3)+height2/5 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth2:dt:xcorrWidth2,xcorr(spikes,xcorrWidth2/dt,'unbiased'))

%% Run Burgess et al. 2007 script, bat mode 1 (0.5 Hz baseline), hide figures
[dt spikes occupancy spikeTimes] = BurgessEtAl2007_bat(0,0,0.5);

% Calculate and plot smoothed rate map
axes('position',[lefts(1) bottoms(2) width height]);
imagesc(posBins,posBins,conv2(gaussian,spikes./(occupancy+eps)));
axis square;
% title('Rate map')

% Add caption
for ind=5:9
  text(-0.25*diff(get(gca,'xlim'))+min(get(gca,'xlim')),linespacing*(ind-4)-captionY2*diff(get(gca,'ylim'))+min(get(gca,'ylim')),caption{ind},'fontsize',9,'horizontalalignment','right')
end

% Add labels
text(-0.85*diff(get(gca,'xlim'))+min(get(gca,'xlim')),0.5*diff(get(gca,'ylim'))+min(get(gca,'ylim')),{'Bat?','f = 0.5 Hz'},'fontsize',12,'horizontalalignment','right')

% Calculate and plot first autocorrelation
axes('position',[lefts(2) bottoms2(2)+height2*2 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth:dt:xcorrWidth,xcorr(spikes,xcorrWidth/dt,'unbiased'))
% title({'Grid cell spike','train autocorrelation'})

% Calculate and plot second autocorrelation
axes('position',[lefts(2) bottoms2(2)+height2/5 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth2:dt:xcorrWidth2,xcorr(spikes,xcorrWidth2/dt,'unbiased'))
% title({'Grid cell spike','train autocorrelation'})

%% Run Burgess et al. 2007 script, bat mode 2, hide figures
[dt spikes occupancy spikeTimes] = BurgessEtAl2007_bat(1,0);

% Calculate and plot smoothed rate map
axes('position',[lefts(1) bottoms(1) width height]);
imagesc(posBins,posBins,conv2(gaussian,spikes./(occupancy+eps)));
axis square;
xlabel('Position (m)')

% Add caption
for ind=10:12
  text(-0.25*diff(get(gca,'xlim'))+min(get(gca,'xlim')),linespacing*(ind-9)-captionY3*diff(get(gca,'ylim'))+min(get(gca,'ylim')),caption{ind},'fontsize',9,'horizontalalignment','right')
end

% Add labels
text(-0.85*diff(get(gca,'xlim'))+min(get(gca,'xlim')),0.5*diff(get(gca,'ylim'))+min(get(gca,'ylim')),{'Bat?','f = 0 Hz'},'fontsize',12,'horizontalalignment','right')

% Calculate and plot first autocorrelation
axes('position',[lefts(2) bottoms2(1)+height2*2 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth:dt:xcorrWidth,xcorr(spikes,xcorrWidth/dt,'unbiased'))
% title({'Grid cell spike','train autocorrelation'})

% Calculate and plot second autocorrelation
axes('position',[lefts(2) bottoms2(1)+height2/5 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth2:dt:xcorrWidth2,xcorr(spikes,xcorrWidth2/dt,'unbiased'))
xlabel('Time (s)')

% "Methods" text
text(-4.1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),-1.45*diff(get(gca,'ylim'))+min(get(gca,'ylim')),methods{1},'fontsize',8,'horizontalalignment','left')


%% Run Burgess 2008 script, not bat mode, hide figures
[dt spikes occupancy spikeTimes] = Burgess2008_bat(0,0);

% Calculate and plot smoothed rate map
axes('position',[lefts(3) bottoms(3) width height]);
imagesc(posBins,posBins,conv2(gaussian,spikes./(occupancy+eps)));
axis square;
title('Rate map')

% Add labels
text(referenceTextX*diff(get(gca,'xlim'))+min(get(gca,'xlim')),referenceTextY*diff(get(gca,'ylim'))+min(get(gca,'ylim')),{'Burgess 2008'},'fontsize',11,'horizontalalignment','center','verticalalignment','middle')


% Calculate and plot first autocorrelation
axes('position',[lefts(4) bottoms2(3)+height2*2 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth:dt:xcorrWidth,xcorr(spikes,xcorrWidth/dt,'unbiased'))
title({'Grid cell spike','train autocorrelation'})

% Calculate and plot second autocorrelation
axes('position',[lefts(4) bottoms2(3)+height2/5 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth2:dt:xcorrWidth2,xcorr(spikes,xcorrWidth2/dt,'unbiased'))

%% Run Burgess 2008 script, bat mode 1 (0.5 Hz baseline), hide figures
[dt spikes occupancy spikeTimes] = Burgess2008_bat(0,0,0.5);

% Calculate and plot smoothed rate map
axes('position',[lefts(3) bottoms(2) width height]);
imagesc(posBins,posBins,conv2(gaussian,spikes./(occupancy+eps)));
axis square;
% title('Rate map')

% Calculate and plot first autocorrelation
axes('position',[lefts(4) bottoms2(2)+height2*2 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth:dt:xcorrWidth,xcorr(spikes,xcorrWidth/dt,'unbiased'))
% title({'Grid cell spike','train autocorrelation'})

% Calculate and plot second autocorrelation
axes('position',[lefts(4) bottoms2(2)+height2/5 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth2:dt:xcorrWidth2,xcorr(spikes,xcorrWidth2/dt,'unbiased'))

%% Run Burgess 2008 script, bat mode 2, hide figures
[dt spikes occupancy spikeTimes] = Burgess2008_bat(1,0);

% Calculate and plot smoothed rate map
axes('position',[lefts(3) bottoms(1) width height]);
imagesc(posBins,posBins,conv2(gaussian,spikes./(occupancy+eps)));
axis square;
xlabel('Position (m)')

% Calculate and plot first autocorrelation
axes('position',[lefts(4) bottoms2(1)+height2*2 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth:dt:xcorrWidth,xcorr(spikes,xcorrWidth/dt,'unbiased'))
% title({'Grid cell spike','train autocorrelation'})

% Calculate and plot second autocorrelation
axes('position',[lefts(4) bottoms2(1)+height2/5 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth2:dt:xcorrWidth2,xcorr(spikes,xcorrWidth2/dt,'unbiased'))
xlabel('Time (s)')

%% Run Hasselmo 2008 script, not bat mode, hide figures
[dt spikes occupancy spikeTimes] = Hasselmo2008_bat(0,0);

% Calculate and plot smoothed rate map
axes('position',[lefts(5) bottoms(3) width height]);
imagesc(posBins,posBins,conv2(gaussian,spikes./(occupancy+eps)));
axis square;
title('Rate map')

% Add labels
text(referenceTextX*diff(get(gca,'xlim'))+min(get(gca,'xlim')),referenceTextY*diff(get(gca,'ylim'))+min(get(gca,'ylim')),{'Hasselmo 2008'},'fontsize',11,'horizontalalignment','center','verticalalignment','middle')


% Calculate and plot first autocorrelation
axes('position',[lefts(6) bottoms2(3)+height2*2 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth:dt:xcorrWidth,xcorr(spikes,xcorrWidth/dt,'unbiased'))
title({'Grid cell spike','train autocorrelation'})

% Calculate and plot second autocorrelation
axes('position',[lefts(6) bottoms2(3)+height2/5 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth2:dt:xcorrWidth2,xcorr(spikes,xcorrWidth2/dt,'unbiased'))

%% Run Hasselmo 2008 script, bat mode 1 (0.5 Hz baseline), hide figures
[dt spikes occupancy spikeTimes] = Hasselmo2008_bat(0,0,0.5);

% Calculate and plot smoothed rate map
axes('position',[lefts(5) bottoms(2) width height]);
imagesc(posBins,posBins,conv2(gaussian,spikes./(occupancy+eps)));
axis square;
% title('Rate map')

% Calculate and plot first autocorrelation
axes('position',[lefts(6) bottoms2(2)+height2*2 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth:dt:xcorrWidth,xcorr(spikes,xcorrWidth/dt,'unbiased'))
% title({'Grid cell spike','train autocorrelation'})

% Calculate and plot second autocorrelation
axes('position',[lefts(6) bottoms2(2)+height2/5 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth2:dt:xcorrWidth2,xcorr(spikes,xcorrWidth2/dt,'unbiased'))

%% Run Hasselmo 2008 script, bat mode 2, hide figures
[dt spikes occupancy spikeTimes] = Hasselmo2008_bat(1,0);

% Calculate and plot smoothed rate map
axes('position',[lefts(5) bottoms(1) width height]);
imagesc(posBins,posBins,conv2(gaussian,spikes./(occupancy+eps)));
axis square;
xlabel('Position (m)')

% Calculate and plot first autocorrelation
axes('position',[lefts(6) bottoms2(1)+height2*2 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth:dt:xcorrWidth,xcorr(spikes,xcorrWidth/dt,'unbiased'))
% title({'Grid cell spike','train autocorrelation'})

% Calculate and plot second autocorrelation
axes('position',[lefts(6) bottoms2(1)+height2/5 width height2]);
spikes = zeros(1,ceil(simdur/dt));
spikes(round(spikeTimes/dt)) = 1;
plot(-xcorrWidth2:dt:xcorrWidth2,xcorr(spikes,xcorrWidth2/dt,'unbiased'))
xlabel('Time (s)')


