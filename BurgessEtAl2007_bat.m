% Burgess, Barry, O'Keefe 2007's abstract oscillatory interference model
% eric zilli - 20111106 - v1.0
% 
% This is the bat implementation of Burgess et al. 2007, accounting for
% Yartsev et al. 2011's report that grid cells in bats that lack theta
% activity in the LFP and in the ISIs/spiking autocorrelation.
% 
% We use precession mechanism A to produce non-baseline-modulated
% firing and precession mechanism C to produce baseline-modulated
% firing with a non-zero baseline. Read about these in the comments
% for the script BurgessEtAl2007_precession.m from which this script
% derived which are included below.
%
% Will the Joker's dastardly deeds destroy Gotham?
% Tune in next .m in the same bat comments, same bat time!
%
% NOTES FROM ORIGINAL SCRIPT (BurgessEtAl2007_precession.m):
% This variation uses equation (6) from the manuscript which is reported
% in the paper to produce proper phase precession.
% Quoting them "This model is similar to the basic model
%   above with n=6 but adds pairs of dendritic oscillators with
%   opposing preferred directions before multiplying the three
%   resulting interference patterns together..."
% 
% Unfortunately, no equation was provided so we can only guess as to
% exactly what that means. In particular, no mention is made
% of a baseline oscillation (or to thresholding), so either:
% A. No baseline oscillation is used in this variation. This seems to
%   produce phase precession but the phase wanders.
% B. The baseline oscillation is multiplied in as a fourth term (along with
%   the three summed pairs of opposing oscillators). This seems to produce
%   less phase wandering because the baseline masks it out, but it causes
%   the cell to not fire on many passes through where fields should be.
% C. The baseline oscillation is included in each of the three sums. This
%   seems to work fairly well.
%
% The variations have not been extensively tested, though, and I'm still
% not quite sure which one he intended.
%
% We also include another solution, subtracting the dendrite phases before
% taking the cosine:
% threshold(cos(phi_0)+cos(phi_1-phi_2))
% This produces correct path integration, but it is not at all clear how
% that could be carried out biologically.
%
% Burgess 2008 gave a more straightforward fix, which is to block (i.e. zero
% out) the output of the VCO when the animal's current head direction
% is more than 90 degrees from the VCO's preferred direction, but to
% allow the VCO to change its frequency regardless of the animal's heading.
%
% See BurgessEtAl2007.m for more information about this model.
%
% Note: This code is not necessarily optimized for speed, but is meant as
% a transparent implementation of the model as described in the manuscript.
%
% This code is released into the public domain. Not for use in skynet.

function [dt spikes occupancy spikeTimes] = BurgessEtAl2007_bat(varargin)

doSubtractPhases = 1;

% This lets us turn bat mode (base freq = 0) on or off
batMode = 0;
if ~isempty(varargin)
  batMode = varargin{1};
end
if batMode
  % See above for the three precession variations
  precessionVariation = 'A'; % 'A', 'B', or 'C'
  baseFreq = 0; % Hz
  % NB. If baseFreq = 0 then baseline oscillation stays as cos(0) = 1
  % so adding it in using 'C' will change the thresholding
else
  % See above for the three precession variations
  precessionVariation = 'C'; % 'A', 'B', or 'C'
  baseFreq = 5; % Hz
end

% second parameter says whether or not to plot figures:
showFigs = 1;
if length(varargin)>1
  showFigs = varargin{2};
end

% third parameter is baseline frequency override
if length(varargin)>2
  baseFreq = varargin{3};
end

% if >0, plots the sheet of activity during the simulation on every livePlot'th step
livePlot = 0;
if showFigs
  livePlot = 0;
end

% if =0, just give constant velocity. if =1, load trajectory from disk
useRealTrajectory = 1;
constantVelocity = 1*[.5; 0*0.5]; % m/s

%% Simulation parameters
dt = .01; % time step, s
simdur = 200; % total simulation time, s
tind = 1; % time step number for indexing
t = 0; % simulation time variable, s
x = 0; % position, m
y = 0; % position, m

%% Model parameters
ncells = 1;
% Directional preference of each dendrite (this also sets the number of dendrites)
dirPreferences = (0:5)*pi/3;
% This will let us add/subtract dendritic values with opposite direction preferences:
if doSubtractPhases
  oppositeDendrites = [1 0 0 -1  0  0;
                       0 1 0  0 -1  0;
                       0 0 1  0  0 -1];
else
  oppositeDendrites = [1 0 0 1 0 0;
                       0 1 0 0 1 0;
                       0 0 1 0 0 1];
end
% Scaling factor relating speed to oscillator frequencies
% NB paper uses 0.05*2pi rad/cm. But we do the conversion to rad later,
% leaving 0.05 1/cm = 5 1/m which produces very tight field spacing. For cosmetic
% purposes for the trajectory we use here, we'll use beta = 2.
beta = 2; % Hz/(m/s) 

if doSubtractPhases
  spikeThreshold = 0.3;
else
  if precessionVariation=='A'
    spikeThreshold = 1; 2.5;
  elseif precessionVariation=='B'
    spikeThreshold = 0.3;
  elseif precessionVariation=='C'
    spikeThreshold = 5;
  else
    error('No such precession option.')
  end
end


%% History variables
speed = zeros(1,ceil(simdur/dt));
curDir = zeros(1,ceil(simdur/dt));
vhist = zeros(1,ceil(simdur/dt));
fhist = zeros(1,ceil(simdur/dt));

%% Firing field plot variables
nSpatialBins = 60;
minx = -1; maxx = 1; % m
miny = -1; maxy = 1; % m
% minx = -0.90; maxx = 0.90; % m
% miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);

spikeTimes = zeros(1,1000);
spikeCoords = zeros(1000,2);
spikePhases = zeros(1,1000);
spikeind = 1;

%% Initial conditions
% Oscillators will start at phase 0:
dendritePhases = zeros(1,length(dirPreferences)); % rad
basePhase = 0; % rad

%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w','name','Activity of one cell');
  if useRealTrajectory
    set(h,'position',[520 378 1044 420])
  end
  gaussian = fspecial('gaussian',[5 5],1); % to smooth rate map
  drawnow
end

%% Possibly load trajectory from disk
if useRealTrajectory
  load data/HaftingTraj_centimeters_seconds.mat;
  % interpolate down to simulation time step
  pos = [interp1(pos(3,:),pos(1,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(2,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(3,:),0:dt:pos(3,end))];
  pos(1:2,:) = pos(1:2,:)/100; % cm to m
  
  vels = [diff(pos(1,:)); diff(pos(2,:))]/dt; % m/s
  x = pos(1,1); % m
  y = pos(2,1); % m
end

%% !! Main simulation loop
fprintf('Simulation starting. Press ctrl+c to end...\n')
while t<simdur
  tind = tind+1;
  t = dt*tind;
  
  % Velocity input
  if ~useRealTrajectory
    v = constantVelocity; % m/s
  else
    v = vels(:,tind); % m/s
  end
  curDir(tind) = atan2(v(2),v(1)); % rad
  speed(tind) = sqrt(v(1)^2+v(2)^2);%/dt; % m/s
  
  x(tind) = x(tind-1)+v(1)*dt; % m
  y(tind) = y(tind-1)+v(2)*dt; % m
    
  % Dendrite frequencies are pushed up or down from the basline frequency
  % depending on the speed and head direction, with a scaling factor beta
  % that sets the spacing between the spatial grid fields.
  % Equation 6:
  dendriteFreqs = baseFreq + beta*speed(tind)*cos(curDir(tind)-dirPreferences).*(cos(curDir(tind)-dirPreferences)>0); % Hz
  
  % Advance oscillator phases
  % Radial frequency (2pi times frequency in Hz) is the time derivative of phase.
  dendritePhases = dendritePhases + dt*2*pi*dendriteFreqs; % rad
  basePhase = basePhase + dt*2*pi*baseFreq; % rad
  
  % Sum opposite oscillations
  if doSubtractPhases
    % works but not clear how this could be carried out biologically
    summedOpposites = cos(oppositeDendrites*dendritePhases'); 
  else
    summedOpposites = oppositeDendrites*cos(dendritePhases');
  end

  % Sum each dendritic oscillation separately with the baseline oscillation
  if precessionVariation=='A'
    dendritesAndBaseline = [summedOpposites];
  elseif precessionVariation=='B'
    dendritesAndBaseline = [summedOpposites; cos(basePhase)];
  elseif precessionVariation=='C'
    dendritesAndBaseline = summedOpposites + cos(basePhase);
  else
    error('No such precession option.')
  end
    
  % Rectify before product
  dendritesAndBaseline = dendritesAndBaseline.*(dendritesAndBaseline>0);

  % Final activity is the product of the oscillations.
  f = prod(dendritesAndBaseline);
  
  % threshold f
  f = f.*(f>0);
  
  % Save for later
  fhist(tind) = f;
      
  % Save firing field information
  if f>spikeThreshold
    if mod(spikeind,1000)==0
      % allocate room for next 1000 spikes
      spikeTimes(spikeind+1000) = 0;
      spikeCoords(spikeind+1000,:) = [0 0];
      spikePhases(spikeind+1000) = 0;
    end
    spikeTimes(spikeind) = t;
    spikeCoords(spikeind,:) = [x(tind) y(tind)];
    spikePhases(spikeind) = basePhase;
    spikeind = spikeind+1;
  end
  if useRealTrajectory
    xindex = round((x(tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((y(tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + double(f>spikeThreshold);
  end
  
  
  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    % We plot a rate map if using the real trajectory, otherwise
    % just the activity and trajectory with phase-coded spikes    
    figure(h);
    subplot(1,2+useRealTrajectory,1);
    plot((0:tind-1)*dt,fhist(1:tind));
    hold on;
    plot([0 tind-1]*dt,[spikeThreshold spikeThreshold],'r')
    title('Activity (blue) and threshold (red)');
    xlabel('Time (s)')
    axis square
    set(gca,'ydir','normal')
    if useRealTrajectory
      subplot(1,2+useRealTrajectory,2);
      imagesc(conv2(gaussian,spikes./(occupancy+eps)));
%       imagesc(spikes./occupancy);
      axis square
      set(gca,'ydir','normal')
      title({'Rate map',sprintf('t = %.1f s',t)})
      subplot(1,2+useRealTrajectory,3);
    else
      subplot(1,2+useRealTrajectory,2);
    end
    plot(x(1,1:10:tind),y(1,1:10:tind));
    hold on;
    if ~isempty(spikeCoords)
      cmap = jet;
      cmap = [cmap((end/2+1):end,:); cmap(1:end/2,:)];
      phaseInds = mod(spikePhases(1:spikeind-1),2*pi)*(length(cmap)-1)/2/pi;
      pointColors = cmap(ceil(phaseInds)+1,:);
      
      scatter3(spikeCoords(1:spikeind-1,1), ...
        spikeCoords(1:spikeind-1,2), ...
        zeros(size(spikeCoords(1:spikeind-1,1))), ...
        30*ones(size(spikeCoords(1:spikeind-1,1))), ...
        pointColors, 'o','filled');
    end
    axis square
%     title({'Trajectory (blue) and',...
%            'spikes (colored by theta phase',...
%            'blues before baseline peak, reds after)'})
    drawnow
  end
end

spikeTimes = spikeTimes(1:spikeind-1);
spikePhases = spikePhases(1:spikeind-1);
spikeCoords = spikeCoords(1:spikeind-1,:);

if showFigs
  figure; hist(diff(spikeTimes(1:spikeind-1)),500); title('Grid cell ISI histogram')
end
