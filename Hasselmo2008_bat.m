% Hasselmo 2008's persistent firing oscillatory interference model
% eric zilli - 20111106 - v1.01
%
% This is the bat implementation of Hasselmo2008, accounting for
% Yartsev et al. 2011's report that grid cells in bats that lack theta
% activity in the LFP and in the ISIs/spiking autocorrelation.
%
% Because this model uses the same precession mechanism as Burgess et al.
% 2007, it naturally works with the baseline frequency set to 0.
% 
% NOTES FROM ORIGINAL SCRIPT:
% In this variation of the classic oscillatory interference model,
% the sinusoidal oscillators are thresholded into a square wave that
% is intended to represent a dense train of spikes. The activities (*)
% of three(**) square wave trains encoding position along 120-degree-spaced
% directions are multiplied together to produce the output, which also
% has a rectangular-pulse shape.
%
% (*) Actually the activities had a nonzero threshold but this was not
% mentioned in the paper.
%
% (**) The paper says 3 but he used 6 oscillators with 60 degree spacings,
% but oscillators only change frequency when the animal is going within
% 90 degrees of their preferred directions, so only 3 of the oscillators
% may change in frequency at any time. Instead of using these oscillators
% directly (see BurgessEtAl2007_precession.m), the differences between the
% phases of oppositely directed pairs of oscillators are used as the
% oscillator phases to produce the main output. This is briefly mentioned
% without much detail on page 1217, bottom left.
%
% The manuscript examines grid cell theta phase precession in the model.
% Phase precession occurs when the spike times of the grid cell occur at
% earlier and earlier phases of a reference oscillation (in vivo it is the
% theta LFP) as an animal passes through a grid field. This model does not
% use a baseline oscillation, though, so we must be careful about how we
% determine the phase of what the reference oscillation would be. Figure
% 3's caption in the manuscript appears to suggest that a fixed frequency
% baseline oscillation is used to provide the phase reference, but in
% reality that plot was generated using one of the oscillators itself
% (whose phase is the sum of two oppositely directed VCOs) as the
% reference. As simulations show, this mechanism only produces precession
% for one direction through the field, and shows procession when the field
% is crossed in the other direction (e.g. set useRealTrajectory=0 and try
% running both positive and negatives directions along the x axis by
% changing the sign of the constant velocity input).
%
% As shown by our simulations below, when the proper baseline reference is
% used, the model does not show proper precession. Instead, when the grid
% cell is active, its steady input train of boxcars causes it to spike
% constantly over the entire range of baseline oscillator phases and often
% for more than one cycle.
% 
% NB. Figure 3 is easily misinterpretted because the phase plots lack
% y axis labels. As a result, the plot appears to show 360 degrees of
% precession, when in reality it is much less than 90 degrees as careful
% comparison of lines 3 and 4 in that figure demonstrate.
%
% Note: This code is not necessarily optimized for speed, but is meant as
% a transparent implementation of the model as described in the manuscript.
%
% This code is released into the public domain. Not for use in skynet.

function [dt spikes occupancy spikeTimes] = Hasselmo2008_bat(varargin)

% Bat mode don't make no difference to this model since it doesn't
% use a baseline.
batMode = 0;
if ~isempty(varargin)
  batMode = varargin{1};
end
if batMode
  baseFreq = 0; % Hz
else
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
livePlot = 0; 20;
if showFigs
  livePlot = 0;
end

% The model doesn't use a baseline, but we can make one up for the purpose
% of defining a baseline phase.
referenceBaselinePhase = 1;

% if =0, just give constant velocity. if =1, load trajectory from disk
useRealTrajectory = 1;
constantVelocity = 1*[.5; 0*0.5]; % m/s

%% Simulation parameters
dt = .02; % time step, s
simdur = 200; % total simulation time, s
tind = 1; % time step number for indexing
t = 0; % simulation time variable, s
x = 0; % position, m
y = 0; % position, m

%% Model parameters
ncells = 1;
% Grid orientation
orientation = 0; % rad
% Directional preference of each VCO (this also sets the number of VCOs)
dirPreferences = (0:5)*pi/3;
% dirPreferences = [0 pi]; 
% Basis vectors for each head direction
H = [cos(dirPreferences+orientation)' sin(dirPreferences+orientation)'];
% This will let us add/subtract dendritic values with opposite direction preferences:
oppositeVCOs = [1 0 0 -1 0 0; 0 1 0 0 -1 0; 0 0 1 0 0 -1];
% oppositeVCOs = [1 -1];
% Scaling factor relating speed to oscillator frequencies
% Paper (e.g. Figure 3) uses P(z) = 0.0193 or P(z) = 0.0048 with units found
% elsewhere in the paper to be cycles/cm = Hz/(cm/s). 
% Pz = 0.48; % ventral; Hz/(m/s) 
Pz = 1.93; % dorsal; Hz/(m/s) 

% Threshold that determines whether the VCO output is 0 or 1
VCOThreshold = 0.5;
% Threshold that determines whether the grid cell output is a spike or not
spikeThreshold = 0;

%% History variables
speed = zeros(1,ceil(simdur/dt));
curDir = zeros(1,ceil(simdur/dt));
vhist = zeros(1,ceil(simdur/dt));
fhist = zeros(1,ceil(simdur/dt));
x = zeros(1,ceil(simdur/dt));
y = zeros(1,ceil(simdur/dt));

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
VCOPhases = zeros(1,length(dirPreferences)); % rad
basePhase = 0; % rad
VCOActivity = zeros(length(dirPreferences)/2,1);

%% Make optional figure of sheet of activity
if livePlot
  figh = figure('color','w');
  if useRealTrajectory
    set(figh,'position',[520 378 1044 420])
  end
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
  speed(tind) = sqrt(v(1)^2+v(2)^2); % m/s
  
  x(tind) = x(tind-1)+v(1)*dt; % m
  y(tind) = y(tind-1)+v(2)*dt; % m
  
  % Project the velocity onto each preferred direction vector
  h = H*v;
  
  % VCOs only change frequency vs baseline if animal is heading within
  % 90 degrees of their preferred directions.
  VCOMask = (cos(curDir(tind)-dirPreferences)>0);
  
  % VCO frequencies are pushed up or down from the baseline frequency
  % depending on the speed and head direction, with a scaling factor Pz
  % that sets the spacing between the spatial grid fields.
  VCOFreqs = baseFreq + Pz*h'.*VCOMask; % Hz
  
  % Advance oscillator phases
  % Radial frequency (2pi times frequency in Hz) is the time derivative of phase.
  VCOPhases = VCOPhases + dt*2*pi*VCOFreqs; % rad
  % Note: baseline oscillation is not used in model, but we use it to
  % determine spike phase (though, in vivo, it is not clear why theta phase
  % would be equal to this phase).
  basePhase = basePhase + dt*2*pi*baseFreq; % rad
    
  % Sum phases of oscillators with opposite direction preferences
  % NB. not clear how this could be carried out biologically
  summedPhases = oppositeVCOs*VCOPhases';
  
  % Take cosine of phases to get current oscillation activity
  summedOpposites = cos(summedPhases);

  % Sum each VCO activation
  % "~VCOActivity" is 1 for those VCOs that were not over threshold
  % on the previous step and is 0 for those VCOs that were. We can implement
  % the refractory period that was not stated in the manuscript by
  % multiplying the thresholded oscillator values on the current step by
  % that vector.
  VCOActivity = heaviside(summedOpposites-VCOThreshold);
  
  % Final activity is the product of the oscillations.
  % Note this rule has some odd features such as positive
  % activity given an even number of negative oscillator sums and
  % the baseline is included separately in each term in the product.
  f = prod(VCOActivity);
%   f = sum(VCOActivity);
  
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
    if referenceBaselinePhase
      spikePhases(spikeind) = basePhase;
    else
      spikePhases(spikeind) = summedPhases(3);
    end
    spikeind = spikeind+1;
  end
  if useRealTrajectory
    xindex = round((x(tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((y(tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + double(f>spikeThreshold);
  end
  
  
  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    if ~useRealTrajectory
      figure(figh);
      subplot(121);
      plot((0:tind-1)*dt,fhist(1:tind));
      set(gca,'ylim',[-0.1 1.1]);
      title('Activity');
      xlabel('Time (s)')
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f s',t))
      subplot(122);
      plot(x(1:tind),y(1:tind))
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
      title({'Trajectory (blue) and',...
             'spikes (colored by theta phase',...
             'blues before baseline peak, reds after)'})
      drawnow
    else
      figure(figh);
      subplot(131);
      plot((0:tind-1)*dt,fhist(1:tind));
      set(gca,'ylim',[-0.1 1.1]);
      hold on;
      title('Activity (blue)');
      xlabel('Time (s)')
      axis square
      set(gca,'ydir','normal')
      subplot(132);
      imagesc(spikes./occupancy);
      axis square
      set(gca,'ydir','normal')
      title({'Rate map',sprintf('t = %.1f s',t)})
      subplot(133);
      plot(x(1:tind),y(1:tind))
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
      title({'Trajectory (blue) and',...
             'spikes (colored by theta phase',...
             'blues before baseline peak, reds after)'})
      drawnow
    end
  end  
end

spikeTimes = spikeTimes(1:spikeind-1);
spikePhases = spikePhases(1:spikeind-1);
spikeCoords = spikeCoords(1:spikeind-1,:);

if showFigs
  figure; hist(diff(spikeTimes(1:spikeind-1)),500); title('Grid cell ISI histogram')
end

