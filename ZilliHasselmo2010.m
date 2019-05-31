% 2d grid simulation using izhikevich simple model neuron network
% eric zilli - 20110922 - v1.01
%
% Zilli and Hasselmo (2010) gave an implementation of an oscillatory
% interference model using more realistic neural models than the
% thresholded sinusoidal oscillators used as VCOs in previous treatments. 
%
% This is a simplification and merging of the script
% SI_simple_model_2d_grid.m (available on ModelDB in the set of
% scripts for Zilli and Hasselmo (2010; PLoS Comp Biol)) into the
% common script framework I've been using for the other models.
%
% The original script and the others in that package on ModelDB are both
% more flexible and noticeably messier.
%
% In particular, this script is set-up for two VCO configurations:
%   Type=1. Single-cell, noise-free VCOs, each simulated as one Izhikevich
%       simple model cell. (Like Figure S1, or Figure 3).
%   Type=2. Multiple, noisy, all-to-all synaptically-coupled Izhikevich simple
%       model cells in each VCO. (Like Figures 6 and 7)
%
% Of the variations described in the original paper, I chose these as
% noise and coupling were the main focus of the manuscript.
%
% The phase precession used by Burgess 2008 (blocking the output of a
% VCO when not within 90 degrees of its preferred direction) is provided
% just for fun. It wasn't in the original paper and it's a bit ad hoc,
% but it might be of interest.   
%
% This code is released into the public domain. Not for use in skynet.

% note to me: this script derives from internal version 1.03

clear all

nVCOs = 2;

% if >0, plots the activity during the simulation on every livePlot'th step
livePlot = 0; 200;

% if =0, just give constant velocity. if =1, load trajectory from disk
useRealTrajectory = 1;
constantVelocity = 1*[.5; 0*0.5]; % m/s

% this blocks the output of VCOs when the animal's current direction is not
% within 90 degrees of a VCO's preferred direction
% (set nVCOs = 6 to use, then update the postsynaptic parameters, but
% finding working ones won't be quick!)
rectifyVCOsForPrecession = 0;

%% Simulation parameters
dt = .0001; % time step, s
simdur = 3*60; % total simulation time, s
tind = 1; % time step number for indexing
t = 0; % simulation time variable, s
x = 0; % position, m
y = 0; % position, m

%% Model parameters
% The parameters for type=2 don't produce a great looking grid, but that's
% at least partially due to noise causing the grid to drif.
type = 1;
if type==1
  ncells = 1;
  useNoise = 0;
  
  load data/simple_model_RS2_FI_Jan09_n1.mat;
  baselineFreq = freqs(round(length(freqs)/2)); % Hz
  
  commonNoiseSTD = 0;
  uniqueNoiseSTD = 100*useNoise;
  g = 0;
  
  % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
  beta = 2; % Hz/(m/s)
  
  % non-gated lif params:
  tau = 40; % ms
  membraneDecay = exp(-1e3*dt/tau);
  postWeight = 0.15;
  baseWeight = 0.8;

  % gated lif params:
  basegateDur = .040; % s
  tau = 25; % (msec)
  gatedmembraneDecay = exp(-1e3*dt/tau);
  weightMult = 1.5;
  gatedpostWeight = weightMult/ncells/nVCOs;
elseif type==2
  ncells = 250;
  useNoise = 1;
  basegateDur = 1;

  load data/simple_model_RS2sn_FI_Jan09_n250.mat;
  baselineFreq = freqs(round(2+length(freqs)/2)); % Hz

  commonNoiseSTD = 0;
  uniqueNoiseSTD = 100*useNoise;
  g = 256;
  beta = 2; % Hz/(m/s)
  
  % non-gated lif params:
  tau = 6; % ms
  membraneDecay = exp(-1e3*dt/tau);
  weightMult = 0.65;
  postWeight= weightMult*1/ncells/(2+nVCOs);
  baseWeight = 2*postWeight;
  
  % gated lif params:
  basegateDur = 0.015; % s
  tau = 5; % (msec)
  gatedmembraneDecay = exp(-1e3*dt/tau);
  gatedpostWeight = 0.0012;
end

% Find input current for baseline oscillator
baseI = currents(find(freqs==baselineFreq));

% Directional preference of each VCO (this also sets the number of VCOs)
dirPreferences = (0:nVCOs-1)*pi/3;

spikeThreshold = 1;

% Grid cell (Izhikevich simple model) parameters
Cf=100; vr=-60; vt=-40; k=0.7; % parameters used for RS
a=0.03; c=-50; d=100; % neocortical pyramidal neurons
vpeak=35; % spike cutoff
b = 2;

%% Initialize values for simulation
VCOvoltage = vr*ones(ncells,nVCOs); % vr*rand(ncells,1);
VCOrecovery = 0*VCOvoltage; % initial values
bVCOvoltage = vr*ones(ncells,1); % vr*rand(ncells,1);
bVCOrecovery = 0*bVCOvoltage; % initial values

VCOSpikes = zeros(ncells,nVCOs);
baseSpikes = zeros(ncells,1);

commonNoise = 0;


%% History variables
speed = zeros(1,ceil(simdur/dt));
curDir = zeros(1,ceil(simdur/dt));
vhist = zeros(1,ceil(simdur/dt));
fhist = zeros(2,ceil(simdur/dt));
f = 0;
post = zeros(2, round(simdur/dt)+1);

x = zeros(1,ceil(simdur/dt));
y = zeros(1,ceil(simdur/dt));

clear VCOSpikeTimes;
VCOSpikeTimes{nVCOs} = []; % implicitly create this struct/class/whatever thing
BaseSpikeTimes = [];
VCOinds = zeros(1,nVCOs);
Baseind = 0;
spikeind = 0;
basegateCount = 0;

%% Firing field plot variables
nSpatialBins = 60;
minx = -0.90; maxx = 0.90; % m
miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);

spikeTimes = [];
spikeCoords = [];
spikePhases = [];

%% Initial conditions
% Oscillators will start at phase 0:
VCOPhases = zeros(1,length(dirPreferences)); % rad
basePhase = 0; % rad

%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w','name','Activity of one cell');
  if useRealTrajectory
    set(h,'position',[520 378 1044 420])
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
  x(1) = pos(1,1); % m
  y(1) = pos(2,1); % m
end

%% !! Main simulation loop
tic
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
    
  % VCO frequencies are pushed up or down from the baseline frequency
  % depending on the speed and head direction, with a scaling factor beta
  % that sets the spacing between the spatial grid fields.
  VCOFreqs = baselineFreq + beta*speed(tind)*cos(curDir(tind)-dirPreferences); % Hz
  
  % Noise shared among all the VCOs can't be corrected through coupling
  if commonNoiseSTD
    commonNoise = commonNoiseSTD*randn;
  end
  
  for vco=1:nVCOs
    oldVCOvoltage = VCOvoltage(:,vco);
    oldVCOrecovery = VCOrecovery(:,vco);

    % Set input current level for each VCO:
    % much faster if we interpolate the FI curve ourself, though we do lose
    % the ability to extrapolate, which will cause an error
    % here with one of lowind or highind being empty
    desiredFreq = VCOFreqs(vco);
    lowind = find(freqs<desiredFreq,1,'last');
    highind = find(freqs>desiredFreq,1,'first');
    proportion = (desiredFreq-freqs(lowind))/(freqs(highind)-freqs(lowind));
    I = currents(lowind) + proportion*(currents(highind)-currents(lowind));
    
    % Update voltage per Izhikevich's simple model
    if ncells>1
      VCOvoltage(:,vco) = oldVCOvoltage + 1e3*dt*(k*(oldVCOvoltage-vr).*(oldVCOvoltage-vt) - oldVCOrecovery + I + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(sum(VCOSpikes(:,vco),1)-VCOSpikes(:,vco)))/Cf;
    else
      VCOvoltage(:,vco) = oldVCOvoltage + 1e3*dt*(k*(oldVCOvoltage-vr).*(oldVCOvoltage-vt) - oldVCOrecovery + I + commonNoise + uniqueNoiseSTD*randn(ncells,1))/Cf;
    end
    VCOrecovery(:,vco) = oldVCOrecovery + 1e3*dt*a*(b*(oldVCOvoltage-vr)-oldVCOrecovery);
    
    % save and reset spikes when VCOvoltage>=vpeak
    VCOSpikes(:,vco) = VCOvoltage(:,vco)>=vpeak;
    VCOrecovery(VCOvoltage(:,vco)>=vpeak,vco) = VCOrecovery(VCOvoltage(:,vco)>=vpeak,vco)+d;

    oldVCOvoltage(VCOvoltage(:,vco)>=vpeak) = vpeak;
    VCOvoltage(VCOvoltage(:,vco)>=vpeak,vco) = c;

    % Save spike times if VCO cells spike
    if any(VCOSpikes(:,vco))
      if mod(VCOinds(vco),1000)==0
        VCOSpikeTimes{vco}(length(VCOSpikeTimes{vco})+1000) = 0;
      end
      VCOinds(vco) = VCOinds(vco)+1;
      VCOSpikeTimes{vco}(VCOinds(vco)) = t;
    end
  end
  VCOspikesum = sum(VCOSpikes,1);
  
  %% baseline VCO network
  oldbVCOvoltage = bVCOvoltage;
  oldbVCOrecovery = bVCOrecovery;
  if ncells>1
    bVCOvoltage = oldbVCOvoltage + 1e3*dt*(k*(oldbVCOvoltage-vr).*(oldbVCOvoltage-vt) - oldbVCOrecovery + baseI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(sum(baseSpikes,1)-baseSpikes))/Cf;
  else
    bVCOvoltage = oldbVCOvoltage + 1e3*dt*(k*(oldbVCOvoltage-vr).*(oldbVCOvoltage-vt) - oldbVCOrecovery + baseI + commonNoise + uniqueNoiseSTD*randn(ncells,1))/Cf;
  end
  bVCOrecovery = oldbVCOrecovery + 1e3*dt*a*(b*(oldbVCOvoltage-vr)-oldbVCOrecovery);
  
  % save and reset spikes when VCOvoltage>=vpeak
  baseSpikes = bVCOvoltage>=vpeak;
  spikesumb = sum(baseSpikes);
  if spikesumb
    if mod(Baseind,1000)==0
      BaseSpikeTimes(Baseind+1000) = 0;
    end
    Baseind = Baseind + 1;
    BaseSpikeTimes(Baseind) = t;
  end  
  bVCOrecovery(bVCOvoltage>=vpeak) = bVCOrecovery(bVCOvoltage>=vpeak)+d;
  oldbVCOvoltage(bVCOvoltage>=vpeak) = vpeak;
  bVCOvoltage(bVCOvoltage>=vpeak) = c;

  % Baseline gating (post row 2) gives much nicer grids:
  % After a baseline oscillator cell spikes, the active oscillator
  % cells are allowed to influence the postsynaptic cell for a duration
  % basegateDur s.
  if spikesumb
    basegateCount = basegateDur/dt;
  else
    basegateCount = basegateCount-1;
  end
  if basegateCount>0
    basegate = 1;
  else
    basegate = 0;
  end

  
  % Save history information:
  if mod(tind,1000)==1
    state(:,tind+1000) = [0 zeros(1,nVCOs)];
  end
  state(:,tind) = [VCOvoltage(1,:)'; bVCOvoltage(1)];
  
  if rectifyVCOsForPrecession
    VCOspikesum = VCOspikesum.*(cos(dirPreferences-curDir(tind))>0);
  end
  
  % Normal integrate-and-fire postsynaptic cell
  post(1,tind) = post(1,tind-1)*membraneDecay + postWeight*sum(VCOspikesum)+baseWeight*spikesumb;
  if post(1,tind)>1
    % Use -1e-10 as our "fired" indicator
    post(1,tind) = -1e-10;
  end
  % Integrate-and-fire postsynaptic cell where the baseline input multiplicatively
  % gates inputs from the other oscillators. That is, if the baseline input
  % occured within the past basegateDur s, active inputs are allowed to
  % affect the cell. Otherwise, they are not. Sort of cheating.
  post(2,tind) = post(2,tind-1)*gatedmembraneDecay + basegate*gatedpostWeight*sum(VCOspikesum);
  if post(2,tind)>1
    % Use -1e-10 as our "fired" indicator
    post(2,tind) = -1e-10;
  end
  
  % Save for later
  fhist(:,tind) = post(:,tind);
  
  % Save firing field information
  if fhist(1,tind)==-1e-10
    if mod(spikeind,1000)==0
      spikeTimes(spikeind+1000) = 0;
      spikeCoords(spikeind+1000,:) = [0 0];
      spikePhases(spikeind+1000) = 0;
    end
    spikeind = spikeind+1;
    spikeTimes(spikeind) = t;
    spikeCoords(spikeind,:) = [x(tind) y(tind)];
    spikePhases(spikeind) = basePhase;
  end
  if useRealTrajectory
    xindex = round((x(tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((y(tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + double(fhist(1,tind)==-1e-10);
  end
  
  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    if ~useRealTrajectory
      figure(h);
      subplot(121);
      plot(fhist(1,1:tind));
      title('Activity');
      xlabel('Time (s)')
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f s',t))
      subplot(122);
      plot(x(1:tind),y(1:tind))
      hold on;
      if ~isempty(spikeCoords)
        plot(spikeCoords(:,1),spikeCoords(:,2),'r.')
      end
      axis square
      title({'Trajectory (blue) and','spikes (red)'})
      drawnow
    else
      figure(h);
      subplot(131);
      plot((0:tind-1)*dt,fhist(1,1:tind));
      hold on;
      plot([0 tind-1]*dt,[spikeThreshold spikeThreshold],'r')
      title('Activity (blue) and threshold (red)');
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
        plot(spikeCoords(:,1),spikeCoords(:,2),'r.')
      end
      axis square
      title({'Trajectory (blue) and',...
             'spikes (red)'})
      drawnow
    end
  end  
end
toc

spikeTimes = spikeTimes(1:spikeind);
spikeCoords = spikeCoords(1:spikeind,:);
spikePhases = spikePhases(1:spikeind);
for vco=1:nVCOs
  VCOSpikeTimes{vco}(VCOSpikeTimes{vco}==0) = [];
end
BaseSpikeTimes(BaseSpikeTimes==0) = [];