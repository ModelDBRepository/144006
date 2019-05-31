% Navratilova, Giocomo, Fellous, Hasselmo, McNaughton 2011's precessing ring attractor
% eric zilli - 20110927 - v1.01
%
% This model is an unbiased ring attractor controlled by velocity. Activity
% in the ring of grid cells is shifted by sets of grid-by-direction
% conjunctive cells: here one North set that shifts the bump in one
% direction and a South set that shifts the bump in the other direction.
% All conjunctive cells receive an 8 Hz sinusoid presumably representing
% theta modulated inhibitory input and receive an input corresponding to
% velocity, but the manuscript does not specify how. We assume the standard
% in which the animal's velocity is projected onto directional vectors for
% each set of conjunctive cells and the input to the cell is proportional
% to that projection.
%
% The grid cells contain slow currents that produce phase precession
% within the network. Precession is measured relative to an 8 Hz sinusoidal
% oscillation injected into the conjunctive cells, the cells responsible
% for shifting the bump of activity among the grid cells.
%
% The grid and conjunctive cells are integrate-and-fire neurons with
% exponential AMPA and GABAA synapses and with a difference-of-exponentials
% NMDA current. The grid but not the conjunctive cells also had mAHP and
% ADP currents (NB these are based on the H current in stellate cells, says
% the manuscript, but stellates are not appropriately modeled as
% integrate-and-fire cells, despite the fact that so many people do that,
% e.g. many of the other models in this package). 
%
% This code was written in part with reference to Zaneta Navratilova's
% original MATLAB scripts due to some information that was not clear
% in the manuscript. A partial list includes
% * It was not clear whether the Gaussian function used for synapses
%     included the normalization coefficient 1/sqrt(2*pi*sigma^2).
%     (It doesn't matter, in part because she did calculated the weights
%      in a different manner and I'm using her method).
% * The 8 mV amplitude theta is modeled as 8*(1+cos) not 8*cos. Also,
%     apparently 8 mV doesn't mean the size of the resulting voltage
%     fluctuations but means the scale of the input, which is treated
%     as a voltage input rather than a current or conductance input.
% * No initial conditions was given (required to form bump in first place).
% * Navratilova's script uses a different NMDA:AMPA ratio, as well as a
%     few other parameters, than the one given in paper.
%
% Possibly with some tweaking the model as described in the manuscript
% would work, but I found it preferable to stick to her original methods.
%
% This code is released into the public domain. Not for use in skynet.

% if >0, plots the sheet of activity during the simulation on every livePlot'th step
livePlot = 20;

% if =0, just give constant velocity. if =1, load trajectory from disk
% if =2, use Navratilova's sinusoidal trajectory from her paper
useRealTrajectory = 2;
constantVelocity = 1*[.5; 0*0.5]; % m/s

%% Navratilova's sinusoidal trajectory
vel_t = 0:10:600000;
lap_freq = 0.5/60000 + 1.7/600000/60000*vel_t;
pos_cm = -150*cos(2*pi*lap_freq.*vel_t) + 150;
vel_cmps = [0 diff(fliplr(pos_cm))*1000/10];
% to scale trajectory inputs into injected currents
mVtocmps=7;
intercept=0.6348;

%% Simulation parameters
dt = 0.5; % time step, ms
simdur = 20e3; % total simulation time, ms
tind = 1; % time step number for indexing
t = 0; % simulation time variable, ms
x = 0; % position, m
y = 0; % position, m

%% Model parameters
nGrids = 100; % number of cells in the grid ring
nConjs = 100; % number of cells per conjunctive ring
dirPreferences = [0 pi]; % rad
nDirs = length(dirPreferences);
nCells = nGrids + nDirs*nConjs; % total number of cells

%% Cell parameters
tau = 10; % grid cell synapse time constant, ms
gH = 0.4; % ADP/AHP conductance  
% reversal potentials:
EH = -20; % mV
ELg = -80; % mV
ELc = -70; % mV
Eleaks = [ELg*ones(nGrids,1); ELc*ones(2*nConjs,1)];
Eexc = 0; % mV
Einh = -80; % mV
% spiking parameters:
spikeThreshold = -54; % mV
postSpikeResetV = -80; % mV
transmissionFailureChance = 0.5;

% AHP/ADP parameters:
% a mAHP time constant, b ADP peak delay, c ADP peak std. dev.
% % from Navratilova's script:
% a = 65; % ms
% b = 130; % ms
% c = 39; % ms
% dorsal:
a = 40; % ms
b = 80; % ms
c = 24; % ms
% % ventral:
% a = 115; % ms
% b = 230; % ms
% c = 69; % ms
% % sandbox:
% a = 115; % ms
% b = 330; % ms
% c = 69; % ms

%% Synaptic currents
% These functions give the current conductances for a time t after a spike
% (where t includes the axonal delay). Thus these require that t-delay >= 0.
 
AMPADecayTau = 10; % ms
AMPADelay = 2; % ms
AMPAKernel = @(t)((t>=0).*exp(-(t)/AMPADecayTau));

NMDARiseTau = 1.98; 2; % ms; Navratilova script uses 1.98 instead of 2
NMDADecayTau = 250; % ms
NMDADelay = 2; % ms
% Navratilova script multiples this by 1.0478
NMDAKernel = @(t)(1.0478*(t>=0).*(exp(-(t)/NMDADecayTau)-(exp(-(t)/NMDARiseTau))));

NMDAConductance = @(V)(1./(1 + exp(-V/16.13)*0.1/3.57)); % param V in mV

GABADecayTau = 10; % ms
GABADelay = 4; % ms
GABAKernel = @(t)((t>=0).*exp(-(t)/GABADecayTau));

%% Intrinsic currents
% These are functions of the last time a cell spiked (t>=0).
mAHPConductance = @(t)(0.5 - 0.5*exp(-t/a));
ADPConductance = @(t)(0.5*exp(-(t-b).^2/(2*c^2)));

%% Weight matrix
% Note: Figure 1B is transposed vs the normal way weight matrices are written
% Params from Table 1
% Note 2: using maxima from Navratilova's script, not the paper, to get
% more accurate weights.
ggWmax = 0.229*3*1.4; %0.0256; % maximum weight of grid->grid connections
ggSigma = 15; % std. of gaussian for grid->grid connections
gcWmax = 0.229*3.7*1.4; %0.0473; % maximum weight of grid->conjunctive connections
gcSigma = 10; % std. of gaussian for grid->conjunctive connections
cgWmax = 0.229*4.3*1.2; %0.0786; % maximum weight of conjunctive->grid connections
cgSigma = 6; % std. of gaussian for conjunctive->grid connections
cgOffset = 11; % directional shift of c->g connections in # of cells
% Using inhibitory weights from Navratilova's script intead of manuscript
baseInhib = 1.6*3*0.8*0.62/nGrids*1.4;
ggWinh = baseInhib*1.85; %0.0617; % g->g inhibitory weight
gcWinh = baseInhib/2; %0.0167; % g->c inhibitory weight
cgWinh = baseInhib/2.5; %0.0133; % c->g inhibitory weight; wrong value given in manuscript
ccWinh = baseInhib/1.5; %0.0222; % c->c inhibitory weight

NMDAtoAMPAratio = 1.7335;
% value in paper:
% NMDAtoAMPAratio = 2;

% These're all done by making a wrapped gaussian bump, swapping the first
% and second halves, then replicating that down the diagonal of a matrix.
% NB. the manuscript says these are gaussians but they are wrapped
% gaussians in Navratilova's script:
gaussian = ggWmax*(normpdf(0:nGrids-1,0.5,ggSigma)+normpdf(1:nGrids,0.5+nGrids,ggSigma));
Wgg = toeplitz(gaussian);
% The manuscript doesn't say, but Figure 1B looks like no self connections
Wgg = Wgg - Wgg(1)*eye(nGrids);

% Grid->conjunctive weights
gaussian = gcWmax*(normpdf(0:nGrids-1,0,gcSigma)+normpdf(0:nGrids-1,nGrids,gcSigma));
Wgc = toeplitz(gaussian);

% Conjunctive->grid weights
gaussian = cgWmax*(normpdf(0:nGrids-1,0,cgSigma)+normpdf(0:nGrids-1,nGrids,cgSigma));
Wc1g = toeplitz(gaussian);
Wc2g = Wc1g;
Wc1g = [Wc1g((1+cgOffset):end,:); Wc1g(1:cgOffset,:)];
Wc2g = [Wc2g((end-cgOffset+1):end,:); Wc2g(1:end-cgOffset,:)];

% Now put them all together
W = [Wgg Wc1g Wc2g; ...
     Wgc zeros(nConjs,2*nConjs); ...
     Wgc zeros(nConjs,2*nConjs)];
   
% Now make the inhibitory weights
Winh = [ggWinh*ones(nGrids,nGrids) cgWinh*ones(nGrids,2*nConjs); ...
        gcWinh*ones(2*nConjs,nGrids) ccWinh*ones(2*nConjs,2*nConjs)];

% Finally a matrix to scale the AMPA synaptic weights into the NMDA weights
% Note: only grid -> grid synapses have NMDA receptors
Wnmda = W.*[NMDAtoAMPAratio*ones(nGrids,nGrids) zeros(nGrids,2*nConjs); ...
        zeros(2*nConjs,nGrids) zeros(2*nConjs,2*nConjs)];

%% History variables
vhist = zeros(nCells,ceil(simdur/dt));

spikesHist = spalloc(nCells,ceil(simdur/dt),ceil(8*nCells*simdur/1e3));

Vs = ELc*ones(nCells,ceil(simdur/dt));

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
% initial input (slightly simplified) from Navratilova's script to 
% ignite the bumps of activity at the start of the simulation
initI = 800*(normpdf(0:nGrids-1,nGrids/2,20)+normpdf(0:nGrids-1,nGrids/2,20));
initI = [initI zeros(1,2*nConjs)]';

% last time of a synaptic event for each synapse
synapticEvents = -9000*ones(nCells);

%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w');
  if useRealTrajectory==1
    set(h,'position',[520 378 1044 420])
  end
  drawnow
end

%% Possibly load trajectory from disk
if useRealTrajectory==1
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
  if useRealTrajectory==0
    v = constantVelocity; % m/s
  elseif useRealTrajectory==1
    v = vels(:,tind); % m/s
  elseif useRealTrajectory==2
    v = [interp1(vel_t, vel_cmps, t) 0];
  end
  curDir(tind) = atan2(v(2),v(1)); % rad
  speed(tind) = sqrt(v(1)^2+v(2)^2);%/dt; % m/s
  
  x(tind) = x(tind-1)+v(1)*dt; % m
  y(tind) = y(tind-1)+v(2)*dt; % m

  %% First compute synaptic currents and ADP/mAHP
  % Max along the rows gives the first maximum, so fliplr spikesHist to
  % get last (most recent) maximum
  [vals LastSpikeInds] = max(fliplr(spikesHist(:,1:tind-1)),[],2);
  % Flip the indices back and convert to times
  LastSpikeTimes = t - dt*(LastSpikeInds-1);
  % Cells that never spiked have a max in spikesHist of 0, set those times to -1000
  LastSpikeTimes(vals==0) = -9000;
  
  % The ADP and mAHP currents have maximal conductance gH, a time-dependent
  % open conductance based on the most recent spike times of each cell, and
  % a reversal potential EH.
  ADPCurrent = gH*ADPConductance(t-LastSpikeTimes).*(Vs(:,tind-1)-EH);
  mAHPCurrent = gH*mAHPConductance(t-LastSpikeTimes).*(Vs(:,tind-1)-EH);
  % Non-grid cells have no ADP/mAHP:
  ADPCurrent(nGrids+1:end) = 0;
  mAHPCurrent(nGrids+1:end) = 0;

  % 100% chance of transmission with GABA, so just grab the most recent spikes
  if (tind-1-GABADelay/dt)>1
    [vals GABASynapticEventInds] = max(fliplr(spikesHist(:,1:tind-1-GABADelay/dt)),[],2);
    % an index of 1 means the spike was GABADelay in the past (and has just turned on)
    % an index of 2 means the spike was GABADelay+dt ago
    % An event time of 0 means the spike was GABADelay ago and the
    % resulting current is coming on right now.
    GABASynapticEventTimes = t - dt*(GABASynapticEventInds-1);
    GABASynapticEventTimes(vals==0) = -9000;
  else
    vals = zeros(nCells,1);
    GABASynapticEventTimes = -9000*ones(nCells,1);
  end
  
  % Each presynaptic spike inhibits all cells with varying synaptic weights
  % as given by Winh and a time-dependent activity by GABAKernel.
  GABASynapticConductance = (Winh*GABAKernel(t-GABASynapticEventTimes));
  GABACurrent = GABASynapticConductance.*(Vs(:,tind-1)-Einh);

  % for NMDA/AMPA, track last event time for each synapse
  if (tind-1-AMPADelay/dt)>1
    synapsesToUpdate = logical((rand(nCells,1)>transmissionFailureChance)*spikesHist(:,tind-1-AMPADelay/dt)');
    % We update a synapse if the transmission succeeded and the presynaptic
    % cell fired AMPADelay ago. Thus when the event is at the current time,
    % the synaptic current is just turning on.
    synapticEvents(synapsesToUpdate) = t;
  else
    vals = zeros(nCells,1);
  end
  
  % these two lines are the slowest part of the script. maybe faster
  % to threshold W and make it sparse?
  AMPASynapticConductance = sum(W.*AMPAKernel(t-synapticEvents),2);
  NMDASynapticConductance = sum(Wnmda.*NMDAKernel(t-synapticEvents),2);

  AMPACurrent = AMPASynapticConductance.*(Vs(:,tind-1)-Eexc);
  NMDACurrent = NMDASynapticConductance.*NMDAConductance(Vs(:,tind-1)).*(Vs(:,tind-1)-Eexc);
  
  %% Detect and reset and spikes from previous time step
  % Notice synaptic currents are calculated before the resetting
  spikes = Vs(:,tind-1)>spikeThreshold;
  spikesHist(:,tind) = spikes;
  Vs(spikes,tind-1) = postSpikeResetV;
    
  %% Set input current
  Iext = zeros(nCells,1);
  % initial input to grid cells:
  if t>200 && t<400
    Iext = Iext + initI;
  end
  Iext(nGrids+1:end) = Iext(nGrids+1:end) + 8*(1+sin(2*pi*8/1e3*t));
  
  % velocity input to conjunctive cells:
  if v(1)>0
    Iext(nGrids+(1:nConjs)) = Iext(nGrids+(1:nConjs)) + v(1)/mVtocmps + intercept;
  else
    Iext(nGrids+nConjs+(1:nConjs)) = Iext(nGrids+nConjs+(1:nConjs)) - v(1)/mVtocmps + intercept;
  end
  
  %% Update voltages
  Vs(:,tind) = Vs(:,tind-1) + dt/tau*(Eleaks-Vs(:,tind-1) - AMPACurrent - NMDACurrent - GABACurrent - ADPCurrent - mAHPCurrent + Iext);

  f = Vs(1,tind);
  % Save firing field information
  if f>spikeThreshold
    spikeTimes = [spikeTimes; t];
    spikeCoords = [spikeCoords; x(tind) y(tind)];
    spikePhases = [spikePhases; 2*pi*8/1e3*t];
  end
  if useRealTrajectory==1
    xindex = round((x(tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((y(tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + double(f>spikeThreshold);
  end
  
  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    figure(h)
    image(2.4*(80+[Vs(201:300,tind)'; Vs(101:200,tind)'; Vs(1:100,tind)']));
    title({'Population activities','(blue low, red high)',sprintf('t = %.1f ms',t)})
    text(-0.1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),1.5,'Conjunctive cells','rotation',90,'horizontalalignment','center')
    text(-0.1*diff(get(gca,'xlim'))+min(get(gca,'xlim')),3,'Grid cells','rotation',90,'horizontalalignment','center')
    text(nGrids/2,1.08*diff(get(gca,'ylim'))+min(get(gca,'ylim')),'Cell #')
  end  
end

%% Raster plot of spiking
figure; [rs cs]=find(spikesHist); plot(cs,rs,'.'); title('spike raster')
xlabel('time (ms)')
