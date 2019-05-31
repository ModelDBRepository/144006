% Gaussier, Banquet, Sargolini, Giovannangeli, Save, Poucet 2007's grid cell model
% eric zilli - 20111110 - v1.0
%
% In this model positional information is first encoded (and path
% integrated) in the firing rate of individual neurons D with each with its
% own preferred direction theta. These cells are leak-less integrators of
% their directional velocity inputs. A cell with preferred direction theta
% starts at firing rate 0 and increases its firing rate in proportion to
% the distance the animal travels along theta (and decreases when moving
% opposite theta where activities seem to be allowed to go negative). Two
% of these cells are used here, with preferred directions at 0 and 60
% degrees.
%
% Each of these cells' activity is discretized into second populations of
% cells E (one E population for each D cell). This is carried out by
% defining a lower and upper limit of a cell's activity (the size of the
% largest area that can be path integrated), say 0 and 10, and a number of
% cells in the discretization population E, say 20. Then the active E cell
% is based directly on the firing rate of the input D cell. If the input
% activity were between 0 and 0.5 then E cell 1 could be active. With input
% between 0.5 and 1, cell 1, etc. up to inputs rates 9.5 to 10 activating E
% cell 20. This reads out the single-valued position from the input cell D
% into a population-level encoding of position E. One one cell in each E
% is active at any time, so the spatial firing pattern of an E cell is
% a single stripe parallel to the preferred direction of its input D cell
% at a distance from the animal's starting point according to the firing
% rate range it encodes from its input D cell.
%
% The stripe-like population-level encoding of position in E is then
% "folded" down into a smaller set of cells M in a manner like the modulo
% operator. When folding down cells to a set of four cells, input cell 1,
% 5, 9, 13, etc. all synaptically project to modulo cell 1. Input cells 2,
% 6, 10, 14, etc. to modulo cell 2, and so forth for the inputs of modulo
% cells 3 and 4. This compresses the stripe-like population code E into a
% modular code M.
%
% These cells M would fire like band or stripe cells used in other models,
% with spatial fields in the form of a series of thin parallel stripes
% across the environment, where the spacing depends on the relative sizes
% of E and M
%
% The activity of the grid cells is then given by assigning each grid
% cell a preferred cell in each modular population and multiplying the
% inputs so the grid cell is only active when both of its inputs are.
%
% The manuscript gives a method for self-organizing the specific spatial
% phases learned by each grid cell in appendix B.2. The model uses a very
% simple one-shot learning mechanism that adds a new grid cell every time
% the current combination of modulo population activities has not been
% previously learned. 
%
% Note that in the manuscript they seemed to discretize the velocity inputs
% themselves, which seemed to produce a large enough path integration error
% that they had to add a place-driven resetting mechanism, but I'm not
% going to discretize the input, so I won't use the resetting (to keep this
% simpler).
%
% One quirk of the paper: they seem to focus on entorhinal layer II pyramidal
% cells rather than the more populous stellate cells (perhaps just out of
% habit of referring to principal cells in cortex as pyramidals).
%
% This code is released into the public domain. Not for use in skynet.

clear all

% if >0, plots the sheet of activity during the simulation on every livePlot'th step
livePlot = 20;

% if =0, uses a line segment as trajectory
% if =1, load trajectory from disk and plot activity only along that
useRealTrajectory = 1;
constantVelocity = 1*[.5; 0*0.5]; % m/s

%% Simulation variables
simdur = 200; % s
dt = .02; % s
t = 0; % s
tind = 1;
spikeind = [1 1 1];

%% Model parameters

preferredDirections = [0 2*pi/3]';
nDirs = length(preferredDirections);
H = [cos(preferredDirections) sin(preferredDirections)];

% Minimum and maximum "activity" of D for discretizing into E:
minD = -2;
maxD = 2;
D = zeros(nDirs,1+ceil(simdur/dt));
% Number of cells used to discretize each D(i)
nE = 32;
E = zeros(nE*nDirs,1+ceil(simdur/dt));
% Number of modulo cells (spacing is nM/nO)
nM = 4;
M = zeros(nM*nDirs,1+ceil(simdur/dt));
% Number of grid cells (one for each pair of inputs from the mod populations)
% TODO should this be larger for nDirs>2?
nO = nM*nM;
O = zeros(nO,1+ceil(simdur/dt));

% tuning factor that controls which grid cells can be active
T = zeros(nO,1);
% vigilance parameter for recruiting new cells
vig = 0.9;

% O threshold for saving a spike for the spikes+trajectory plot
spikeThreshold = 0.8;


%% Synaptic weight matrices
% Weights from modulo populations to grid cells start random here
Wmo = .01*rand(nO,nDirs*nM);
% Modulo weights from discretized activity cells to modulo populations
% gotta be a better way of doing this!
Wem = [];
for i=1:nDirs
  Wem = blkdiag(Wem, repmat(eye(nM),1,nE/nM));
end

%% Firing field plot variables
watchCell = 4; % which cell's spatial activity will be plotted?
nSpatialBins = 60;
minx = -0.90; maxx = 0.90; % m
miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);
spikeTimes = [];
spikeCoords = [];
spikeCoords2 = [];
spikeCoords3 = [];


%% Trajectory
if useRealTrajectory
  load data/HaftingTraj_centimeters_seconds.mat;
  % interpolate down to simulation time step
  stopAt = find(pos(3,:)>simdur,1,'first');
  pos = [interp1(pos(3,:),pos(1,:),0:dt:pos(3,stopAt));
         interp1(pos(3,:),pos(2,:),0:dt:pos(3,stopAt));
         interp1(pos(3,:),pos(3,:),0:dt:pos(3,stopAt))];
  pos(1:2,:) = pos(1:2,:)/100; % cm to m
  vels = [diff(pos(1,:)); diff(pos(2,:))]/dt; % m/s
  x = pos(1,:); % m
  y = pos(2,:); % m
else
  'oops not done yet'
end

if livePlot
  figh = figure('position',[520 463 929 335]);
end

while t<simdur
  tind = tind+1;
  t = tind*dt;

  %% Integrate velocity inputs along preferredDirections matrix H
%   D = cumsum(H*vels)*dt; % do all at once
  D(:,tind) = D(:,tind-1) + H*vels(:,tind)*dt;

  %% Discretize activity of cells D(i) into populations E1 and E2
  for i=1:nDirs
    E(ceil(interp1([minD maxD],[1 nE],D(i,tind)))+(i-1)*nE,tind) = 1;
  end

  %% Project E onto the modulo populations
  M(:,tind) = Wem*E(:,tind);
  
  % Make a copy of the modulo inputs for each cell, for comparing to
  % its weights.
  Ii = repmat(M(:,tind)',nO,1);
  % Activation rule (equation B.4)
  thresh = @(x)(x.*double(x>0));
  O(:,tind) = T.*thresh(1-sum(abs(Wmo-Ii),2)/(nDirs*nM));
  
  %% Apply local competition:
  Oprime = O(:,tind).*(O(:,tind)>vig);
  % TODO local competition, they don't give a value for d_max... but it
  % works without local competition anyhow.

  if max(O(:,tind))<vig
    % "The learning and recruitment of a new neuron is triggered only if
    % the activity of the winner neuron is inferior to a vigilance
    % parameter vig."
    newCell = find(T==0,1,'first');
    
    if isempty(newCell)
      continue
    end
    
    T(newCell) = 1;
    Oprime(newCell,tind) = 1;
    
    % This is equation B.7. It is not clear whether only the newly
    % recruited neuron learns or whether all cells do. It sounds like
    % only the winner.
    % Winner only:
    Wmo(newCell,:) = Wmo(newCell,:) + (Ii(newCell,:)-Wmo(newCell,:)).*Oprime(newCell,tind);
%     % All cells:
%     Wmo = Wmo + (Ii-Wmo).*O(:,tind);
  end
  
  % Save firing field information
  if O(watchCell,tind)>spikeThreshold
    if mod(spikeind(1),1000)==1
      spikeTimes(spikeind(1)+1000) = 0;
      spikeCoords(spikeind(1)+1000,:) = [0 0];
    end
    spikeTimes(spikeind(1)) = t;
    spikeCoords(spikeind(1),:) = [x(tind) y(tind)];
    spikeind(1) = spikeind(1)+1;
  end
  if O(watchCell+1,tind)>spikeThreshold
    if mod(spikeind(2),1000)==1
      spikeCoords2(spikeind(2)+1000,:) = [0 0];
    end
    spikeCoords2(spikeind(2),:) = [x(tind) y(tind)];
    spikeind(2) = spikeind(2)+1;
  end
  if O(watchCell+2,tind)>spikeThreshold
    if mod(spikeind(3),1000)==1
      spikeCoords3(spikeind(3)+1000,:) = [0 0];
    end
    spikeCoords3(spikeind(3),:) = [x(tind) y(tind)];
    spikeind(3) = spikeind(3)+1;
  end
  if useRealTrajectory
    xindex = floor((x(tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = floor((y(tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + O(watchCell,tind);
  end
  
  
  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
      figure(figh);
      subplot(131);
      imagesc(O)
      title('Activity');
      xlabel('Time (s)')
      ylabel('Grid cell #')
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f s',t))
      
      subplot(132);
      imagesc(Wmo);
      axis square
      set(gca,'ydir','normal')
      title({'Weights from modulo pop.','to grid cells O'});
      xlabel('Modulo cell #')
      ylabel('Grid cell #')

      subplot(133);
      plot(pos(1,1:tind),pos(2,1:tind))
      hold on;
      if ~isempty(spikeCoords2)
        plot(spikeCoords2(1:spikeind(2)-1,1),spikeCoords2(1:spikeind(2)-1,2),'g.')
      end
      if ~isempty(spikeCoords3)
        plot(spikeCoords3(1:spikeind(3)-1,1),spikeCoords3(1:spikeind(3)-1,2),'m.')
      end
      if ~isempty(spikeCoords)
        plot(spikeCoords(1:spikeind(1)-1,1),spikeCoords(1:spikeind(1)-1,2),'r.')
      end
      axis square
      title({'Trajectory (blue) and','spikes of 3 grid cells'})
      drawnow
      
  end  

end
