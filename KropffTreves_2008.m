% Kropff and Treves 2008's grid cell model
% eric zilli - 20111010 - v1.0
% 
% In this developmental model of grid cell spatial firing, the hexagonally-
% arrayed fields of a grid cell are explained by the grid cell receiving
% input from place cells that have fields hexagonally arranged.
%
% This comes about through a slow learning process whereby a grid cell
% first develops fields at random locations (the positions where the
% initially-random synaptic inputs from place to grid are strongest), and
% then those fields drift slightly, trying to settle into a stable pattern
% (which happens to be a hexagonal grid).
%
% The drift occurs through a sort of respulsive effect that fields have on
% each other via the adaptation current of the grid cells. When a grid
% cell becomes active in its firing, the adaptation current comes on
% and decreases the cell's firing rate. Through the plasticity rule,
% place cells that fire while the grid cell's rate is decreased will be
% weakened in their influence on that cell. If the animal passes through
% a field in straight lines at all directions, this will weaken the
% input from place cells in all directions around a field.
%
% Eventually the adaptation wears off and the synaptic weights of active
% place cells onto the grid cell are once again strengthened, allowing
% fields to stably persist at that distance/time.
%
% This script was written with reference to Bailu Si's c++ code that he
% sent me. I could not implement the model based on the manuscript, and
% there were many, many details in his code that were not in the
% manuscript, but I can't be sure whether they were part of the original
% model and omitted from the manuscript or are idiosyncratic from his
% implementation, so I just point out the bits along the way that were
% missing from the manuscript so you'll know where I got each bit of
% information.
% 
% For this model it seems like we need thorough coverage of the environment
% so instead of options for constant velocity or using a real trajectory,
% I just replicate the random walk the authors used.
%
% Note: This code is not necessarily optimized for speed, but is meant as
% a transparent implementation of the model as described in the manuscript.
%
% This code is released into the public domain. Not for use in skynet.

% No "real" trajectory to use here. The model
% learns by randomly walking through a square
% environment. The Environment is 20 by 20
% "units" where 1 unit can be thought of as 5 cm.

%% Simulation parameters
dt = 1; % time step, s
simdur = 5e5; % total simulation time, s
tind = 1; % time step number for indexing
t = 0; % simulation time variable, s
x = 0; % position, m
y = 0; % position, m

%% Model parameters
% number of grid cells (thresholded saturating neurons)
NmEC = 25;
% number of place cells (inputs to grid cells)
N = 14;
N1 = N*N; % they say 200, but if the env is square, using a square number instead
% growth rate for onset of adaptation
b1 = 0.1;
% growth rate of inactivation of adaptation
b2 = b1/3;
% peak (saturating) value of grid cell activation
phiSat = 30;
% target mean activity (dynamic normalization tries to attain this mean level)
a0 = 0.1*phiSat;
% target sparsense (dynamic normalization aims for a specific sparseness too)
s0 = 0.3;
% learning rate for plasticity from place to grid cells
epsilon = 1e-5;
% place field standard deviation
placeSTD = 1; % cm
% max weight to a grid cell (not in paper?)
maxWeight = 0.08;
% for weighted running average of input and postsynaptic activities:
runningAverageWeight = 0.3;

% animal velocity (important but not specified in paper or in the authors'
% personal communication to me--just guessing based on their Figure 1A)
vel = 0.08;


% No values given for these in manuscript:
b3 = 0.1;
b4 = 0.0005;

%% Model variables
J = rand(NmEC, N1);
% Normalize weights so the sum of the squares is 1
normFactor = repmat(sum(J,2),1,N1);
J = maxWeight*J./normFactor;

% Gain g and threshold theta are dynamically set, but need starting values:
g = 1;
theta = 1;

%% History variables
% speed = zeros(1,ceil(simdur/dt));
curDir = zeros(1,ceil(simdur/dt));

%% Firing field plot variables
nSpatialBins = 60;
minx = 0; maxx = 20; % cm
miny = 0; maxy = 20; % cm
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);

spikeTimes = [];
spikeCoords = [];
spikePhases = [];

%% Place cell properties
% From Bailu's code:
cellPositions = (mod(0:(N-1),N)+0.5)*maxx/N;
[X Y] =ndgrid(cellPositions,cellPositions);
% [X Y] = ndgrid(linspace(minx,maxx,sqrt(N1)),linspace(miny,maxy,sqrt(N1)));
placeCenters = [reshape(X,1,[]); reshape(Y,1,[])];

%% Initial conditions
% Output of transfer function:
phi = zeros(1,NmEC);
% Activation variable:
r_act = ones(1,NmEC);
% Inactivation variable:
r_inact = ones(1,NmEC);
% Input vector:
r = zeros(1,N1);
% Input from place cells to grid cells:
h = zeros(1,NmEC);

meanphi = phi;
meanr = r;

%% !! Main simulation loop
fprintf('Simulation starting. Press ctrl+c to end...\n')
while t<simdur
  tind = tind+1;
  t = dt*tind;
  
  % Velocity input
  % In Bailu's script, he starts with an std of 0.12, then tries 0.742 if
  % that hits a wall, then tries 1.1 times that amount until the animal
  % is no longer walking through a wall. I'll just use a big one from the
  % start.
  randDir = curDir(tind-1) + 0.12*randn;
  while (x(tind-1)+cos(randDir))<minx || (x(tind-1)+cos(randDir))>maxx || ...
      (y(tind-1)+sin(randDir))<miny || (y(tind-1)+sin(randDir))>maxy
    randDir = curDir(tind-1) + 2*pi*randn;
  end
  
  curDir(tind) = randDir; % rad

  x(tind) = x(tind-1)+vel*cos(randDir); %v(1)*dt; % m
  y(tind) = y(tind-1)+vel*sin(randDir); %v(2)*dt; % m

  % Place cell activation at this location
  squareDists = (x(tind)-placeCenters(1,:)).^2 + (y(tind)-placeCenters(2,:)).^2;
  % Not mentioned in the paper, but phiSat also multiplies the place cell
  % activity:
  r = phiSat*exp(-squareDists/2/placeSTD^2);
  
  % Update adaptation variables
  r_act = r_act + b1*(h - r_inact - r_act);
  r_inact = r_inact + b2*(h - r_inact);

  % Place to grid input
  h = (J*r')';

  % Iteratively set the threshold and gain on this time step to produce
  % the desired mean activity and sparseness
  a = -100;
  s = -100;
  g = 1;
  theta = 0;
  iter = 0;
  % None of this is mentioned in the original paper:
  % use relative tolerance of 10% (from Bailu's code)
  while abs(a-a0)>(0.1*a0) || abs(s-s0)>(0.1*s0)
    % Transfer function (confusingly phi is the name of both a function and
    % a variable in the manuscript -- we merge them into one here.)
    thresholdedActivity = double((r_act-theta).*((r_act-theta)>0));
    % pi/2 inside atan not mentioned in paper
    phi = phiSat*2/pi*atan(pi/2*g*(thresholdedActivity));

    % average activity:
    a = sum(phi)/NmEC;
    
    % sparseness (Bailu's script floors this at 1e-6)
    if (NmEC*sum(phi.^2))>1e-6
      s = sum(phi)^2/(NmEC*sum(phi.^2));
    else
      s = 1e-6;
    end

    theta = theta + b3*(a-a0);
    % manuscript does not mention the squared mean here, I don't think:
    g = g + b4*(g+0.01)*a^2*(s-s0);

    % keep threshold below max activation variable
    if theta > max(r_act)
      theta = r_act - 0.1;
    end
    % limits on gain:
    if g > 0.5
      g = 0.5;
    elseif g < 0.01
      g = 0.01;
    end
    
    % Bailu's code gives up after 50 rounds
    iter = iter+1;
    if iter>50
      break
    end
  end
  
  % Update mean input and grid cell activities for learning
%   % Paper suggests temporal mean like:
%   meanphi = (meanphi*(tind-1)+phi)/tind;
%   meanr = (meanr*(tind-1)+r)/tind;
  % But Bailu's code does this as a sort of weighted running average:
  meanphi = meanphi*(1-runningAverageWeight) + runningAverageWeight*phi;
  meanr = meanr*(1-runningAverageWeight) + runningAverageWeight*r;
  
  % Learning
  J = J + epsilon*(phi'*r - meanphi'*meanr);
  
  % No negative weights, not mentioned in paper
  J(J<0) = 0;
  % No weights >1, not mentioned in paper
  J(J>1) = 1;
  
  % Normalize weights so the sum of the squares is 1
  normFactor = repmat(sum(J,2),1,N1);
%   normFactor = repmat(sqrt(sum(J.^2,2)),1,N1);
  J = maxWeight*J./normFactor;
  
  % Save for later
  fhist(tind) = phi(1);
end

%% Draw a few examples of weights from place cells onto grid cells
figure;
subplot(131)
imagesc(reshape(J(5,:),14,14))
axis equal
subplot(132)
imagesc(reshape(J(15,:),14,14))
axis equal
title('Weights from place cells onto 3 grid cells')
subplot(133)
imagesc(reshape(J(25,:),14,14))
axis equal

