% Blair, Gupta, and Zhang 2008's ring attractor oscillatory interference model
% eric zilli - 20110908 - v1.0
%
% Here we carry out the extension to the model that the authors suggest
% so that it can handle 2D trajectories. 
%
% See BlairEtAl2008.m for more.
%
% This code is released into the public domain. Not for use in skynet.

% if =0, uses a line segment as trajectory
% if =1, load trajectory from disk and plot activity only along that
useRealTrajectory = 1;

%% Model parameters
f0 = 7; % reference oscillator freq, Hz
nCellsPerRing = 6;
% phase offsets for each cell in the ring
phi = (1:nCellsPerRing)*2*pi/nCellsPerRing;
% inverse gain factors for each ring attractor velocity input:
ilambda = [0 1/60 1/60];
ringPreferredDirections = [0 0 pi/3];
H = [cos(ringPreferredDirections(1)) sin(ringPreferredDirections(1));
     cos(ringPreferredDirections(2)) sin(ringPreferredDirections(2));
     cos(ringPreferredDirections(3)) sin(ringPreferredDirections(3))];

%% Trajectory
% we're just using 1D trajectories here:
if useRealTrajectory
  nTrajectorySamples = 10e3;
  % compute activity at each position in trajectory
  % load trajectory from disk:
  load data/HaftingTraj_centimeters_seconds.mat;
  % package the trajectory the way we want it:
  x = pos(1:2,1:nTrajectorySamples); % cm
  dT = pos(3,2)-pos(3,1); % s
  clear pos;
else
  nTrajectorySamples = 10e3;
  x = [linspace(0,1000,nTrajectorySamples); linspace(0,0,nTrajectorySamples)]; % cm
  dT = 0.002; % s
end

%% Calculate velocity inputs
% vels is 3--by--nTrajectorySamples-1
% diff(x,1,2) is a 3--by--nTrajectorySamples-1 matrix of the x (first row)
% and y (second row) velocities. H*diff(x,1,2) projects those two rows
% onto the preferred direction vectors in H and so row 1 is now the
% velocity along direction H(1,:) and row 2 is the velocity along H(2,:), etc.
vels = H*diff(x,1,2)/dT; % cm/s

%% Oscillator frequencies are a baseline plus velocity times gain
% Note F is nTrajectorySamples-1--by--3 because vels' is nTrajectorySamples-1--by--3
% and lambda is 1-by-3 (repmat'd to nTrajectorySamples-1--by--3).
F = f0 + repmat(ilambda,size(vels,2),1).*vels'; % Hz

%% Integrate frequencies to get phase offsets for each ring:
alpha = 2*pi*cumsum(F)*dT; % nTrajectorySamples-1--by--2, rad
% Take each ring's base phase offset and make one copy for each cell in each ring
alpha = [repmat(alpha(:,1),1,nCellsPerRing) repmat(alpha(:,2),1,nCellsPerRing) repmat(alpha(:,3),1,nCellsPerRing)];

%% Calculate ring cell activities at each point in time
% This takes the base phase offset and adds the appropriate value from
% phi to find the phase offset for each individual cell in a ring.
theta = (sin(alpha + repmat(phi,nTrajectorySamples-1,3))+1)/2;

%% Grid output is the product of three ring cell activities
% We'll use the first cell in each ring
gridActivity = theta(:,1).*theta(:,nCellsPerRing+1).*theta(:,2*nCellsPerRing+1);

%% Plot the activity and trajectory
figure; plot(x(1,:),x(2,:),'b');
hold on;
plot(x(1,gridActivity>0.5),x(2,gridActivity>0.5),'r.')
title('Trajectory (blue) and spikes (red)')