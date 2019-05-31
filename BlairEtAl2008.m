% Blair, Gupta, and Zhang 2008's ring attractor oscillatory interference model
% eric zilli - 20110908 - v1.0
%
% This manuscript deals with both grid and place cells, but this
% implementation focuses only on the grid cell component.
%
% (Notice this is the 1D version of the model, so the output is just a
%   series of blue bumps. Use BlairEtAl2008_2D.m for a 2D grid as output).
%
% Blair et al. 2008 overcomes one of the limitations of the Blair et al.
% 2007 model in that the existance of the grid pattern is no longer assumed
% in the inputs. Instead, this model draws inspiration from oscillatory
% interference models(*) and creates the (1D) repetitive spatial pattern from
% a pair of oscillators. (At the very end they sketch out a 2D version
% using one reference oscillator and multiple ring attractors that
% integrate velocity along directions at 60 degree angles, as in other
% 2D interference models, see BlairEtAl2008_2D.m).
%
% (*) That's sort of putting words in their mouth, but ultimately it does
% seem like ring attractor models originally came out of interference
% models if you go back to O'Keefe 1991, Burgess et al. 1992 and then
% forward to Touretzky et al. 1993 and Redish et al. 1994 who started using
% so-called "sinusoidal arrays" that look like the first steps toward
% ring attractor models, at least in the medial temporal lobe. Perhaps they
% go back much father!
%
% Like Blair et al. 2007, this model does not actually simulate the
% path integration. Instead, on each time step the ring attractor
% activities are set to the values that would occur if path integration had
% been carried out (calculated by translating the trajectory velocity
% signal to a frequency signal and integrating the frequencies to yield the
% phase of each attractor oscillator as a function of time).
%
% As a result, this script has little to do and the model amounts to
% literally multiplying together cosine gratings to make the hexagonal
% pattern.
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
% inverse gain factors for ring attractor i's velocity input:
ilambda = [0 1/60];

%% Trajectory
% we're just using 1D trajectories here:
if useRealTrajectory
  % compute activity at each position in trajectory
  % load trajectory from disk:
  load data/HaftingTraj_centimeters_seconds.mat;
  % package the trajectory the way we want it:
  % this is basically the "simulation duration" for this script:
  nTrajectorySamples = length(pos(1,:));
  x = pos(1,1:nTrajectorySamples); % cm
  dT = pos(3,2)-pos(3,1); % s
  clear pos;
else
  nTrajectorySamples = 10e3;
  x = linspace(0,1000,nTrajectorySamples); % cm
  dT = 0.002; % s
end

%% Calculate velocity inputs
vel = diff(x)/dT; % cm/s

%% Oscillator frequencies are a baseline plus velocity times gain
% Note F is nTrajectorySamples-by-2 because vel(:) is nTrajectorySamples-by-1
% and lambda is 1-by-2.
F = f0 + vel(:)*ilambda; % Hz

%% Integrate frequencies to get phase offsets for each ring:
alpha = 2*pi*cumsum(F)*dT; % nTrajectorySamples-1--by--2, rad
% Take each ring's base phase offset and make one copy for each cell in each ring
alpha = [repmat(alpha(:,1),1,nCellsPerRing) repmat(alpha(:,2),1,nCellsPerRing)];

%% Calculate ring cell activities at each point in time
% This takes the base phase offset and adds the appropriate value from
% phi to find the phase offset for each individual cell in a ring.
theta = (sin(alpha + repmat(phi,nTrajectorySamples-1,2))+1)/2;

%% Grid output is the product of two ring cell activities
% We'll use the first cell in each ring
gridActivity = theta(:,1).*theta(:,nCellsPerRing+1);

%% Plot the activity and trajectory
figure; plot(x(1:end-1),gridActivity)
title('Activity vs. position')
xlabel('1D Position (cm)')
ylabel('Activity')
