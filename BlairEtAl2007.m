% Blair, Welday, and Zhang 2007's Moire oscillatory interference model
% eric zilli - 20110907 - v1.0
%
% This model produce the entorhinal grid pattern by assuming
% that elsewhere in the brain, grid patterns already exist but at
% much higher spatial frequencies (much smaller field spacings). The
% model then suggest that the interference pattern produced by two such
% "theta grids" can produce a larger scale grid. The theta grids in this
% implementation have a fixed relationship to position, so the activity
% needn't be updated in a manner to produce path integration: the grid cell
% activity at each location along a trajectory can all be calculated
% directly in one step.
%
% Because of this, this script does not perform path integration but rather
% plots either the total spatial pattern or the firing that
% would occur along the specified trajectory (which is just a subset of the
% overall pattern). The output is a grid cell "firing rate" that is the
% thresholded sum of the activity of the input theta grids. The manuscript
% uses two passes of a convolution filter to smooth/blur the fields, but
% they do not specify which filter was used, so we use a simple box-car
% filter.
%
% Since the model already assumes the grids exist, perhaps the primary
% contribution of the manuscript is the demonstration of how, if one
% assumes one of the theta grids is the generator of LFP theta, theta
% phase precession may or may not occur depending on how the interfering
% minigrids are chosen (i.e. the scale vs. rotation rules).
%
% Note: This code is not necessarily optimized for speed, but is meant as
% a transparent implementation of the model as described in the manuscript.
%
% This code is released into the public domain. Not for use in skynet.

% Note: to generate the BlairEtAl2007 component of my Figure 1, set:
% useRealTrajectory = 0; useLengthScalingRule = 1; alpha = 0.1;
% nSpatialBins = 1000; minx = -3; maxx = 18; miny = -3; maxy = 18; % cm


% if =0, plot whole 2D pattern.
% if =1, load trajectory from disk and plot activity only along that
useRealTrajectory = 0;

% if =0, uses rotation scaling rule (fixed phase)
% if =1, uses length scaling rule (precessing)
useLengthScalingRule = 1;

%% Model parameters
mu_M = 4; % threshold for cell to fire
c = [0; 0]; % center of one grid field; ignored right now
theta = 0; % orientation of hexagonal pattern
% for length scaling rule:
alpha = 0.05; % difference in scale of theta grids, dorsoventral values: 0.1429 >= lalpha >= 0.0667
% for rotation scaling rule:
phii = 1; % one-half difference in orientation of two theta grids, degrees
lambda = 5; % baseline spacing of theta grids, cm
omegaDirs = theta + [0 pi/3 2*pi/3]; % orientation of component gratings
omega = lambda*[cos(omegaDirs); sin(omegaDirs)]';

%% Gain function g parameters
a = 0.3;
b = -3/2;
g = @(x)(exp(a*(x-b))-1); % gain function

%% Activation function for single grid
% r is a 2-by-n vector of n coordinates
G = @(r)(g(sum(cos(omega*r))));

%% Microgrid rotation matrix for product scaling rule
R = @(phi)([cos(phi) sin(phi); -sin(phi) cos(phi)]);

%% Firing field plot variables
nSpatialBins = 1000;
minx = -100; maxx = 100; % cm
miny = -100; maxy = 100; % cm

%% Compute spatial firing
if useRealTrajectory
  % compute activity at each position in trajectory
  % load trajectory from disk:
  load data/HaftingTraj_centimeters_seconds.mat;
  % package the trajectory the way we want it:
  r = [pos(1,:); pos(2,:)];
  clear pos;
else
  % compute activity at every position in environment
  nSpatialBins = 1000;
  minx = 0; maxx = 100; % cm
  miny = 0; maxy = 100; % cm
  % This is the "environment" to use for my Figure 1:
%   nSpatialBins = 1000;
%   minx = -3; maxx = 18; % cm
%   miny = -3; maxy = 18; % cm
  xs = linspace(minx,maxx,nSpatialBins);
  ys = linspace(miny,maxy,nSpatialBins);
  [X,Y] = meshgrid(xs,ys);
  xs = reshape(X,1,[]);
  ys = reshape(Y,1,[]);
  r = [xs; ys];
  clear X Y;
end

if useLengthScalingRule
  M = G(r) + G(r*(1+alpha)) - mu_M;
else
  M = G(R(-phii)*r) + G(R(phii)*r) - mu_M;
end

% Thresholding:
M = M.*(M>0);

if useRealTrajectory
  figure; plot(r(1,:),r(2,:));
  hold on;
  plot(r(1,M>0), r(2,M>0),'r.')
  title('Spikes (red) and trajectory (blue)')
else
  M = reshape(M,nSpatialBins,[]);
  figure; imagesc(M); title('Raw firing fields');
  figure; imagesc(conv2(M,ones(6)/36)); title('Boxcar-smoothed firing fields');
  
  % No thresholding so fields smoothly fade out:
%   mM = 64*ones(size(M,1),size(M,2),3);
%   scalefact = 64/G([0; 0]);
%   mM(:,:,1) = mM(:,:,1)-scalefact*reshape(G(r*(1+alpha)),nSpatialBins,[]);
%   mM(:,:,2) = mM(:,:,2)-scalefact*reshape(G(r),nSpatialBins,[]);
%   mM(:,:,3) = mM(:,:,3)-scalefact*reshape(G(r),nSpatialBins,[]);
%   figure; image(mM/64); title('Microgrids (red and green)')
  
  % Threshold colors in plot to make fields clearer:
  mM = 64*ones(size(M,1),size(M,2),3);
  scalefact = 64;
  mM(:,:,1) = mM(:,:,1)-scalefact*(reshape(G(r*(1+alpha))>1,nSpatialBins,[]));
  mM(:,:,2) = mM(:,:,2)-scalefact*(reshape(G(r)>1,nSpatialBins,[]));
  mM(:,:,3) = mM(:,:,3)-scalefact*(reshape(G(r)>1,nSpatialBins,[]));
  figure; image(mM/64); title('Microgrids (red and green)')

end
