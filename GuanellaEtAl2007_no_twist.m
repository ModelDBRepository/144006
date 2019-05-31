% Guanella, Kiper, and Verschure 2007's attractor model, non-twisted version
% eric zilli - 20110829 - v1.0
% 
% (NB. This is not a toroidal, rectangular grid model. It is a toroidal
%  *hexagonal* grid. If you want a rectangular grid, set R = eye(2);).
%
% In this version of the Guanella et al. 2007 model, we demonstrate that
% strictly the twisted torus connectivity is not necessary to produce a
% hexagonal grid. A standard torus would normally produce a rectangular
% grid, but by skewing that grid along one direction, we can change it
% into a hexagonal grid.
%
% This idea would work for both standard continuous attractor models where
% conjunctive cells are used to update the grid as well as this less
% common version of the model where the weights update dynamically.
%
% In the present case, the skew is introduced somewhat abstractly by
% multiplying the input velocities along the x and y directions by a matrix
% that leaves the x velocity alone but skews the y velocity by 60 degrees
% (pi/3 radians). To see how this is done, consider the relationship
% between the movements of a bump on a square of cells with the movements
% of the animal relative to a rhombus of four adjacent grid fields.
% 
%      Space               Cells
%    ___________         __________
%   /          /        |          |
%  /   --->   /         |   --->   |
% /__________/          |__________|
%
% As shown above, when the animal moves left in space, we want the bump
% to move left on the space of cells. Thus we will leave velocities along
% the x direction unchanged. But consider
%
%      Space               Cells
%    ___________         __________
%   /  7       /        | ^        |   (Ok, neither of those looks very
%  /  /       /         | |        |    arrow-y. Hope you get the picture!)
% /__________/          |__________|
%
% When the animal moves up and to the left in space (along a direction
% from one field to a nearest neighbor field), we want the bump of activity
% in the cells to move "up".
%
% To find the velocity transformation that accomplishes this, it is easier
% to work in the reverse direction: from sheet-of-cell directions to
% spatial directions. We will represent the transformation as matrix
% multiplication. We want a matrix X such that, first X*[1; 0] = [1; 0],
% i.e. a velocity of 1 along the x direction transforms into a velocity
% of 1 along the x direction, and second X*[0; 1] = [cos(60); sin(60)]
% Clearly X is [1 cos(60); 0 sin(6)]. To reverse the direction of this
% transformation, we just invert X. Below we call this matrix R = inv(X).
%
% Of course, that's a fairly abstract way of going about things, but it was
% the easiest way to demonstrate the mechanism.
%
% Neurally what is going on here is that there is a disagreement between
% the directional velocity that a cell receives as input and the direction
% of its asymmetric synaptic effect on the other cells. The cells receiving
% left and right velocities have the "correct" output directions, but the
% cell with an "upward" output direction are receiving velocities along
% a 60 degree angle rather than the expected 90 degree angle. 
%
% This code is released into the public domain. Not for use in skynet.

% if >0, plots the sheet of activity during the simulation on every livePlot'th step
livePlot = 100;

% if =0, just give constant velocity. if =1, load trajectory from disk
useRealTrajectory = 1;
constantVelocity = 1*[.0005; 0*0.0005]; % m/s

%% Network/Weight matrix parameters
Nx = 10; % number of cells in x direction
Ny = 10; 9; % number of cells in y direction
ncells = Nx*Ny; % total number of cells in network
% grid spacing is approx 1.02 - 0.48*log2(alpha), pg 236
alpha = 30; % input gain, unitless
beta = 0; % input direction bias (i.e. grid orientation), rad
sigma = 0.24; % exponential weight std. deviation
I = 0.3; % peak synaptic strength
T = 0.05; % shift so tail of exponential weights turn inhibitory
tau = 0.8; % relative weight of normalized vs. full-strength synaptic inputs

%% Simulation parameters
dt = 20; % time step, ms
simdur = 5*59000; % total simulation time, ms
stabilizationTime = 80; % no-velocity time for pattern to form, ms
tind = 0; % time step number for indexing
t = 0; % simulation time variable, ms
v = [0; 0]; % velocity, m/ms

%% Initial conditions
A = rand(1,ncells)/sqrt(ncells); % activation of each cell

%% Firing field plot variables
watchCell = round(ncells/2)-round(Ny/2); % which cell's spatial activity will be plotted?
nSpatialBins = 60;
minx = -0.90; maxx = 0.90; % m
miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);
spikeCoords = zeros(2000,2);
spikei = 0;

% Directional input matrix (see above for explanation)
R = double(inv([cos(0) cos(pi/3); sin(0) sin(pi/3)]));

%% Make x a 2-by-ncells vector of the 2D cell positions on the neural sheet
x = ((1:Nx) - 0.5)/Nx;
y = ((1:Ny) - 0.5)/Ny;
% y = sqrt(3)/2*((1:Ny) - 0.5)/Ny;
[X,Y] = meshgrid(x,y);
% x's first row is the x coordinates and second row the y coordinates
x = [reshape(X,1,[]); reshape(Y,1,[])];

%% Weight matrix variables
% We compute the weight matrix in one big vectorized step, so we need
% to eventually make a big matrix where entry i,j is the distance between
% cells i and j. To do this, we'll make four big matrices (that we reshape
% into vectors for later). We will calculate the distance from i to j
% along the X axis and Y axis separately, so we need the x coordinates for
% each cell i, ix, as well as the x coordinates for each cell j, jx, and
% similarly the y axes. The i and j matrices must have the coordinates
% arranged in different directions (jx has the same x coordinate in each
% column and ix the same coordinate in each row--thus j is the input
% column and i the output). Then ix-jx calculates each pairwise distance
% of x coordinates, and similarly iy-jy the y coordinate differences.
[jx,ix] = meshgrid(x(1,:),x(1,:));
[jy,iy] = meshgrid(x(2,:),x(2,:));
jx = reshape(jx,1,[]);
ix = reshape(ix,1,[]);
jy = reshape(jy,1,[]);
iy = reshape(iy,1,[]);
W = ones(ncells);

%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w');
  drawnow
end

%% Possibly load trajectory from disk
if useRealTrajectory
  load data/HaftingTraj_centimeters_seconds.mat;
  % our time units are in ms so:
  pos(3,:) = pos(3,:)*1e3;
  % interpolate down to simulation time step
  pos = [interp1(pos(3,:),pos(1,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(2,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(3,:),0:dt:pos(3,end))];
  pos(1:2,:) = pos(1:2,:)/100; % cm to m
  vels = [diff(pos(1,:)); diff(pos(2,:))]/dt; % m/s
end

%% Simulation
fprintf('Simulation starting. Press ctrl+c to end...\n')
while t<simdur
  tind = tind+1;
  t = dt*tind;
  
  % Velocity input
  if t<stabilizationTime
    if useRealTrajectory
      v = vels(:,tind); % m/s
    else
      v = [0; 0]; % m/s
    end
  else
    if useRealTrajectory
      v = vels(:,tind); % m/s
    else
      v = constantVelocity; % m/s
    end
  end
%   curDir(tind) = atan2(v(2),v(1)); % rad
%   speed(tind) = sqrt(v(1)^2+v(2)^2);%/dt; % m/s

  %% Generate new weight matrix for current velocity
  
  % to change the grid orientation, this model rotates the velocity input
  v = R*v;
  
  % Compute the pairwise distances of cells with the second cell shifted
  % in each of seven directions, then for each point pick the smallest
  % distance out of the seven shifted points.
  clear squaredPairwiseDists;
  squaredPairwiseDists = (ix-jx+0+alpha*v(1)).^2 + (iy-jy+0+alpha*v(2)).^2;
  squaredPairwiseDists(2,:) = (ix-jx-1+alpha*v(1)).^2 + (iy-jy-1+alpha*v(2)).^2;
  squaredPairwiseDists(3,:) = (ix-jx-1+alpha*v(1)).^2 + (iy-jy+0+alpha*v(2)).^2;
  squaredPairwiseDists(4,:) = (ix-jx-1+alpha*v(1)).^2 + (iy-jy+1+alpha*v(2)).^2;
  squaredPairwiseDists(5,:) = (ix-jx+0+alpha*v(1)).^2 + (iy-jy-1+alpha*v(2)).^2;
  squaredPairwiseDists(6,:) = (ix-jx+0+alpha*v(1)).^2 + (iy-jy+1+alpha*v(2)).^2;
  squaredPairwiseDists(7,:) = (ix-jx+1+alpha*v(1)).^2 + (iy-jy-1+alpha*v(2)).^2;
  squaredPairwiseDists(8,:) = (ix-jx+1+alpha*v(1)).^2 + (iy-jy+0+alpha*v(2)).^2;
  squaredPairwiseDists(9,:) = (ix-jx+1+alpha*v(1)).^2 + (iy-jy+1+alpha*v(2)).^2;
  squaredPairwiseDists = min(squaredPairwiseDists);
  squaredPairwiseDists = reshape(squaredPairwiseDists,ncells,ncells)';

  % Weights have an excitatory center that peaks at I-T and if T>0, the
  % weights are inhibitory for sufficiently high distances; specifically,
  % for distance > sigma*sqrt(-log(T/I)).
  W = I*exp(-squaredPairwiseDists/sigma^2) - T;

  % Synaptic input
  B = A*W';
  
  % Activity based on the synaptic input.
  % Notice B/sum(A) is equivalent to (A/sum(A))*W', so the second
  % term is tau times the synaptic inputs that would occur if the total
  % activity were normalized to 1. The first term is (1-tau) times
  % the full synaptic activity. tau is between 0 and 1 and weights
  % whether the input is completely normalized (tau=1) or completely
  % "raw" or unnormalized (tau=0).
  A = (1-tau)*B + tau*(B/sum(A));

  % Zero out negative activities
  A(A<0) = 0;
  
  % Save firing field information
  if useRealTrajectory
    if A(watchCell)>0.1
      spikei = spikei+1;
      if spikei==size(spikeCoords,1)
        % allocate next 1000 spikes:
        spikeCoords(spikei+1000,:) = [0 0];
        spikeCoords(spikei,:) = [pos(1,tind) pos(2,tind)];
      else
        spikeCoords(spikei,:) = [pos(1,tind) pos(2,tind)];
%         spikeCoords = [spikeCoords; pos(1,tind) pos(2,tind)];
      end
    end
    xindex = round((pos(1,tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((pos(2,tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + A(watchCell)>0.1;
  end

  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    if ~useRealTrajectory
      figure(h);
      set(h,'name','Activity of sheet of cells on brain''s surface');
      imagesc(reshape(A,Ny,Nx));
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f ms',t))
      drawnow
    else
      figure(h);
      subplot(131);
      imagesc(reshape(A,Ny,Nx));
      axis square
      title('Population activity')
      set(gca,'ydir','normal')
      subplot(132);
      imagesc(spikes./occupancy);
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f ms',t))
      subplot(133);
      plot(pos(1,1:tind),pos(2,1:tind));
      hold on;
      if spikei>0
        plot(spikeCoords(1:spikei,1),spikeCoords(1:spikei,2),'r.')
      end
      axis square
      drawnow
    end
  end
end
