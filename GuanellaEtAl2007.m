% Guanella, Kiper, and Verschure 2007's twisted torus attractor model
% eric zilli - 20110829 - v1.01
%
% Running the model first creates a bump of activity on the sheet of cells.
% The bump will shift and wrap around the edges as the velocity input
% changes the set of active cells, carrying out path integration
% This goes on for a while and you may want to ctrl+c your way out of the
% simulation rather than let it run.
%
% The model is a network of recurrently connected abstract cells whose
% activity on each time step is solely a function of its input on that time
% step. The activity is a weighted average of the actual synaptic input to
% the cell with the synaptic input that would occur if the presynaptic
% activity were normalized to 1.
% 
% Each cell sends both excitatory and inhibitory connections to certain
% other cells in the network and cells only get this synaptic input. The
% synaptic weights are a function of distance between the cells on the
% twisted torus. The function is an exponential decay that is shifted
% downward, so each cell excites one local group of cells and inhibits
% cells in the network at higher distances.
% 
% The velocity inputs do not go to the cells in this model, but rather to
% the synaptic weights. The weight matrix changes dynamically to reflect
% the current velocity input: the offset direction of the excitatory bump
% in each cell's output weights shifts as the head direction changes. As
% the velocity changes, the size of that directional shift in the synapses
% changes proportionally. Thus each synapse is constantly changing its
% value in both magnitude and sign (excitation vs. inhibition).
% 
% For example, when the animal is moving north, the excitatory output of
% each cell shifts to the "north" direction in the sheet of cells so that
% the pattern drives a new pattern to appear, but a bit to the north. The
% faster the animal moves, the larger the shift in the "north" direction
% and so the faster the bump shifts.
% 
% When the animal is not moving, each cell's excitatory bump is centered
% on itself so the pattern does not move either.
%
% Note: Some parameters have been changed from the text to accomodate
% our real animal trajectory input. Also, the parameters of this model
% are not in terms of the time scale of the simulation, so changing the
% time step dt will change the grid spacing.
%
% Note: This code is not necessarily optimized for speed, but is meant as
% a transparent implementation of the model as described in the manuscript.
%
% This code is released into the public domain. Not for use in skynet.


% % To save the variables for FigureAttractorWeights:
% % run after a simulation:
% W = Gu07W([0 0],90);
% Gu07_full_netsyn = reshape((W*A')',Ny,[]);
% Gu07_full_act = reshape(A',Ny,[]);
% Gu07_bump_netsyn = Gu07_full_netsyn;
% Gu07_bump_act = Gu07_full_act;
% [v i] = max(A);
% A1 = zeros(size(A));
% A1(i) = 1;
% Gu07_n1_netsyn = reshape((W*A1')',Ny,[]);
% Gu07_n1_act = reshape(A1',Ny,[]);
% Wnorth = Gu07W([0 -.005],90);
% Gu07_north_netsyn = reshape((Wnorth*A')',Ny,[]);
% Gu07_north_act = reshape(A',Ny,[]);
% save Gu07_WeightFigure_vars.mat Gu07_full_netsyn Gu07_full_act Gu07_bump_netsyn Gu07_bump_act Gu07_n1_netsyn Gu07_n1_act Gu07_north_netsyn Gu07_north_act
% figure; imagesc(Gu07_bump_netsyn)
% figure; imagesc(Gu07_bump_act)
% figure; imagesc(Gu07_n1_netsyn)
% figure; imagesc(Gu07_n1_act)
% figure; imagesc(Gu07_full_netsyn)
% figure; imagesc(Gu07_full_act)
% figure; imagesc(Gu07_north_netsyn)
% figure; imagesc(Gu07_north_act)

% if >0, plots the sheet of activity during the simulation on every livePlot'th step
livePlot = 10;

% if =0, just give constant velocity. if =1, load trajectory from disk
useRealTrajectory = 1;
constantVelocity = [0.0005; 0.0000]; % m/s

%% Network/Weight matrix parameters
Nx = 10; % number of cells in x direction
Ny = 9; % number of cells in y direction
ncells = Nx*Ny; % total number of cells in network
% grid spacing is approx 1.02 - 0.48*log2(alpha), pg 236
alpha = 30; % input gain, unitless
beta = 0; % input direction bias (i.e. grid orientation), rad
sigma = 0.24; % exponential weight std. deviation
I = 0.3; % peak synaptic strength
T = 0.05; % shift so tail of exponential weights turn inhibitory
tau = 0.8; % relative weight of normalized vs. full-strength synaptic inputs

%% Simulation parameters
dt = 20; 0.5; % time step, ms
simdur = 200e3; 3*59000; % total simulation time, ms
stabilizationTime = 80; % no-velocity time for pattern to form, ms
tind = 0; % time step number for indexing
t = 0; % simulation time variable, ms
v = [0; 0]; % velocity, m/ms

%% Initial conditions
A = rand(1,ncells)/sqrt(ncells); % activation of each cell
Ahist = zeros(1,1+ceil(simdur/dt));

%% Firing field plot variables
watchCell = round(ncells/2-Ny/2); % which cell's spatial activity will be plotted?
nSpatialBins = 60;
minx = -0.90; maxx = 0.90; % m
miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);
spikeCoords = zeros(2000,2);
spikeind = 1;

% Directional input matrix
R = [cos(beta) -sin(beta); sin(beta) cos(beta)];

%% Make x a 2-by-ncells vector of the 2D cell positions on the neural sheet
x = ((1:Nx) - 0.5)/Nx;
y = sqrt(3)/2*((1:Ny) - 0.5)/Ny;
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
% column and ix the same coordinate in each row). Then ix-jx calculates
% each pairwise distance of x coordinates, and similarly iy-jy.
[jx,ix] = meshgrid(x(1,:),x(1,:));
[jy,iy] = meshgrid(x(2,:),x(2,:));
jx = reshape(jx,1,[]);
ix = reshape(ix,1,[]);
jy = reshape(jy,1,[]);
iy = reshape(iy,1,[]);
W = ones(ncells);

% This function can generate the W as needed, given a 2D velocity v and
% the number of cells
Gu07W = @(v,ncells)(I*exp(-reshape(min([(ix-jx+0+alpha*v(1)).^2 + (iy-jy+0+alpha*v(2)).^2; (ix-jx-0.5+alpha*v(1)).^2 + (iy-jy+sqrt(3)/2+alpha*v(2)).^2; (ix-jx-0.5+alpha*v(1)).^2 + (iy-jy-sqrt(3)/2+alpha*v(2)).^2; (ix-jx+0.5+alpha*v(1)).^2 + (iy-jy+sqrt(3)/2+alpha*v(2)).^2; (ix-jx+0.5+alpha*v(1)).^2 + (iy-jy-sqrt(3)/2+alpha*v(2)).^2; (ix-jx-1+alpha*v(1)).^2 + (iy-jy+0+alpha*v(2)).^2; (ix-jx+1+alpha*v(1)).^2 + (iy-jy+0+alpha*v(2)).^2]),ncells,ncells)'/sigma^2) - T);

%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w');
  drawnow
end

%% Possibly load trajectory from disk
if useRealTrajectory
  load data/HaftingTraj_centimeters_seconds.mat;
  % our time units are in ms so:
  pos(3,:) = pos(3,:)*1e3; % s to ms
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

  %% Generate new weight matrix for current velocity
  
  % to change the grid orientation, this model rotates the velocity input
  v = R*v;
  
  % Compute the pairwise distances of cells with the second cell shifted
  % in each of seven directions, then for each point pick the smallest
  % distance out of the seven shifted points.
  clear squaredPairwiseDists;
  squaredPairwiseDists = (ix-jx+0+alpha*v(1)).^2 + (iy-jy+0+alpha*v(2)).^2;
  squaredPairwiseDists(2,:) = (ix-jx-0.5+alpha*v(1)).^2 + (iy-jy+sqrt(3)/2+alpha*v(2)).^2;
  squaredPairwiseDists(3,:) = (ix-jx-0.5+alpha*v(1)).^2 + (iy-jy-sqrt(3)/2+alpha*v(2)).^2;
  squaredPairwiseDists(4,:) = (ix-jx+0.5+alpha*v(1)).^2 + (iy-jy+sqrt(3)/2+alpha*v(2)).^2;
  squaredPairwiseDists(5,:) = (ix-jx+0.5+alpha*v(1)).^2 + (iy-jy-sqrt(3)/2+alpha*v(2)).^2;
  squaredPairwiseDists(6,:) = (ix-jx-1+alpha*v(1)).^2 + (iy-jy+0+alpha*v(2)).^2;
  squaredPairwiseDists(7,:) = (ix-jx+1+alpha*v(1)).^2 + (iy-jy+0+alpha*v(2)).^2;
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

  % Save activity of one cell for nostalgia's sake
  Ahist(tind) = A(1,1);
  
  % Zero out negative activities
  A(A<0) = 0;
  
  % Save firing field information
  if useRealTrajectory
    if A(watchCell)>0
      if spikeind==size(spikeCoords,1)
        % allocate space for next 1000 spikes:
        spikeCoords(spikeind+1000,:) = [0 0];
        spikeCoords(spikeind+1,:) = [pos(1,tind) pos(2,tind)];
      else
        spikeCoords(spikeind+1,:) = [pos(1,tind) pos(2,tind)];
      end
      spikeind = spikeind+1;
    end
    xindex = round((pos(1,tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((pos(2,tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + A(watchCell);
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
      title({sprintf('t = %.1f ms',t),'Rate map'})
      subplot(133);
      plot(pos(1,1:tind),pos(2,1:tind));
      hold on;
      if ~isempty(spikeCoords)
        plot(spikeCoords(2:spikeind,1),spikeCoords(2:spikeind,2),'r.')
      end
      title({'Trajectory (blue)','and spikes (red)'})
      axis square
      drawnow
    end
  end
end
