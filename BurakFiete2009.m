% Burak and Fiete 2009's continuous attractor model
% eric zilli - 20110812 - v1.0
% 
% A run of the model has three clear phases:
%   * First, the recurrent weight matrix is loaded from disk or generated.
%       This step is a slooow, but progress will be printed to the console.
%       NB My synaptic matrix is thresholded at a tiny value and stored as
%       a sparse matrix. This differs from the manuscript but seems fine.
%   * Next, a figure appears showing the intially-random activity on the
%       sheet of grid cells. This will sparkle and fade away into an
%       irregular but generally hexagonal pattern of bumps of activity
%       on the sheet of cells.
%   * A slow drift in the grid will become apparent as the velocity input
%       slowly shifts the set of active cells, carrying out a path
%       integration function. This goes on for a while and you may want
%       to ctrl+c your way out of the simulation rather than let it run.
%
% The model is a network of recurrently connected abstract cells, each of
% whose firing rate is given by its synaptic input, and its synaptic
% currents have a single time constant.
% 
% Each cell sends strong inhibitory connections to certain other cells in
% the network and each cell also gets two other inputs:
%   * a global input of I=1 that causes spontaneous activity in
%       non-inhibited cells, and
%   * a time-varying "head direction" sort of input that indicates how fast
%       the animal is moving along a particular direction.
% 
% The network is divided up into 2-by-2 blocks and each cell in a block is
% assigned to prefer one of N, S, E, or W. An N cell gets strongest input
% (cosine shaped as a function of movement direction) when the animal is
% moving north and the magnitude of the input is proportional to velocity.
% 
% The synaptic output from any particular cell produces a ring of
% inhibition on other cells in the sheet of cells that is slightly offset
% in the direction of the cell's preference. That is, an N cell inhibits
% a ring of cells that is centered slightly "north" of the cell. "North"
% on the sheet of cells really means any direction, as long as it is
% consistent among all the other "north" cells and perpendicular to the
% shift directions of the "east" and "west" cells and so forth, if that
% makes any sense.
% 
% With the ring of inhibition shifted slightly forward, the cell will
% slightly inhibit itself and cells near it, but cells in the non-inhibited
% space ahead in the center of the ring will increase in activity. In this 
% way, making "north" cells active causes the set of active cells to shift 
% slightly along the sheet of cells in some particular direction. This slow
% shifting is what you see in the slowly moving fields when this script is
% run.
% 
% When the animal is not moving, all of the directional cells get the same
% input more or less and cancel out in their effects.
%
% Note: This code is not necessarily optimized for speed, but is meant as
% a transparent implementation of the model as described in the manuscript.
%
% This code is released into the public domain. Not for use in skynet.

% % To save the variables for FigureAttractorWeights:
% % run after a simulation:
% Bu09_full_netsyn = reshape((W*s')',sqrt(ncells),[]);
% Bu09_full_act = reshape(s',sqrt(ncells),[]);
% st = reshape(s,sqrt(ncells),[]);
% sbump = zeros(size(st));
% [v i] = max(s);
% [r c] = ind2sub([sqrt(ncells) sqrt(ncells)],i);
% sbump(c+(-8:8),r+(-8:8)) = st(r+(-8:8),c+(-8:8)); % grab a 17x17 window around the peak as one "bump"
% Bu09_bump_netsyn = reshape((W*reshape(sbump',[],1))',sqrt(ncells),[]);
% Bu09_bump_act = reshape(sbump',sqrt(ncells),[]);
% s1 = zeros(size(s));
% s1(i) = 1;
% Bu09_n1_netsyn = reshape((W*s1')',sqrt(ncells),[]);
% Bu09_n1_act = reshape(s1',sqrt(ncells),[]);
% snorth = (W*s')' + A.*(1+alpha*dirVects'*[0 -2]')';
% Bu09_north_netsyn = reshape((W*(snorth.*(snorth>0))')',sqrt(ncells),[]);
% Bu09_north_act = reshape((snorth.*(snorth>0))',sqrt(ncells),[]);
% figure; imagesc(Bu09_bump_netsyn)
% figure; imagesc(Bu09_bump_act)
% figure; imagesc(Bu09_n1_netsyn)
% figure; imagesc(Bu09_n1_act)
% figure; imagesc(Bu09_full_netsyn)
% figure; imagesc(Bu09_full_act)
% figure; imagesc(Bu09_north_netsyn)
% figure; imagesc(Bu09_north_act)
% save Bu09_WeightFigure_vars.mat Bu09_full_netsyn Bu09_full_act Bu09_bump_netsyn Bu09_bump_act Bu09_n1_netsyn Bu09_n1_act Bu09_north_netsyn Bu09_north_act


% if >0, plots the sheet of activity during the simulation on every livePlot'th step
livePlot = 40;

% if =0, the aperiodic network (no wrap-around weights, and a decaying
% envelope function) is used. if =1, periodic network (toroidally wrap-around
% weights, constant envelope function) is used.
usePeriodicNetwork = 1;

% Note: this simulation is quite slow (at least on my laptop) so I couldn't
% test this option very well and it may not work. Better to use a constant
% velocity to be safe.
% if =0, just give constant velocity. if =1, load trajectory from disk
useRealTrajectory = 1;
constantVelocity = 0*[0.5; 0]; % m/s


% generating W is very time consuming, wise to pre-load it if you already
% have one generated and, if you have to generate your own, to save it
% for reuse later
useCurrentW = 0; % if W exists in namespace, use that one instead of loading/generating
loadWIfPossible = 1;
saveW = 1; % save W to disk after generating (can be >1 Gb! e.g. if wSparseThresh=0)

%% Cell parameters
tau = 10; % grid cell synapse time constant, ms
if ~useRealTrajectory
  alpha = 0.10315; % input gain
else
  alpha = 50; % input gain
end

%% Network/Weight matrix parameters
ncells = 128*128; % total number of cells in network
a = 1; % if >1, uses local excitatory connections
lambda = 13; % approx the periodicity of the lattice on the neural sheet
beta = 3/(lambda^2);
% This sets how much wider the inhibitory region is than the excitatory
% region and must be large enough that it prevents the population activity
% from merging into one large bump of activity. The manuscript says 1.05*beta
% but Yoram Burak (personal communication, 2011) says that is a typo and
% the correct value is 1.02*beta, but neither of those values work here
% for some reason and I have to use something higher like 1.1*beta.
% gamma = 1.02*beta;
gamma = 1.1*beta;

% threshold for plotting a cell as having spiked
spikeThresh = 0.1;

%% Simulation parameters
dt = 0.5; % time step, ms
simdur = 100e3; % total simulation time, ms
stabilizationTime = 100; % no-velocity time for pattern to form, ms
tind = 0; % time step number for indexing
t = 0; % simulation time variable, ms

%% Initial conditions
s = rand(1,ncells); % activation of each cell

%% Firing field plot variables
watchCell = ncells/2-sqrt(ncells)/2; % which cell's spatial activity will be plotted?
nSpatialBins = 60;
minx = -0.90; maxx = 0.90; % m
miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);
spikeCoords = [];

%% Create 2-by-ncells preferred direction vector (radians)
dirs = [0 pi/2; pi 3*pi/2];
dirs = repmat(dirs,sqrt(ncells)/2,sqrt(ncells)/2);
dirs = reshape(dirs,1,[]);
dirVects = [cos(dirs); sin(dirs)];

%% Make x a 2-by-ncells vector of the 2D cell positions on the neural sheet
% plot as: figure; imagesc(reshape(cellDists,sqrt(ncells),sqrt(ncells)))
% x = linspace(-sqrt(ncells)/2,sqrt(ncells)/2,sqrt(ncells)); % what the paper suggests
x = (0:(sqrt(ncells)-1))-(sqrt(ncells)-1)/2; % what the authors actually used in their code
[X,Y] = meshgrid(x,x);
x = [reshape(X,1,[]); reshape(Y,1,[])];
cellSpacing = Y(2)-Y(1); % sets length of field shift in recurrent connections
ell = 2*cellSpacing; % offset of center of inhibitory output
cellDists = sqrt(x(1,:).^2 + x(2,:).^2); % distance from (0,0) for A below

%% Weight matrix
wSparseThresh = -1e-6; -1e-8;
if ~(useCurrentW && exist('W','var'))
  if usePeriodicNetwork
    fname = sprintf('data/W_Bu09_torus_n%d_l%1g.mat',ncells,ell);
  else
    fname = sprintf('data/W_Bu09_aperiodic_n%d_l%1g.mat',ncells,ell);
  end
  if loadWIfPossible && exist(fname,'file')
    fprintf('Attempting to load pre-generated W...\n')
    load(fname);
    fprintf('+ Loaded pre-generated W. Using a = %2g, lambda = %2g, beta = %2g, gamma = %2g, ell = %d\n',a,lambda,beta,gamma,ell)
  else
    fprintf('Generating new W. This may take a while. Notifications at 10%% intervals...\n')

    % Define inputs weights for each neuron i one at a time
    W = [];
    for i=1:ncells
      if mod(i,round(ncells/10))==0
        fprintf('Generating weight matrix. %d%% done.\n',round(i/ncells*100))
      end
      if usePeriodicNetwork
        clear squaredShiftLengths;
        % We follow Guanella et al 2007's approach to the periodic distance function
        shifts = repmat(x(:,i),1,ncells) - x - ell*dirVects;
        squaredShiftLengths(1,:) = shifts(1,:).^2 + shifts(2,:).^2;
        shifts = repmat(x(:,i),1,ncells) - x - sqrt(ncells)*[ones(1,ncells); zeros(1,ncells)] - ell*dirVects;
        squaredShiftLengths(2,:) = shifts(1,:).^2 + shifts(2,:).^2;
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[ones(1,ncells); zeros(1,ncells)] - ell*dirVects;
        squaredShiftLengths(3,:) = shifts(1,:).^2 + shifts(2,:).^2;
        shifts = repmat(x(:,i),1,ncells) - x - sqrt(ncells)*[zeros(1,ncells); ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(4,:) = shifts(1,:).^2 + shifts(2,:).^2;
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[zeros(1,ncells); ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(5,:) = shifts(1,:).^2 + shifts(2,:).^2;
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[ones(1,ncells); ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(6,:) = shifts(1,:).^2 + shifts(2,:).^2;
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[-1*ones(1,ncells); ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(7,:) = shifts(1,:).^2 + shifts(2,:).^2;
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[ones(1,ncells); -1*ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(8,:) = shifts(1,:).^2 + shifts(2,:).^2;
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[-1*ones(1,ncells); -1*ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(9,:) = shifts(1,:).^2 + shifts(2,:).^2;
        
        % Select respective least distances:
        squaredShiftLengths = min(squaredShiftLengths);
      else
        shifts = repmat(x(:,i),1,ncells) - x - ell*dirVects;
        squaredShiftLengths = shifts(1,:).^2 + shifts(2,:).^2;
      end
      temp = a*exp(-gamma*squaredShiftLengths) - exp(-beta*squaredShiftLengths);
      temp(temp>wSparseThresh) = 0;
      W = [W; sparse(temp)];
    end

    if saveW
      save(fname,'W','a','lambda','beta','gamma','ell','-v7.3');
    end
  end
end

%% Define envelope function
if usePeriodicNetwork
  % Periodic
  A = ones(size(cellDists));
else
  % Aperiodic
  % plot as: figure; imagesc(reshape(A,sqrt(ncells),sqrt(ncells)))
  R = sqrt(ncells)/2; % radius of main network in cell-position units
  a0 = sqrt(ncells)/32; % envelope fall-off rate
  dr = sqrt(ncells)/2; % diameter of non-tapered region, in cell-position units
  A = exp(-a0*(((cellDists)-R+dr)/dr).^2);
  nonTaperedInds = find(cellDists < (R-dr));
  A(nonTaperedInds) = 1;
  % zero out if >R? paper doesn't say, but it doesn't seem to matter
  % zeroInds = find(cellDists > R);
  % A(zeroInds) = 0;
end

%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w','name','Activity of sheet of cells on brain''s surface');
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
    v = [0; 0]; % m/s
  else
    if ~useRealTrajectory
      v = constantVelocity; % m/s
    else
      v = vels(:,tind); % m/s
    end
  end
  curDir(tind) = atan2(v(2),v(1)); % rad
  speed(tind) = sqrt(v(1)^2+v(2)^2);%/dt; % m/s

  
  % Feedforward input
  B = A.*(1+alpha*dirVects'*v)';
  
  % Total synaptic driving currents
  sInputs = (W*s')' + B;

  % Synaptic drive only increases if input cells are over threshold 0
  sInputs = sInputs.*(sInputs>0);

  % Synaptic drive decreases with time constant tau
  s = s + dt*(sInputs - s)/tau;

  % Save firing field information (average value of s in each spatial bin)
  if useRealTrajectory
    if s(watchCell)>spikeThresh
      spikeCoords = [spikeCoords; pos(1,tind) pos(2,tind)];
    end
    xindex = round((pos(1,tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((pos(2,tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + s(watchCell);
  end

  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    if ~useRealTrajectory
      figure(h);
      set(h,'name','Activity of sheet of cells on brain''s surface');
      imagesc(reshape(s,sqrt(ncells),sqrt(ncells)));
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f ms',t))
      drawnow
    else
      figure(h);
      subplot(131);
      imagesc(reshape(s,sqrt(ncells),sqrt(ncells)));
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
        plot(spikeCoords(:,1),spikeCoords(:,2),'r.')
      end
      axis square
      title({'Trajectory (blue)','and spikes (red)'})
      drawnow
    end
  end
end
