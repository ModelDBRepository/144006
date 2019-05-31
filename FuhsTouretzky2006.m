% Fuhs and Touretzky 2006's continuous attractor model
% eric zilli - 20110822 - v1.0
% 
% A run of the model has three clear phases:
%   * First, the recurrent weight matrix is loaded from disk or generated.
%       This step takes a moment, but progress will be printed to the
%       console.
%   * Next, a figure appears showing the intially-random activity on the
%       sheet of grid cells. This will erode away into an
%       irregular but generally hexagonal pattern of bumps of activity
%       on the sheet of cells.
%   * A slow drift in the grid will become apparent at t = 200 ms as the
%       velocity input slowly shifts the set of active cells, carrying out
%       a path integration function. This goes on for a while and you may
%       want to ctrl+c your way out of the simulation rather than let it
%       run.
%
% When running the script keep in mind that Burak and Fiete 2006 reported
% this model does not correctly path integrate due to uncontrollable
% rotations of the activity pattern. I have found the same result.
%
% The model is a network of recurrently connected integrate-and-fire cells
% whose firing rate is given by its synaptic input.
% 
% The synaptic outputs of a cell have two components. First, each cell
% projects with an excitatory or inhibitory connection to cells in rings
% centered on it. The cell excites other cells that are very nearby,
% inhibits cells in a ring around those, excites cells in a ring a bit
% farther out, and so forth, sending alternating rings of excitation and
% inhibition. Each cell also sends an asymmetric output to a small group of
% cells very nearby the cell but offset from itself. This inhibitory output
% is used to shift the activity pattern on the grid in the opposite
% direction (see below).
% 
% The network is divided up into 2-by-2 blocks and each cell in a block is
% assigned to prefer one of N, S, E, or W. An N cell gets strong input
% (a somewhat gaussian-shaped function of direction) when the animal is
% moving north and the magnitude of the input increases with speed.
% 
% The synaptic output from any particular cell produces a bump of
% inhibition on other cells in the sheet of cells that is slightly offset
% opposite the direction of the cell's preference. That is, an N cell
% inhibits a bump of cells that is centered slightly "south" of the cell.
% "South" on the sheet of cells really means any direction, as long as it
% is consistent with all the other "south" cells and perpendicular to the
% shift directions of the "east" and "west" cells and so forth, if that
% makes any sense.
% 
% With the bump of inhibition shifted slightly backward, the cell will 
% slightly inhibit itself and cells near it, but cells in the non-inhibited
% space ahead will increase in activity. In this way, making "north" cells
% active causes the set of active cells to shift slightly along the sheet
% of cells in some particular direction. This slow shifting is what you see
% in the slowly moving fields when this script is run.
% 
% When the animal is not moving, all of the directional cells get the same
% input more or less and cancel out in their effects.
%
% The manuscript Fuhs and Touretzky (2006) also gave a developmental model
% for learning the symmetric connections needed in this model. That is
% implemented separately in FuhsTouretzky2006_development.m.
%
% IMPORTANTLY: This implementation differs from the original paper's in
% that my Phi function is just an approximation to the complicated function
% they used. An interested reader could simply replace the Phi function
% with the original (and tweak a few other parameters back to the
% manuscript values) and recover the original model. I really should do
% this at some point.
%
% This implementation also differs in that I just used a gaussian for the
% asymmetric component of the synaptic weights.
% I just wanted a quick and dirty version!
%
% This code is released into the public domain. Not for use in skynet.


% % To save the variables for FigureAttractorWeights
% % run after a simulation of constant vel = 0:
% Fu06_full_netsyn = reshape((W*f')',n,[]);
% Fu06_full_act = reshape(f',n,[]);
% ft = reshape(f,n,[]);
% fbump = zeros(size(ft));
% [v i] = max(f);
% [r c] = ind2sub([n n],i);
% fbump(r+(-4:4),c+(-4:4)) = ft(r+(-4:4),c+(-4:4)); % grab an 9x9 window around the peak as one "bump"
% Fu06_bump_netsyn = reshape((W*reshape(fbump,[],1))',n,[]);
% Fu06_bump_act = reshape(fbump,n,[]);
% f1 = zeros(size(f));
% f1(i) = 1;
% Fu06_n1_netsyn = reshape((W*f1')',n,[]);
% Fu06_n1_act = reshape(f1',n,[]);
% velInput = 1/2 + 2*0.5.*(exp(-sin((pi/2-dirs)/2).^2/0.245^2) - 1/4);
% fnorth = f + velInput;
% Fu06_north_netsyn = reshape(W*(fnorth)',n,[]);
% Fu06_north_act = reshape(fnorth',n,[]);
% save data/Fu06_WeightFigure_vars.mat Fu06_full_netsyn Fu06_full_act Fu06_bump_netsyn Fu06_bump_act Fu06_n1_netsyn Fu06_n1_act Fu06_north_netsyn Fu06_north_act
% % figure; imagesc(Fu06_bump_netsyn)
% % figure; imagesc(Fu06_bump_act)
% % figure; imagesc(Fu06_n1_netsyn)
% % figure; imagesc(Fu06_n1_act)
% % figure; imagesc(Fu06_full_netsyn)
% % figure; imagesc(Fu06_full_act)
% % figure; imagesc(Fu06_north_netsyn)
% % figure; imagesc(Fu06_north_act)

% if >0, plots the sheet of activity during the simulation on every livePlot'th step
livePlot = 100;

% if =0, just give constant velocity. if =1, load trajectory from disk
useRealTrajectory = 1;
constantVelocity = .1*[.5; 0*0.5]/1e3; % m/ms

% Weight matrix options
useCurrentW = 0; % if W exists in namespace, use that one instead of loading/generating
loadWIfPossible = 0;
saveW = 0; % save W to disk after generating

%% Cell parameters
tau = 10; % grid cell synapse time constant, ms

%% Network/Weight matrix parameters
n = 162; % originally 61x61, but an even number lets us perfectly tile with preferred directions
ncells = n*n; % total number of cells in network, 
omega = 0.67; % wave function frequency (rad/cell?)
alphaSym = 1; 0.5; % annuli amplitudes, changed from original
alphaAsym = -1.5; % offset inhibition amplitude
phi1 = 7/4; 2.55; % first cycle cutoff, changed from original
sigmaGamma = 13.375; % weight fadeout annulus
beta = 0.75; 1.5; % asymmetric offset, changed from original
% gain on velocity input (not in original model)
% higher velGain = smaller spacing but too high and things will transiently fall apart
velGain = 1000; 5000;

%% Simulation parameters
dt = 1; 0.5; % time step, ms
simdur = 120e3; % total simulation time, ms (20e3 ms = 20 s)
stabilizationTime = 200; % no-velocity time for pattern to form, ms
tind = 0; % time step number for indexing
t = 0; % simulation time variable, ms
x = 0; % position, cm
y = 0; % position, cm

%% History variables
speed = zeros(1,ceil(simdur/dt));
curDir = zeros(1,ceil(simdur/dt));
vhist = zeros(1,ceil(simdur/dt));
fhist = zeros(1,ceil(simdur/dt));

%% Firing field plot variables
watchCell = 1970; ncells/2-round(n/2); % which cell's spatial activity will be plotted?
nSpatialBins = 60;
minx = -0.90; maxx = 0.90; % m
miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);

spikeCoords = [];

%% Initial conditions
V = rand(1,ncells); % activation of each cell
f = zeros(1,ncells); % initial firing rate

%% Create 2-by-ncells preferred direction vector (radians)
dirs = [0 pi/2; pi 3*pi/2];
dirs = repmat(dirs,round(n/2),round(n/2));
dirs = reshape(dirs,1,[]);
dirVects = [cos(dirs); sin(dirs)];

%% Make x a 2-by-ncells vector of the 2D cell positions on the neural sheet
% plot as: figure; imagesc(reshape(cellDists,n,n))
x = linspace(-round(n/2),round(n/2),n);
[X,Y] = meshgrid(x,x);
x = [reshape(X,1,[]); reshape(Y,1,[])];
ell = beta/omega; % offset of center of inhibitory output
cellDists = sqrt(x(1,:).^2 + x(2,:).^2); % distance from (0,0)
% plot as: figure; imagesc(reshape(circleMask,n,[]))
circleMask = cellDists<round(n/2);
% plot as: figure; imagesc(reshape(gamma,n,[]))
gamma = exp(-(cellDists/sigmaGamma).^4);

%% Generate or load weight matrix
% Symmetric weight function (simplified version of the one Fuhs and Touretzky use)
Phi = @(x)(0.5*exp(-.25*abs(x)).*cos(2*pi*abs(x)/7)); % looks roughly the same as theirs, but faster initial decay
if ~(useCurrentW && exist('W','var'))
  fname = sprintf('data/W_FT2006_n%d_l%2g.mat',ncells,ell);
  if loadWIfPossible && exist(fname,'file')
    fprintf('Attempting to load pre-generated W...\n')
    load(fname);
    fprintf('+ Loaded pre-generated W.\n')
  else
    fprintf('Generating new W. This may take a while. Notifications at 20%% intervals...\n')

    % Define inputs weights onto each neuron i one at a time
    W = zeros(ncells);
    tic
    for i=1:ncells
      if mod(i,round(0.2*ncells))==0
        fprintf('Generating weight matrix. %d%% done.\n',round(i/ncells*100))
        toc
      end
      
      % no weights onto cells outside of main circle
      if cellDists(i)>round(n/2)
      W(i,:) = zeros(ncells,1);
        continue
      end
      
      shifts = repmat(x(:,i),1,ncells) - x - ell*dirVects;
      asymDists = sqrt(shifts(1,:).^2 + shifts(2,:).^2);
      wasym = alphaAsym*exp(-(omega*asymDists/1.5).^2)/2;
      
      shifts = repmat(x(:,i),1,ncells) - x;
      symDists = sqrt(shifts(1,:).^2 + shifts(2,:).^2);
      wsym = alphaSym*Phi(omega*symDists);

      W(i,:) = circleMask.*gamma.*(wsym + wasym);
    end

    if saveW
      save(fname,'W','-v7.3');
    end
  end
end

%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w','name','Activity of sheet of cells on brain''s surface');
  set(h,'position',[520 378 1044 420])
  drawnow
end

%% Possibly load trajectory from disk
if useRealTrajectory
  load data/HaftingTraj_centimeters_seconds.mat;
  pos(3,:) = pos(3,:)*1e3; % s to ms
  % interpolate down to simulation time step
  pos = [interp1(pos(3,:),pos(1,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(2,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(3,:),0:dt:pos(3,end))];
  pos(1:2,:) = pos(1:2,:)/100; % cm to m
  vels = [diff(pos(1,:)); diff(pos(2,:))]/dt; % m/ms
  x = pos(1,1); % m
  y = pos(2,1); % m
end
  
%% !! Main simulation loop
fprintf('Simulation starting. Press ctrl+c to end...\n')
while t<simdur
  tind = tind+1;
  t = dt*tind;

  % Velocity input
  if ~useRealTrajectory
    if t>stabilizationTime
      v = constantVelocity; % m/ms
    else
      v = [0; 0]; % m/s
    end
  else
    v = vels(:,tind); % m/ms
  end
  v = velGain*v;
  curDir(tind) = atan2(v(2),v(1)); % rad
  speed(tind) = sqrt(v(1)^2+v(2)^2);%/dt; % m/ms
  
  x = x+v(1)*dt;
  y = y+v(2)*dt;
  
  % to plot velocity input function:
  %   s = linspace(0,1,200); phi = linspace(0,2*pi,300); [S Phi] = meshgrid(s,phi); v = 1/2 + 2*S.*(exp(-sin((pi-Phi)/2).^2/0.245^2) - 1/4); figure; imagesc(phi,s,v'); set(gca,'ydir','normal'); xlabel('Direction'); ylabel('Normalized speed'); title('Velocity input')
  velInput = 1/2 + 2*speed(tind).*(exp(-sin((curDir(tind)-dirs)/2).^2/0.245^2) - 1/4);
  
  % Synaptic drive decreases with time constant tau
  V = V + dt*((W*f')' + velInput - V + 0.2*randn(size(V)))/tau;

  vhist(tind) = V(ncells/2-round(n/2));
  
  % sqrt transfer function to get firing rate:
  f = circleMask.*sqrt(V.*(V>0));

  fhist(tind) = f(ncells/2-round(n/2));
  
  if useRealTrajectory
    if f(watchCell)>0
      spikeCoords = [spikeCoords; pos(1,tind) pos(2,tind)];
    end
    xindex = round((pos(1,tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((pos(2,tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + f(watchCell);
  end
  
  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    if ~useRealTrajectory
      figure(h);
      imagesc(reshape(f,n,n));
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f ms',t))
      drawnow
    else
      figure(h);
      subplot(131);
      imagesc(reshape(f,n,n));
      axis square
      set(gca,'ydir','normal')
      subplot(132);
      imagesc(spikes./occupancy);
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f ms',t))
      subplot(133);
      plot(pos(1,1:tind),pos(2,1:tind));
      hold on;
      if ~isempty(spikeCoords)
        plot(spikeCoords(:,1),spikeCoords(:,2),'r.')
      end
      axis square
      drawnow
    end
  end
end
