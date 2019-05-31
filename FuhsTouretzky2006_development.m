% Fuhs and Touretzky (2006)'s developmental model
% eric zilli - 20111214 - v1.0
%
% This script demonstrates Fuhs and Touretzky (2006)'s developmental model
% of the formation of symmetric rings of output from a cell onto the other
% cells in the rectangular sheet of cells.
%
% In the model, sinusoidal gratings called "wave packets" flow across the
% sheet of cells, each beginning with a random orientationg and moving
% from one side to the other, in a direction perpendicular to the
% orientation of its stripes.
%
% The BCM-like learning rule associates a cell to the cells coactive with
% it, so as the wave packets pass any individual cell, the sinusoidal
% grating gets lightly stamped into its outputs. As gratings of random
% orientations keep passing over it, its synaptic weights contain the
% combination of gratings at every orientation, which eventually produces
% symmetric rings of synaptic output.
%
% The desired ring pattern should start to appear once 20-30 packets
% have passed (depending on the random orientations they had).
%
% This code is released into the public domain. Not for use in skynet.

% if >0, plots the sheet of activity and weights during the simulation on
% every livePlot'th step
livePlot = 100;

% Weight matrix options
useCurrentW = 0; % if W exists in namespace, use that one instead of loading/generating
loadWIfPossible = 0;
saveW = 0; % save W to disk after generating

%% Network/Weight matrix parameters
n = 62;
ncells = n*n; % total number of cells in network

% frequency of wave packets
kappa = 9*pi/31;
wavePerPacket = 3;
% tonic firing rate of cells during development (wave packets push activity
% from 0 to 2)
ftonic = 1;

%% Simulation parameters
dt = 1; % time step, ms
npackets = 1000;
simdur = 1e3; % time per packet simulation, ms
t = 0; % simulation time variable, ms
tind = 0; % time step number
x = 0; % position, cm
y = 0; % position, cm


%% Initial conditions
f = zeros(1,ncells); % initial firing rate


%% Make x a 2-by-ncells vector of the 2D cell positions on the neural sheet
% plot as: figure; imagesc(reshape(cellDists,n,n))
x = linspace(-round(n/2),round(n/2),n);
[X,Y] = meshgrid(x,x);
x = [reshape(X,1,[]); reshape(Y,1,[])];
cellDists = sqrt(x(1,:).^2 + x(2,:).^2); % distance from (0,0)

% convert to polar coordinates:
% get angles of cells around (0,0)
theta = atan2(x(2,:),x(1,:));

%% Optionally load weight matrix
fname = sprintf('data/W_FT2006_dev_n%d.mat',ncells);
if ~(useCurrentW && exist('W','var'))
  if loadWIfPossible && exist(fname,'file')
    fprintf('Attempting to load pre-generated W...\n')
    load(fname);
    fprintf('+ Loaded pre-generated W.\n')
  else
    W = zeros(ncells);
  end
end

%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w','name','Activity of sheet of cells on brain''s surface');
  set(h,'position',[520 378 1044 420])
  drawnow
end
  
%% !! Learn weight matrix (main loop)
fprintf('Development starting. Press ctrl+c to end...\n')
for packet=1:npackets
  t = 0;
  tind = -40;
  
  % random direction for this packet:
  packetDir = 2*pi*rand;
  
  % learning rate decreases from 2e3 to 1e5 "units?"
  % but they don't really tell us how. I'm guessing they mean
  % it slowed exponentially toward 1e5 rather than growing exponentially
  % toward it, and since we don't know the growth rate, I'll guess:
  tauw = 2e3 + (1e5-2e3)*(2./(1+exp(-packet/200))-1);
  
  %% Time within each packet
  while t<simdur
    t = dt*tind;
    tind = tind+1;
    
    % packet phase of each cell:
    phases = t - kappa*cellDists.*sin(theta - packetDir);

    % firing rate of each cell as a function of packet phase:
    f = ftonic + sin(phases);
    f(phases<0) = ftonic;
    f(phases>2*pi*wavePerPacket) = ftonic;
    
    % end this run if the packet has left the cells
    if all(phases>0) && all(all(f==ftonic))
      break
    end
    
    % update weight matrix
    % we use their approximation that mean(f) = ftonic
    W = W + dt/tauw*((f' - ftonic)*f - W);
    

    if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
        figure(h);
        subplot(121);
        image(32*reshape(f,n,n));
        title(sprintf('Activity of sheet of cells; t = %.1f ms, packet %d',t,packet))
        axis square
        set(gca,'ydir','normal')
        subplot(122);
        imagesc(reshape(W(:,ncells/2-round(n/2)),n,n));
        axis square
        set(gca,'ydir','normal')
        title({'Weights from center cell onto all other cells'})
        drawnow
    end
  end
end

%% Optionally save weight matrix
if saveW
  save(fname,'W','-v7.3');
end
