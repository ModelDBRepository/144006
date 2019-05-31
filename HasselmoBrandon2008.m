% Hasselmo and Brandon 2008's cyclical persistent firing grid cell model
% eric zilli - 20111202 - v1.0
%
% Hasselmo and Brandon 2008 discussed a number of grid cell models. Here
% we implement the cyclical persistent firing model in their section 5.1.
% 
% Some important information needed to simulate the model was omitted from
% the paper, but as the model was also presented in a more complex form
% than necessary, we will implement a simpler model and first show how the
% two are related.
%
% Equation (5) in the manuscript look like this:
% 
% dhplus/dt = -omega^2*V.*posDistanceMoved.^(3/2)
% dhminus/dt = omega^2*V.*negDistanceMoved.^(3/2)
% dV/dt = hplus.*sqrt(posDistanceMoved) - hminus.*sqrt(negDistanceMoved)
%
% These are vector equations with one hplus, one hminus, and one V term
% for each head direction. posDistanceMoved and negDistanceMoved are the
% distances the animal moved along each preferred direction on the past
% time step (distance, not velocity as the paper says). Only one is nonzero
% at any time, depending on which direction the animal is moving.
%
% Before making a simplified version of this model, we first explain how
% it works (which was not described in the paper (*)). First, let's assume
% the animal has moved in a positive direction so negDistanceMoved=0 and,
% in fact, we'll just ignore hminus altogether for now. Then the model
% looks like this:
% 
% dhplus/dt = -omega^2*V.*posDistanceMoved.^(3/2)
% dV/dt = hplus.*sqrt(posDistanceMoved)%
% The equations are linear in hplus and V, so the eigenvalues of the system
% will perfectly describe its behavior. By standard methods, we find
% that the Jacobian of the system is
%
% A = [0 -omega^2*posDistanceMoved.^(3/2)]
%     [sqrt(posDistanceMoved)           0]
% 
% The eigenvalues of this matrix are clearly +/- omega*posDistanceMoved*i,
% so the system is simply an oscillator with frequency
% omega*posDistanceMoved.
%
% While standing still (posDistanceMoved=0) the frequency is 0, so the
% system also stands still (the values of hplus and V are constant). If
% the animal takes a small step, the frequency of the system momentarily
% increases in proportion to the distance moved, so the values of hplus and
% V change in proportion to distance moved. V will therefore reach a
% maximum at positions arranged in parallel bands over the environment
% (like the stripe cells of Mhatre et al 2010).
%
% With this insight, we can design a simpler and more transparent model
% that has the same eigenvalues and so behaves the same way. For example,
%
% dh/dt = -omega*V.*distMoved
% dV/dt = omega*h.*distMoved
%
% is one possibility. And indeed, this is the one we will implement.
%
% (*) I thank Michael Hasselmo for explaining the intention of
% his original equations so I could work out this simplified version.
%
% This code is released into the public domain. Not for use in skynet.

% if >0, plots activity during the simulation on every livePlot'th step
livePlot = 200;

% if =0, just give constant velocity. if =1, load trajectory from disk
useRealTrajectory = 1;
constantVelocity = 1*[.5; 0*0.5]; % m/s

%% Simulation parameters
dt = .02; % time step, s
simdur = 200; % total simulation time, s
tind = 1; % time step number for indexing
t = 0; % simulation time variable, s
x = 0; % position, m
y = 0; % position, m

%% Model parameters
omega = 300;

% Grid orientation
orientation = 0; % rad
% Directional preference of each VCO (this also sets the number of VCOs)
dirPreferences = (0:2)*pi/3;
% Basis vectors for each head direction
H = [cos(dirPreferences+orientation)' sin(dirPreferences+orientation)'];

% Threshold that determines whether the grid cell output is a spike or not
spikeThreshold = 0.25;

%% History variables
speed = zeros(1,ceil(simdur/dt));
curDir = zeros(1,ceil(simdur/dt));
vhist = zeros(1,ceil(simdur/dt));
fhist = zeros(1,ceil(simdur/dt));
x = zeros(1,ceil(simdur/dt));
y = zeros(1,ceil(simdur/dt));

%% Firing field plot variables
nSpatialBins = 60;
minx = -0.90; maxx = 0.90; % m
miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);

spikeTimes = [];
spikeCoords = [];

%% Initial conditions
% Oscillators will start at phase 0:
h = zeros(length(dirPreferences),1);
V = ones(length(dirPreferences),1);

%% Make optional figure of sheet of activity
if livePlot
  figh = figure('color','w');
  if useRealTrajectory
    set(figh,'position',[520 378 1044 420])
  end
  drawnow
end

%% Possibly load trajectory from disk
if useRealTrajectory
  load data/HaftingTraj_centimeters_seconds.mat;
  % interpolate down to simulation time step
  pos = [interp1(pos(3,:),pos(1,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(2,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(3,:),0:dt:pos(3,end))];
  pos(1:2,:) = pos(1:2,:)/100; % cm to m
  vels = [diff(pos(1,:)); diff(pos(2,:))]/dt; % m/s
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
    v = constantVelocity; % m/s
  else
    v = vels(:,tind); % m/s
  end
  curDir(tind) = atan2(v(2),v(1)); % rad
  speed(tind) = sqrt(v(1)^2+v(2)^2); % m/s
  
  x(tind) = x(tind-1)+v(1)*dt; % m
  y(tind) = y(tind-1)+v(2)*dt; % m
  
  % Project the velocity onto each preferred direction vector and multiply
  % by the time step to get the distance moved along each preferred dir.
  distMoved = H*v*dt;
  
  % NB We're using foward Euler here, but in general forward Euler is not
  % save to use with linear systems like this. It seems like forward Euler
  % is worst when the real part of the eigenvalues are small but nonzero,
  % but our real parts are 0 here so it seems ok. Just be careful though!
  % The midpoint method is better for linear systems (feel free to email
  % me for code solving a linear system with the midpoint method if you're
  % interested).
  h = h + dt*-omega*V.*distMoved;
  V = V + dt*omega*h.*distMoved;
  
  g = prod(V);
    
  % Save for later
  ghist(tind) = g;
      
  % Save firing field information
  if g>spikeThreshold
    spikeTimes = [spikeTimes; t];
    spikeCoords = [spikeCoords; x(tind) y(tind)];
  end
  if useRealTrajectory
    xindex = round((x(tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((y(tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + double(g>spikeThreshold);
  end
  
  
  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    if ~useRealTrajectory
      figure(figh);
      subplot(121);
      plot((0:tind-1)*dt,ghist(1:tind));
      set(gca,'ylim',[-0.1 1.1]);
      title('Activity');
      xlabel('Time (s)')
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f s',t))
      subplot(122);
      plot(x(1:tind),y(1:tind))
      hold on;
      if ~isempty(spikeCoords)
        plot(spikeCoords(:,1),spikeCoords(:,2),'r.')
      end
      axis square
      title({'Trajectory (blue) and','spikes (red)'})
      drawnow
    else
      figure(figh);
      subplot(131);
      plot((0:tind-1)*dt,ghist(1:tind));
      set(gca,'ylim',[-0.1 1.1]);
      hold on;
      title('Activity (blue)');
      xlabel('Time (s)')
      axis square
      set(gca,'ydir','normal')
      subplot(132);
      imagesc(spikes./occupancy);
      axis square
      set(gca,'ydir','normal')
      title({'Rate map',sprintf('t = %.1f s',t)})
      subplot(133);
      plot(x(1:tind),y(1:tind))
      hold on;
      if ~isempty(spikeCoords)
        plot(spikeCoords(:,1),spikeCoords(:,2),'r.')
      end
      axis square
      title({'Trajectory (blue) and','spikes (red)'})
      drawnow
    end
  end  
end

