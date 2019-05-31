% Mhatre, Gorchetchnikov, Grossberg 2010 model
% eric zilli - 20111102 - v1.0
%
% Mhatre et al. 2010 gave a self-organizing map model of grid cells.
% The input to the grid cells is provided by unbiased ring attractors.
% The firing activity of a cell in a ring is a set of parallel stripes over
% the environment. Each ring has a different orientation, but the same
% spatial scale. Each cell in a ring has a systematically offset spatial
% phase.
%
% These cells project onto a set of grid cells that form the
% self-organizing map. The grid cells compete for activity so that only one
% can have a high activity at any time. Plasticity is dependent on
% activity, so the weights from the active stripe cell inputs to that
% highly active grid cell will be increased in strength.
% 
% Mhatre et al. gave a geometric argument that a pair of stripe cells will
% be co-active most often if the angle of the two is a multiple of 60
% degrees (well, they say +/-60 degrees +/- 180 degrees which, well, what
% does that mean? I don't know. Could be a typo, I guess.) It seems to be
% true that the network results in strong synapses onto a grid cell from
% inputs close to a 60 degree angle.
%
% After the simulation the rate maps, spatial autocorrelations, and spatial
% power spectral densities are plotted. I included the latter because I
% think might be a better way to measure the global gridness, spacing,
% and orientation.
%
% I am not able to reproduce their results using spacings other than 20
% degrees between the orientations of each ring (even setting the
% stripe spacing to 30 cm and b0 to 1.125, per their Figure 8). Possibly
% some other value was also changed to allow the model to learn with more
% closely spaced orientations.
%
% I was unable to successfully implement the model based on the original
% manuscript, so this script was written in part based on source code from
% the authors themselves (and Dr. Praveen Pilly who is carrying on
% development of this model). In retrospect, part of the problem was that
% the model is dependent on the specific trajectory used (or some aspect of
% it, e.g. the mean velocity) but they did not mention this in the paper so
% I spent a lot of time using the wrong trajectory. Another issue was that
% they failed to give units for their parameters.
%
% This code is released into the public domain. Not for use in skynet.

%% Simulation variables
% they ran the same 600 s trajectory 5 times
simdur = 600; % s
trajdt = 0.020; % sampling rate for output; s (20 ms is actual sampling rate)
dt = 0.002; % s
ntrials = 5; 10;
tind=1; % index for current time step

% save history of all variables?
% faster and less memory if not
saveHistory = 0;

% if saving history, optionally downsample it
% store value if mod(tind,historySkip)==0
historySkip = 40;
histind = 1;

% if 1, don't clear W
continueSimulation = 0;

%% Trajectory
% This model is trajectory-dependent so we must use the same trajectory
% they did, which came from Sargolini et al. 2006:
% Careful, there are some NaNs in there.
load data/11207-21060501+02_t6c1

% Upsample trajectory if needed:
pos_x = interp1(t(~isnan(x1)),x1(~isnan(x1)),0:dt:t(end));
pos_y = interp1(t(~isnan(y1)),y1(~isnan(y1)),0:dt:t(end));

%% Model variables
% number of ECII cells, which one hopes will become grid cells:
nECII = 5;

% stripe angles:
stripePhis = deg2rad(-80:20:80); % as in paper, 9 unique orientations
% stripePhis = deg2rad(-80:10:80); % they say this works, but it doesn't for me
nstripeOrientations = length(stripePhis);
nstripePhases = 4;
nstripes = nstripeOrientations*nstripePhases;

% Set the orientation of each stripe input
angles = stripePhis;
stripePhis = [];
for phase=1:nstripeOrientations
  stripePhis = [stripePhis angles(phase)*ones(1,nstripes/nstripeOrientations)];
end

% They didn't bother to give us the units, but they are:
A=3; % decay rate, 1/s
B=1; % excitatory reversal potential, voltage units
D=1.5; % negative of inhibitory reversal potential, voltage units
alpha=17.5; % self-excitatory feedback gain coefficient, unitless
p=1.5; % inhibitory connection strength
lambda_g = 10; % not in paper, but multiplies a few equations -- some sort of rate parameters? 1/s
lambda_w=0.025; % learning rate, wrong value in paper, 1/s
eta=0.4; % habituative transmitter rate, 1/s
beta=0.2; % scales influence of input and feedback on rate of habituation, unitless

spacing = 20; % cm
b0 = 0.0884*spacing; % cm; value from Praveen's script

if ~continueSimulation
  % weight matrix from stripes to ECII
  W = 0.1*rand(nECII,nstripes);
%   figure; imagesc(W) % plot weights before learning
end

% ECII inhibitory matrix
P = p*(ones(nECII)-eye(nECII));

%% Firing field plot variables
nSpatialBins = 40;
minx = -50; maxx = 50; % cm
miny = -50; maxy = 50; % cm
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins,nSpatialBins,nECII);

%% Precompute stripe activities
% Project the trajectory velocity onto the preferred direction vectors
H = [cos(deg2rad(angles')) sin(deg2rad(angles'))];
projTraj = H*[pos_x; pos_y];
% Make a copy of the directional displacement for each spatial phase
projTraj = repmat(projTraj, [1 1 nstripePhases]);
% Rearrange dimensions
projTraj = permute(projTraj, [1 3 2]);
% Spatial phases:
mus = repmat(spacing*((1:nstripePhases)-1)/nstripePhases, [nstripeOrientations 1 length(pos_x)]);
% Now evaluate the stripe cells at each point along the trajectory:
strip_fun = @(mu,sigma,lambda,pos)(exp(-(min(mod(pos-mu,lambda),lambda-mod(pos-mu,lambda)).^2)/(2*sigma^2)));
g = strip_fun(mus, b0, spacing, projTraj);

%% Main simulation loop
for trial=1:ntrials
  t=0; % current time in simulation, s
  tind=1; % index for current time step
  
  if saveHistory
    x = zeros(nstripes,round(simdur/dt)+1);
    v = zeros(nECII,round(simdur/dt)+1);
    z = zeros(nECII,round(simdur/dt)+1);
    pos = zeros(2,round(simdur/dt)+1);
  else
    x = zeros(nstripes,1);
    v = zeros(nECII,1);
    z = zeros(nECII,1);
    pos = zeros(2,1);
  end
  
  spikeTimes = [];
  spikeCoords = [];
  
  fprintf('\nStarting trial %d/%d\n',trial,ntrials)
  
  tic
  while t<simdur
    t = t+dt;
    tind = tind+1;
    
    if mod(t,round(simdur)/10)<=dt
      disp(sprintf('t = %d, %g elapsed',t,toc));
    end
    
    pos = [pos_x(tind) pos_y(tind)];
    
    oldx = x(:);
    oldv = v(:);
    oldz = z(:);
    
    % update stripe activity for current location:
    x = reshape(g(:,:,tind),[],1);

    % voltage
    oldvsquared = oldv.^2;
    v = v + lambda_g*dt*(-A*oldv + (B-oldv).*((W*oldx) + alpha*oldv.^2).*oldz - (D+oldv).*(P*(oldv'.^2)'));

    % "habituative transmitter gate"
    z = z + lambda_g*dt*(eta*((1-oldz) - beta*oldz.*((W*oldx) + alpha*oldv.^2).^2));

    % learning rule: "competitive instar"
    W = W + dt*lambda_w*repmat(oldv.^2,1,nstripes).*((oldx*(2-sum(W,2)'))' - W.*repmat(sum(oldx)-oldx',nECII,1));
    % that rule is a little odd: each synapse knows the true value of all
    % the inputs to the cell
    
    % Save spiking for ratemap
    xindex = floor((pos(1)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = floor((pos(2)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex,:) = squeeze(spikes(yindex,xindex,:)) + oldvsquared;
    
    if saveHistory && mod(tind,historySkip)==0
      histind = histind+1;
      poshist(:,histind) = pos;
      xhist(:,histind) = x;
      vhist(:,histind) = v;
      zhist(:,histind) = z;
    end
    
  end
  toc
  
end

%% Plot weight matrix at end. 60 degree differences should be noticeable.
figure;
imagesc(W);
title('Stripes to ECII weight matrix');
ylabel('ECII cell');
xlabel('Stripe input (grouped by phase)');
set(gca,'xtick',(0:nstripePhases:nstripes)+nstripePhases/2);

%% Plot where voltage of each cell exceeded a threshold
if saveHistory
  figure('name','trajectory and ECII cell firing');
  for i=1:5
    subplot(3,2,i);
    plot(pos_x(1,:),pos_y(1,:),'b-',pos_x(1,vhist(i,:)>0.15),pos_y(1,vhist(i,:)>0.15),'r.');
    axis equal
  end
end

%% Plot rate map, spatial autocorrelation, and spatial power spectral density
% matrix to normalize 2d autocorrelation to yield the unbiased estimator
xcorr2NormMatrix = [1:length(occupancy) (length(occupancy)-1):-1:1];
xcorr2NormMatrix = xcorr2NormMatrix'*xcorr2NormMatrix;
for i=1:5
  figure('position',[520 378 560 225],'name',['Cell ' num2str(i)]);
  subplot(1,3,1);
  gaussWindow = fspecial('gaussian',[5 5],1);
  smoothedRateMap = conv2(spikes(:,:,i)./(occupancy+eps),gaussWindow,'same');
  imagesc(smoothedRateMap);
  title('rate map');
  axis equal;
  xlim([0 40]); ylim([0 40])
  subplot(1,3,2);
  unbiasedSpatialAutoCorr = xcorr2(smoothedRateMap)./xcorr2NormMatrix;
  imagesc(unbiasedSpatialAutoCorr);
  axis equal;
  xlim([0 79]); ylim([0 79])
  title('spatial autocorrelogram');
  subplot(1,3,3);
  spatialPSD = abs(fftshift(fft2(unbiasedSpatialAutoCorr)));
  imagesc(spatialPSD);
  title({'spatial power spectral','density amplitude'});
  axis equal;
  xlim([0 79]); ylim([0 79])
end

