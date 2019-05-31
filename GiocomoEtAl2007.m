% Giocomo, Zilli, Fransen, and Hasselmo 2007's temporal interference model
% eric zilli - 20110824 - v1.01
% 
% This model was a slight variation on the model originally described in
% the Burgess, Barry, Jeffery, O'Keefe poster "A Grid & Place Cell Model
% of Path Integration Utilizing Phase Precession Versus Theta."
%
% The only change to the model itself was in the way the frequencies of the
% active oscillators are set (which changes the way grid spacing relates to
% the parameters of the model). Let s(t) and phi(t) be the animal's speed
% and direction at time t, phi the preferred direction of an oscillator,
% w_active its frequency, w_baseline the baseline frequency, and beta a
% scaling factor relating velocities to frequency changes. Originally the
% active oscillator frequency was
% w_active = w_baseline + beta*s(t)*cos(phi(t) - phi)
%
% In Giocomo et al. 2007, Mike modified it to be
% w_active = w_baseline + w_baseline*beta*s(t)*cos(phi(t) - phi)
%
% Grid field spacing in these models depends on how the frequency
% difference w_active - w_baseline is related to the animal's directional
% speed s(t)*cos(phi(t) - phi). The change to these equations therefore
% results in field spacing changing from depending on beta alone in the
% original equation to w_baseline*beta in the second equation.
%
% Notice that if we had two separate simulations, we could set beta in the
% simulation of the first equation equal to the value of w_baseline*beta in
% the second equation, and the two equations would always give the same
% result. 
%
% Thus the models are identical except the field spacing in Giocomo et al.
% 2007 depends on the baseline frequency. 
%
% This change was made because Lisa Giocomo's data showed a gradient of
% frequencies along the dorsoventral axis of entorhinal cortex from say 8
% Hz to 4 Hz. First, all along the DV axis, theta frequency is around 8 Hz
% so the difference between the frequencies of the cells near threshold and
% theta increases in the ventral direction. In this model, a larger
% frequency difference produces smaller field spacing, but ventral grid
% cells have larger field spacings, so the theta oscillation could not play
% the role of the baseline frequency in this model (e.g. this model will
% not precess relative to theta).
%
% The assumption then had to be that the 8 to 4 Hz gradient was the
% baseline frequency. But in the original model, the field spacing is not
% related to baseline frequency, so the model was modified to be consistent
% with this data. By making the difference in frequency depend on the
% baseline frequency, a lower baseline frequency directly caused larger
% field spacing assuming beta remained unchanged.
%
% Burgess et al. 2007 pointed out that Lisa's data might instead correspond
% to a gradient of beta values, making the baseline dependency introduced
% here unnecessary.
%
% Note: This code is not necessarily optimized for speed, but is meant as
% a transparent implementation of the model as described in the manuscript.
%
% This code is released into the public domain. Not for use in skynet.

% if >0, plots the sheet of activity during the simulation on every livePlot'th step
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
ncells = 1;
% Basline maintains a fixed frequency
baseFreq = 6.42; % dorsal, Hz
% baseFreq = 4.23; % ventral, Hz
% Directional preference of each dendrite (this also sets the number of dendrites)
dirPreferences = [0 2*pi/3 4*pi/3];
% Scaling factor relating speed to oscillator frequencies
beta = 0.385; % Hz/(m/s) 
spikeThreshold = 1.8;


%% History variables
speed = zeros(1,ceil(simdur/dt));
curDir = zeros(1,ceil(simdur/dt));
vhist = zeros(1,ceil(simdur/dt));
fhist = zeros(1,ceil(simdur/dt));

%% Firing field plot variables
nSpatialBins = 60;
minx = -0.90; maxx = 0.90; % m
miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);

spikeTimes = [];
spikeCoords = [];
spikePhases = [];

%% Initial conditions
% Oscillators will start at phase 0:
dendritePhases = zeros(1,length(dirPreferences)); % rad
basePhase = 0; % rad

%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w','name','Activity of one cell');
  if useRealTrajectory
    set(h,'position',[520 378 1044 420])
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
  speed(tind) = sqrt(v(1)^2+v(2)^2);%/dt; % m/s
  
  x(tind) = x(tind-1)+v(1)*dt; % m
  y(tind) = y(tind-1)+v(2)*dt; % m
    
  % Dendrite frequencies are pushed up or down from the basline frequency
  % depending on the speed and head direction, with a scaling factor
  % baseFreq*beta that sets the spacing between the spatial grid fields.
  dendriteFreqs = baseFreq + baseFreq*beta*speed(tind)*cos(curDir(tind)-dirPreferences); % Hz

  % Advance oscillator phases
  % Radial frequency (2pi times frequency in Hz) is the time derivative of phase.
  dendritePhases = dendritePhases + dt*2*pi*dendriteFreqs; % rad
  basePhase = basePhase + dt*2*pi*baseFreq; % rad
  
  % Sum each dendritic oscillation separately with the baseline oscillation
  dendritePlusBaseline = cos(dendritePhases) + cos(basePhase);
    
  % Final activity is the product of the oscillations.
  % Note this rule has some odd features such as positive
  % activity given an even number of negative oscillator sums and
  % the baseline is included separately in each term in the product.
  f = prod(dendritePlusBaseline);
  
  % threshold threshold f
  f = f.*(f>0);
  
  % Save for later
  fhist(tind) = f;
  
  % Save firing field information
  if f>spikeThreshold
    spikeTimes = [spikeTimes; t];
    spikeCoords = [spikeCoords; x(tind) y(tind)];
    spikePhases = [spikePhases; basePhase];
  end
  if useRealTrajectory
    xindex = round((x(tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((y(tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + double(f>spikeThreshold);
  end
  
  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    if ~useRealTrajectory
      figure(h);
      subplot(121);
      plot(fhist(1:tind));
      title('Activity');
      xlabel('Time (s)')
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f s',t))
      subplot(122);
      plot(x(1:tind),y(1:tind))
      hold on;
      if ~isempty(spikeCoords)
        cmap = jet;
        cmap = [cmap((end/2+1):end,:); cmap(1:end/2,:)];
        phaseInds = mod(spikePhases,2*pi)*(length(cmap)-1)/2/pi;
        pointColors = cmap(ceil(phaseInds)+1,:);
  
        scatter3(spikeCoords(:,1), ...
                 spikeCoords(:,2), ...
                 zeros(size(spikeCoords(:,1))), ...
                 30*ones(size(spikeCoords(:,1))), ...
                 pointColors, ...
                 'o','filled');
      end
      axis square
      title({'Trajectory (blue) and',...
             'spikes (colored by theta phase',...
             'blues before baseline peak, reds after)'})
      drawnow
    else
      figure(h);
      subplot(131);
      plot((0:tind-1)*dt,fhist(1:tind));
      hold on;
      plot([0 tind-1]*dt,[spikeThreshold spikeThreshold],'r')
      title('Activity (blue) and threshold (red)');
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
        cmap = jet;
        cmap = [cmap((end/2+1):end,:); cmap(1:end/2,:)];
        phaseInds = mod(spikePhases,2*pi)*(length(cmap)-1)/2/pi;
        pointColors = cmap(ceil(phaseInds)+1,:);
  
        scatter3(spikeCoords(:,1), ...
                 spikeCoords(:,2), ...
                 zeros(size(spikeCoords(:,1))), ...
                 30*ones(size(spikeCoords(:,1))), ...
                 pointColors, ...
                 'o','filled');
      end
      axis square
      title({'Trajectory (blue) and',...
             'spikes (colored by theta phase',...
             'blues before baseline peak, reds after)'})
      drawnow
    end
  end  
end
