% Burgess 2008's oscillatory interference models
% eric zilli - 20110909 - v1.0
% 
% Burgess 2008 described a wide chunk of theory relating to oscillatory
% interference grid cell models, both reviewing the previous work on
% these models and pointing out many simple variations on the basic model.
% 
% The variations include:
%   * Additive vs. multiplicative interactions of the active
%       oscillators/VCOs
%   * Sinusoidal or "punctate" (spike-like) oscillator output
%   * Combining oscillators with or without leaky integration
%   * Using baseline modulation or not (showing baseline modulation is
%       necessary, at least when the oscillations are summed rather than
%       multiplied, cf. Hasselmo2008.m)
%   * Directional or nondirectional VCOs (whether VCOs output when
%       preferred direction is more than 90 degrees from current heading)
% and we provide options to simulate each of these variations.
%
% Importantly: the grid cell threshold depends on the exact combination of
% options used, and while I started tryin to pre-program them all, there
% are just too many! It is very easy to find the threshold needed for
% any particular case though:
%   * First, set the simulation duration to 10 seconds or so and run it.
%   * Then do figure; plot(fhist) if it isn't already plotted and select
%       an appropriate threshold value a bit below the highest peaks sin
%       that plot but well above the lowest peaks.
%
% The manuscript suggests one other variation where the frequency
% of the baseline oscillation is the average of the frequencies
% of the active VCOs over all preferred directions. It is not clear though
% how one oscillator could take the average of the frequencies of a number
% of other oscillators and use that as its own frequency.
%
% The manuscript includes a precession model using 6 directional VCOs
% (see others in BurgessEtAl2007_precession.m and Hasselmo2008.m). The fix
% isn't so plausible, though. The directional fix used here blocks the
% output from a VCO onto a grid cell if the animal's heading is not within
% 90 degrees of the VCO's preferred direction, but the VCO still integrates
% velocity even when its output is blocked. In theory this could be
% implemented, e.g.
%   A. as a cell that has nonlinear limit cycle subthreshold oscillations
%       that transition into spiking at a frequency equal to the baseline
%       frequency, but it seems unlikely that such an oscillator could
%       have its frequency correctly controlled in the post-spike period,
%       so it may not work.
%   B. the output from this cell could pass through an intermediate cell
%       that is inhibited by the opposite direction head direction input,
%       but adding in extra ad hoc circuiry to produce phase precession
%       is not particularly compelling
%
% This manuscript also fixes a problem with the Burgess et al. 2007 model
% where the product over the pairs of oscillators was being
% thresholded as a whole, rather than its individual components
% being thresholded as is done here. Though really it is only a problem in
% as much as it lets two negative oscillations multiply together into a
% positive value, which does not seem biologically plausible.
%
% Finally, this manuscript presents the model's equations in a form
% suitable for arbitrary input trajectories (though Giocomo et al. 2007
% was probably historically the earliest paper to do so).
%
% Note: This code is not necessarily optimized for speed, but is meant as
% a transparent implementation of the model as described in the manuscript.
%
% This code is released into the public domain. Not for use in skynet.

% if >0, plots the sheet of activity during the simulation on every livePlot'th step
livePlot = 200;

% are the oscillators combined with a sum or product?
oscillatorInteraction = 'sum'; % 'sum' or 'product'
% are the active VCOs or baseline sinusoidal or spike/delta-shaped?
VCOOutput = 'sine'; % 'sine' or 'spike'
baselineOutput = 'sine'; % 'sine' or 'spike'

% does the grid cell reflect only its immediate inputs or does it integrate
% them over time? integration makes more sense with
% oscillatorInteraction = 'sum' but if it is 'prod', we'll multiply the
% oscillations first and let the result be leakily-integrated. NB that
% though Burgess states that leaky integration is a reasonable assumption
% for entorhinal neurons, it is specifically untrue for stellate cells,
% which are famously resonant and so have a more complex response than
% simple integration.
integrateGridInputs = 0; % if =1, grid cell is a leaky integrator; if =0, grid cell only responds to immediate inputs
baselineModulation = 1; % if =0, no baseline oscillation is used; if =1, a baseline is used
directionalVCOs = 0; % if =0, VCOs output to grid cell for any heading; if =1, VCOs only output if heading is within 90 degrees of their preferred direction

nVCOs = 3;

% Threshold to determine when the grid cell has spiked
spikeThreshold = 3.5;

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
% Baseline maintains a fixed frequency
baseFreq = 6; % Hz
% Directional preference of each VCO (this also sets the number of VCOs)
dirPreferences = (0:nVCOs-1)*pi/3;
% Scaling factor relating speed to oscillator frequencies
% NB paper uses 0.05*2pi rad/cm. But we do the conversion to rad later,
% leaving 0.05 Hz/cm = 5 Hz/m which produces very tight field spacing. For
% cosmetic purposes for the trajectory we use here, we'll use beta = 2.
beta = 2; % Hz/(m/s) 
T = 0.010; % s, time constant for integration of oscillations if used
C = 1/0.5; % normalization coefficient. the integral works out to be pi*(2n-3)!!/(2^(n-1)*(n-1)!) by my reckoning
if integrateGridInputs && dt>(0.05*T)
  warning('Your time step dt should be much smaller than the time constant T!')
  dt = 0.05*T;
  fprintf('Using dt = %g, but frankly it should be smaller!\n',dt)
end

%% History variables
speed = zeros(1,ceil(simdur/dt));
curDir = zeros(1,ceil(simdur/dt));
vhist = zeros(1,ceil(simdur/dt));
fhist = zeros(1,ceil(simdur/dt));
f = 0;

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
VCOPhases = zeros(1,length(dirPreferences)); % rad
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
  speed(tind) = sqrt(v(1)^2+v(2)^2); % m/s
  
  x(tind) = x(tind-1)+v(1)*dt; % m
  y(tind) = y(tind-1)+v(2)*dt; % m
    
  % VCO frequencies are pushed up or down from the baseline frequency
  % depending on the speed and head direction, with a scaling factor beta
  % that sets the spacing between the spatial grid fields.
  VCOFreqs = baseFreq + beta*speed(tind)*cos(curDir(tind)-dirPreferences); % Hz
  
  % Advance oscillator phases
  % Radial frequency (2pi times frequency in Hz) is the time derivative of phase.
  VCOPhases = VCOPhases + dt*2*pi*VCOFreqs; % rad
  basePhase = basePhase + dt*2*pi*baseFreq; % rad
  
  % Compute active oscillation values from current phase
  if strcmp(VCOOutput,'sine')
    VCOs = cos(VCOPhases);
  elseif strcmp(VCOOutput,'spike')
    VCOs = (1/2+cos(VCOPhases)/2).^50;
  end
  
  if directionalVCOs
    % "Directional" VCOs do not output unless they
    % are within 90 degress of the VCO's preferred direction
    VCOs = VCOs.*(cos(curDir(tind)-dirPreferences)>0);
  end
  
  % Compute baseline oscillation from current phase
  if strcmp(baselineOutput,'sine')
    baseline = baselineModulation*cos(basePhase);
  elseif strcmp(baselineOutput,'spike')
    baseline = baselineModulation*(1/2+cos(VCOPhases)/2).^50;
  end
  
  % Sum each dendritic oscillation separately with the baseline oscillation
  VCOsPlusBaseline = VCOs + baseline;
  
  % Threshold individual oscillator pairs
  VCOsPlusBaseline = VCOsPlusBaseline.*(VCOsPlusBaseline>0);
  
  % Final activity is the sum or product of the thresholded oscillations.
  if strcmp(oscillatorInteraction,'sum')
    act = sum(VCOsPlusBaseline);
  elseif strcmp(oscillatorInteraction,'product')
    act = prod(VCOsPlusBaseline);
  end
  
  if integrateGridInputs
    % final output f is a leaky integration of the combined oscillations
    f = f + dt*(-f/T + act*C);
  else
    % final output f is simply the current value of the combined oscillations
    f = act;
  end
  
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
