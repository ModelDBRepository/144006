% figure demonstrating read-out mechanisms in interference models
% eric zilli - 20111105 - v1.0
%
% Read-out rules that have been used:
% * thresholdlinear[(base+a1)*(base+a2)*(base+a3)] Burgess et al. 2007
% * thresholdlinear[(base+a1)+(base+a2)+(base+a3)] Burgess et al. 2007
% * (s1+s2) Blair et al. 2007
% * s1*s2 Gaussier et al. 2007
% * threshold(a1>0.5)*threshold(a2>0.5)*threshold(a3>0.5) Hasselmo 2008
% * (1D) (base+a1) Burgess 2008 
% * (1D) base*a1 Blair et al. 2008, Burgess 2008
% * thresholdlinear(base+a1)*thresholdlinear(base+a2) Burgess 2008 
% * base*(a1+a2+a3) Burgess 2008 
% * (a1+a2+a3) Burgess 2008 
% * threshold((base+a1+a2)>2) Zilli and Hasselmo 2010
% * threshold(base)*(thresholdlinear(a1)+thresholdlinear(a2)) Zilli and Hasselmo 2010
% * -base - a1 - a2 Zilli and Hasselmo 2010, Welday et al 2011
% * threshold(s1+s2+s3+...) Mhatre et al. 2010
%
% base indicates the baseline oscillation, a1-a3 active oscillations,
% s1-s3 stripe cell activities, thresholdlinear() is the function
% thresholdlinear(x) = 0 if x<=0 and thresholdlinear(x) = x if x>0
% threshold() is the Heaviside function:
% threshold(x)=0 if x<=0 and threshold(x) = 1 if x>0, except for
% Zilli and Hasselmo (2010) where threshold indicates the thresholded
% spiking activity of a single neuron.

%% Variables for plotting oscillations
plotDur = 2*1.5; % s
dt = 0.001;
tvect = dt:dt:plotDur;

baseFreq = 5; % Hz
activeFreqs = 5 + [0.5 -0.5 1];
startPhase = [0 pi pi 0];

% Oscillations:
base = cos(2*pi*baseFreq*tvect + startPhase(1));
a1 = cos(2*pi*activeFreqs(1)*tvect + startPhase(2));
a2 = cos(2*pi*activeFreqs(2)*tvect + startPhase(3));
a3 = cos(2*pi*activeFreqs(3)*tvect + startPhase(4));

%% Evaluate readout rules for sinusoids
Bu07a = heaviside((base+a1).*(base+a2).*(base+a3)).*(base+a1).*(base+a2).*(base+a3); % his eq (5)
Bu07b = heaviside((base+a1)+(base+a2)+(base+a3))+(base+a1)+(base+a2)+(base+a3); % his eq (5)
Bu08e =  (a1+a2+a3); % implied in Fig 7 left
Ha08 = heaviside((a1-0.5)).*heaviside((a2-0.5)).*heaviside((a3-0.5));
Bu08a = (base+a1); % his Fig 1 a (and b)
Bu08b = (base.*a1); % his Fig 1 c (and d); he seemed to threshold it
% Bl08 = base.*a1; % blair et al. 2008 is same rule as burgess 2008b
Bu08c = heaviside((base+a1)).*(base+a1).*heaviside((base+a2)).*(base+a2); % his eq (6) r(t) with n=2
Bu08d = (0.5+0.5*base).*(a1+a2+a3); % his equation 9
Zi10a = heaviside(((base+a1+a2)-2));
Zi10b = heaviside(base-0.5).*(a1+a2);
Zi10c = - base - a1 - a2;
% Zi10c = - base - a1 - a2 + 1.5;
% Zi10b = heaviside(base).*(heaviside(a1).*a1+heaviside(a2).*a2);
% Ga07 = s1.*s2;
% Mh10 = heaviside(s1+s2+s3);

%% Pull variables and labels together to plot in a loop
signals = [base; a1; a2; a3; Bu07a; Bu07b; Bu08e; Ha08; Bu08a; Bu08b; Bu08c; Bu08d; Zi10a; Zi10c];
% signalNames = {'base', 'a1', 'a2', 'a3', 'Bu07', 'Bu08e', 'Ha08', 'Bu08a', 'Bu08b, Bl08', 'Bu08c', 'Bu08d', 'Zi10a', 'Zi10b'};
% signalNames = {'base', 'a1', 'a2', 'a3', 'Burgess et al. 2007b', 'Burgess et al. 2007a', 'Blair et al. 2007', 'Hasselmo 2008', 'Burgess 2008a', {'Blair et al. 2008, Burgess 2008b'}, 'Burgess 2008c', 'Burgess 2008d', 'Zilli et al. 2010a', 'Zilli et al. 2010b','Welday et al. 2011, Zilli et al. 2010c'};
signalNames = {'base', 'a1', 'a2', 'a3', 'Burgess et al. 2007', 'Burgess et al. 2007', 'Burgess et al. 2008', 'Hasselmo 2008', 'Burgess 2008', {'Blair et al. 2008, Burgess 2008'}, 'Burgess 2008', 'Burgess 2008', 'Zilli et al. 2010', 'Welday et al. 2011, Zilli et al. 2010'};
shortEqs = {'R[(b + a_1)(b + a_2)(b + a_3)]', ...
            'R[(b + a_1)+(b + a_2)+(b + a_3)]', ...
            'a_1 + a_2 + a_3', ...
            'H(a_1 - 0.5)H(a_2 - 0.5)H(a_3 - 0.5)', ...
            'b + a_1', ...
            'ba_1', ...
            'R(b + a_1)R(b + a_2)', ...
            '(b+1)(a_1 + a_2 + a_3)/2', ...
            'H(b + a_1 + a_2 - 2)', ...
            '-b - a_1 - a_2'};

permutation = [1 2 3 4 9 10 11 7 12 6 5 8 13];

%% Optionally re-order signals before plotting
tsignals = signals;
for i=1:length(permutation)
  signals(i,:) = tsignals(permutation(i),:);
end
tsignalNames = signalNames;
for i=1:length(permutation)
  signalNames{i} = tsignalNames{permutation(i)};
end
tshortEqs = shortEqs;
for i=1:(length(permutation)-4)
  shortEqs{i} = tshortEqs{permutation(i+4)-4};
end


%% Figure options
nrows = size(signals,1)-3;
leftMargin = 0.05;
bottomMargin = 0.06;
width = 0.89/2;
height = .75/nrows;
lefts = leftMargin + [0 0.5];
bottoms = bottomMargin + 1.04*linspace(0,1,nrows+2);

set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontSize', 8)
set(0,'defaultTextFontName', 'Arial')

% size on paper:
widthOnPaper = 18.0; % cm
heightOnPaper = 13.5; % cm

figure('units','centimeters','position',[1 1 widthOnPaper heightOnPaper],'color','w');
set(gcf, 'renderer', 'painter')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [widthOnPaper heightOnPaper]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 widthOnPaper heightOnPaper]);

%% Colors of oscillations
% % cmykish
% colors = [0.8 0.75 0.3; 0.7 0.45 0.98; 0.1 0.7 0.6];
% % blue/orange/yellow
% colors = [0.3 0.2 0.8; 0.8 0.3 0.15; 0.85 0.9 0.3];
% standard rgb
colors = [1 0 0; 0 0.5 0; 0 0 1];
% muted rgb
colors = [.9 0 0; 0 0.4 0; 0 0 0.9];

%% Optional line to separate oscillations from read-outs
% lineY = 0.925;
% annotation('line',[.05 .45],[lineY lineY],'color',[0.6 0.7 0.6])
% annotation('line',[.58 .97],[lineY lineY],'color',[0.6 0.7 0.6])

%% Plot sinusoid column
% figure('name','Read-out','color','w','position',[181 146 451 649])
for ind=4:size(signals,1)
  axes('position',[lefts(1) bottoms(size(signals,1)-ind+1) width height]);
  if ind==4
    set(gca,'colororder',colors,'nextplot','replacechildren')
    plot(tvect,signals(2:4,:)')
    hold on;
    plot(tvect,signals(1,:),'k');
    text(-0.075*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
         0.96*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
         '(A)',...
         'FontSize',10,...
         'FontWeight','bold',...
         'HorizontalAlignment','center')
  else
    plot(tvect,signals(ind,:),'k');
  end
  set(gca,'box','off')
  
  % Clean up yticks and y limits of plots
  if min(signals(ind,:))<0
    set(gca,'ylim',1.1*[min(signals(ind,:)) max(signals(ind,:))])
  else
    set(gca,'ylim',[.91*min(signals(ind,:)) 1.1*max(signals(ind,:))])
  end
  if ind==7 % Bu08c
    set(gca,'ytick',[0 4])
  elseif ind==8 % Bu08e
    set(gca,'ytick',[-3 0 3])
    ylim([-3.2 3])
  elseif ind==9 % Bu08d
    ylim([-2.1 4.1])
    set(gca,'ytick',[-2 0 2 4])
%     set(gca,'ytick',[-1 1 3])
  elseif ind==11 % Bu07b
    set(gca,'ytick',[0 8])
  elseif ind==12 % Ha08
    set(gca,'ytick',[0 1])
  elseif ind==13 % Zi10a
    set(gca,'ytick',[0 1])
  elseif ind==14 % Zi10c
    set(gca,'ytick',[-3 0 3])
  end
%   if min(get(gca,'ylim'))==0
%     ylim([-0.05 max(get(gca,'ylim'))])
%   end
  
  % Turn off most tick labels and xlabel the bottom plot "Time (s)"
  if ind~=size(signals,1)
    set(gca,'xticklabel',[])
  else
    text(0.5*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
         -.63*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
         'Time (s)',...
         'FontSize',9,...
         'HorizontalAlignment','center')
  end
  if ind~=4
    text(2,0.9*diff(get(gca,'ylim'))+min(get(gca,'ylim')),shortEqs{ind-4},'fontsize',8,'horizontalalign','center')
  end
end


%% Define EPSP oscillations
synKernel = exp(-tvect/0.060); % 10 ms decay
baseTimes = 0:1/baseFreq:max(tvect);
a1Times = 0:1/activeFreqs(1):max(tvect);
a2Times = 0:1/activeFreqs(2):max(tvect);
a3Times = 0:1/activeFreqs(3):max(tvect);

tvect = tvect+max(tvect);

base = zeros(1,length(tvect));
base(ceil(baseTimes/dt)+1) = 1;
base = conv(base,synKernel);
base = base(1:length(tvect));

a1 = zeros(1,length(tvect));
a1(ceil(a1Times/dt)+1) = 1;
a1 = conv(a1,synKernel);
a1 = a1(1:length(tvect));

a2 = zeros(1,length(tvect));
a2(ceil(a2Times/dt)+1) = 1;
a2 = conv(a2,synKernel);
a2 = a2(1:length(tvect));

a3 = zeros(1,length(tvect));
a3(ceil(a3Times/dt)+1) = 1;
a3 = conv(a3,synKernel);
a3 = a3(1:length(tvect));

%% Evaluate readout rules for EPSPs
Bu07a = heaviside((base+a1).*(base+a2).*(base+a3)).*(base+a1).*(base+a2).*(base+a3); % his eq (5)
Bu07b = heaviside((base+a1)+(base+a2)+(base+a3))+(base+a1)+(base+a2)+(base+a3); % his eq (5)
Bu08e =  (a1+a2+a3); % their eq (1), though it wasn't really a model per se
Ha08 = heaviside((a1-0.5)).*heaviside((a2-0.5)).*heaviside((a3-0.5));
Bu08a = (base+a1); % his Fig 1 a (and b)
Bu08b = (base.*a1); % his Fig 1 c (and d); he seemed to threshold it
% Bl08 = base.*a1; % blair et al. 2008 is same rule as burgess 2008b
Bu08c = heaviside((base+a1)).*(base+a1).*heaviside((base+a2)).*(base+a2); % his eq (6) r(t) with n=2
Bu08d = (0.5+0.5*base).*(a1+a2+a3); % his equation 9
Zi10a = heaviside(((base+a1+a2)-2));
% Zi10b = heaviside(base-0.5).*(a1+a2);
Zi10c = - base - a1 - a2;
% Zi10c = - base - a1 - a2 + 1.5;

signals = [base; a1; a2; a3; Bu07a; Bu07b; Bu08e; Ha08; Bu08a; Bu08b; Bu08c; Bu08d; Zi10a; Zi10c];
tsignals = signals;
for i=1:length(permutation)
  signals(i,:) = tsignals(permutation(i),:);
end

%% Plot EPSPs and readout
for ind=4:size(signals,1)
  axes('position',[lefts(2) bottoms(size(signals,1)-ind+1) width height]);
  if ind==4
    set(gca,'colororder',colors,'nextplot','replacechildren')
    plot(tvect,signals(2:4,:)')
    hold on;
    plot(tvect,signals(1,:),'k');
    set(gca,'ytick',[0 0.5 1])
  else
    plot(tvect,signals(ind,:),'k');
  end
  set(gca,'box','off')
  
  % Clean up yticks and y limits of plots
  if ind>4
    if min(signals(ind,:))<0
      set(gca,'ylim',1.1*[min(signals(ind,:)) max(signals(ind,:))])
    elseif min(signals(ind,:))<.2
      set(gca,'ylim',[0 1.1*max(signals(ind,:))])
    else
      set(gca,'ylim',[.91*min(signals(ind,:)) 1.1*max(signals(ind,:))])
    end
  else
    ylim([-0.01 1.05])
  end
  if ind==14
    ylim([-3 0.75]);
  end
  if ind==4
    set(gca,'ytick',[0 1]);
  elseif ind==5
    set(gca,'ytick',[0 2]);
  elseif ind==6
    set(gca,'ytick',[0 1]);
  elseif ind==7
    set(gca,'ytick',[0 4]);
  elseif ind==8
    set(gca,'ytick',[0 3]);
  elseif ind==9
    set(gca,'ytick',[0 3]);
  elseif ind==10
    set(gca,'ytick',[0 6]);
    ylim([0 7.5])
  elseif ind==11
    set(gca,'ytick',[0 8]);
  elseif ind==12
    set(gca,'ytick',[0 1]);
  elseif ind==13
    set(gca,'ytick',[0 1]);
  elseif ind==14
    set(gca,'ytick',[-3 0]);
  end
  
  % Panel label
  if ind==4
    text(-0.085*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
      .92*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
      '(B)',...
      'FontSize',10,...
      'FontWeight','bold',...
      'HorizontalAlignment','center');
  end
  
  % Turn off most tick labels and xlabel the bottom plot "Time (s)"
  if ind~=size(signals,1)
    set(gca,'xticklabel',[])
  else
    text(0.5*diff(get(gca,'xlim'))+min(get(gca,'xlim')),...
         -0.63*diff(get(gca,'ylim'))+min(get(gca,'ylim')),...
         'Time (s)',...
         'FontSize',9,...
         'HorizontalAlignment','center')
  end
  if ind>=5
    text(4,0.9*diff(get(gca,'ylim'))+min(get(gca,'ylim')),signalNames{ind},'fontsize',8,'horizontalalign','center')
  end
end
