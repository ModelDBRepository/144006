% draw a ring attractor into given axes
% eric zilli - 20111114 - v1.0
% drawSimpleRing_0a(10,1,0.1,0.01,1,repmat([0.5 0.5 0.5],10,1))
function drawSimpleRing_0a(nDirs,nRings,ringRadii,cellWidth,highlightDir,ringFillColors,drawArrow)

% % Which position will be drawn as most active or will have its synapses
% % illustrated?
% highlightDir = ceil(nDirs/4)+1;

dirs = (0:nDirs-1)/nDirs*2*pi;

ringEdgeColors = [0 0 0]; % [0.7 0.7 0.7; 0 0 0; 0.7 0.7 0.7];
% ringFillColors = [1 1 1; 1 1 1; 1 1 1];

nSymOffsets = 0; 4;
symLineColor = repmat([0 0 0],nSymOffsets,1);
symLineWidth = linspace(1.0,0.5,nSymOffsets);

nAsymOffsets = 0; 4;
asymLineColor = repmat([0 0 0],nAsymOffsets,1);
asymLineWidth = linspace(1,0.25,nAsymOffsets);


for ring=1:nRings
  for d=1:nDirs
    lowerLeftX = ringRadii(ring)*cos(dirs(d))-cellWidth/2;
    lowerLeftY = ringRadii(ring)*sin(dirs(d))-cellWidth/2;
    rectangle('Position',[lowerLeftX,lowerLeftY,cellWidth,cellWidth],...
              'Curvature',[1 1],...
              'FaceColor',ringFillColors(1+mod(d-highlightDir,nDirs),:));%,...
  end
  %% Draw symmetric connections within one ring
  if ring==(length(ringRadii)+1)/2
    for offset=1:nSymOffsets
      arc = cellWidth/2+ringRadii(ring)+offset/2*ringRadii(ring)/2*sin(linspace(0,pi,50));
      angles = linspace(dirs(highlightDir),dirs(highlightDir)+2*pi*offset/nDirs,50);
      h = polar(angles,arc);
      set(h,'Color',symLineColor(offset,:),'LineWidth',symLineWidth(offset));
      % find slope of last bit of synaptic line so we can draw two lines at 25 degrees to make an arrow
      arrowLen = 1.25*cellWidth;
      arrowAngle = deg2rad(25);
      endX1 = arc(end-1)*cos(angles(end-1));
      endY1 = arc(end-1)*sin(angles(end-1));
      endX2 = arc(end)*cos(angles(end));
      endY2 = arc(end)*sin(angles(end));
      endAngle = atan2(endY2-endY1,endX2-endX1);
      line([endX2 endX2+arrowLen*cos(endAngle+pi-arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle+pi-arrowAngle)],'LineWidth',symLineWidth(offset),'Color',symLineColor(offset,:))
      line([endX2 endX2+arrowLen*cos(endAngle-pi+arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle-pi+arrowAngle)],'LineWidth',symLineWidth(offset),'Color',symLineColor(offset,:))
      
      angles = linspace(dirs(highlightDir),dirs(highlightDir)-2*pi*offset/nDirs,50);
      h = polar(angles,arc);
      set(h,'Color',symLineColor(offset,:),'LineWidth',symLineWidth(offset));
      % find slope of last bit of synaptic line so we can draw two lines at 25 degrees to make an arrow
      arrowLen = 1.25*cellWidth;
      arrowAngle = deg2rad(25);
      endX1 = arc(end-1)*cos(angles(end-1));
      endY1 = arc(end-1)*sin(angles(end-1));
      endX2 = arc(end)*cos(angles(end));
      endY2 = arc(end)*sin(angles(end));
      endAngle = atan2(endY2-endY1,endX2-endX1);
      line([endX2 endX2+arrowLen*cos(endAngle+pi-arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle+pi-arrowAngle)],'LineWidth',symLineWidth(offset),'Color',symLineColor(offset,:))
      line([endX2 endX2+arrowLen*cos(endAngle-pi+arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle-pi+arrowAngle)],'LineWidth',symLineWidth(offset),'Color',symLineColor(offset,:))
    end
  end
  %% Draw motion arrow outside outer ring
  if ring==nRings && drawArrow
    % Note, drawing these with patch to avoid discontinuities at ends of lines
    hold on;
    arrowColor = [0 0 0];
    arrowWidth = 0.5;
    
    % draw arrow stem
    arc = (ringRadii(ring)+1.5*cellWidth)*ones(1,50);
    angles = linspace(dirs(highlightDir)-2*pi/nDirs/2,dirs(highlightDir)+2*pi/nDirs/2,50);
    h = patch([arc arc(end:-1:1)].*cos([angles angles(end:-1:1)]),[arc arc(end:-1:1)].*sin([angles angles(end:-1:1)]),'k');
    set(h,'faceColor',arrowColor,'LineWidth',arrowWidth);
    
    % draw arrow point outer line
    arc = (ringRadii(ring)+1.5*cellWidth)*ones(1,50)+cellWidth/3*linspace(0,1,50);
    angles = linspace(dirs(highlightDir)+2*pi/nDirs/2,dirs(highlightDir)+1/3*2*pi/nDirs/2,50);
    h = patch([arc arc(end:-1:1)].*cos([angles angles(end:-1:1)]),[arc arc(end:-1:1)].*sin([angles angles(end:-1:1)]),'k');
    set(h,'faceColor',arrowColor,'LineWidth',arrowWidth);

    % draw arrow point inner line
    arc = (ringRadii(ring)+1.5*cellWidth)*ones(1,50)-cellWidth/3*linspace(0,1,50);
    angles = linspace(dirs(highlightDir)+2*pi/nDirs/2,dirs(highlightDir)+1/3*2*pi/nDirs/2,50);
    h = patch([arc arc(end:-1:1)].*cos([angles angles(end:-1:1)]),[arc arc(end:-1:1)].*sin([angles angles(end:-1:1)]),'k');
    set(h,'faceColor',arrowColor,'LineWidth',arrowWidth);
  end

  %% Draw asymmetric connections from one ring to another
  if ring==3
    offset = 0;
    arc = linspace((ringRadii(ring)-ringRadii(ring-1)),0,50) + cellWidth/2 + ringRadii(ring-1) + offset/2*ringRadii(ring)*sin(linspace(0,pi,50));
    angles = linspace(dirs(highlightDir), dirs(highlightDir) + 2*pi*(offset-1)/nDirs, 50);
    h = polar(angles,arc);
    set(h,'Color',asymLineColor(1,:),'LineWidth',1.5);
    % find slope of last bit of synaptic line so we can draw two lines at 25 degrees to make an arrow
    arrowLen = 1.25*cellWidth;
    arrowAngle = deg2rad(25);
    endX1 = arc(end-1)*cos(angles(end-1));
    endY1 = arc(end-1)*sin(angles(end-1));
    endX2 = arc(end)*cos(angles(end));
    endY2 = arc(end)*sin(angles(end));
    endAngle = atan2(endY2-endY1,endX2-endX1);
    line([endX2 endX2+arrowLen*cos(endAngle+pi-arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle+pi-arrowAngle)],'LineWidth',1.0,'Color',asymLineColor(1,:))
    line([endX2 endX2+arrowLen*cos(endAngle-pi+arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle-pi+arrowAngle)],'LineWidth',1.0,'Color',asymLineColor(1,:))

    for offset=1:nAsymOffsets
      %% Counter-clockwise going arrows:
      if offset==1
        arc = linspace((ringRadii(ring)-ringRadii(ring-1)),cellWidth*2.25,50) + cellWidth/2+ringRadii(ring-1)+offset/2*ringRadii(ring)/3*sin(linspace(pi,1.5*pi,50));
        angles = pi/2+.2*linspace(1,0,50).^3.*sin(linspace(0,2*pi,50)); %[linspace(dirs(highlightDir),dirs(highlightDir)+2*pi*(offset)/nDirs,25) linspace(dirs(highlightDir)+2*pi*(offset)/nDirs,dirs(highlightDir),25)];
      else
        arc = linspace((ringRadii(ring)-ringRadii(ring-1)),0,50) + cellWidth/2+ringRadii(ring-1)+offset/2*ringRadii(ring)/3*sin(linspace(0,pi,50));
        angles = linspace(dirs(highlightDir),dirs(highlightDir)+2*pi*(offset-1)/nDirs,50);
      end
      h = polar(angles,arc);
      set(h,'Color',asymLineColor(offset,:),'LineWidth',asymLineWidth(offset));
      % find slope of last bit of synaptic line so we can draw two lines at 25 degrees to make an arrow
      arrowLen = 1.25*cellWidth;
      arrowAngle = deg2rad(25);
      endX1 = arc(end-1)*cos(angles(end-1));
      endY1 = arc(end-1)*sin(angles(end-1));
      endX2 = arc(end)*cos(angles(end));
      endY2 = arc(end)*sin(angles(end));
      endAngle = atan2(endY2-endY1,endX2-endX1);
      line([endX2 endX2+arrowLen*cos(endAngle+pi-arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle+pi-arrowAngle)],'LineWidth',asymLineWidth(offset),'Color',asymLineColor(offset,:))
      line([endX2 endX2+arrowLen*cos(endAngle-pi+arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle-pi+arrowAngle)],'LineWidth',asymLineWidth(offset),'Color',asymLineColor(offset,:))
      
      %% Clockwise going arrows:
      arc = linspace((ringRadii(ring)-ringRadii(ring-1)),0,50) + cellWidth/2+ringRadii(ring-1)+offset/2*ringRadii(ring)/3*sin(linspace(0,pi,50));
      angles = linspace(dirs(highlightDir),dirs(highlightDir)-2*pi*(offset+1)/nDirs,50);
      h = polar(angles,arc);
      set(h,'Color',asymLineColor(offset,:),'LineWidth',asymLineWidth(offset));
      % find slope of last bit of synaptic line so we can draw two lines at 25 degrees to make an arrow
      arrowLen = 1.25*cellWidth;
      arrowAngle = deg2rad(25);
      endX1 = arc(end-1)*cos(angles(end-1));
      endY1 = arc(end-1)*sin(angles(end-1));
      endX2 = arc(end)*cos(angles(end));
      endY2 = arc(end)*sin(angles(end));
      endAngle = atan2(endY2-endY1,endX2-endX1);
      line([endX2 endX2+arrowLen*cos(endAngle+pi-arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle+pi-arrowAngle)],'LineWidth',asymLineWidth(offset),'Color',asymLineColor(offset,:))
      line([endX2 endX2+arrowLen*cos(endAngle-pi+arrowAngle)],[endY2 endY2+arrowLen*sin(endAngle-pi+arrowAngle)],'LineWidth',asymLineWidth(offset),'Color',asymLineColor(offset,:))
    end
  end
end
