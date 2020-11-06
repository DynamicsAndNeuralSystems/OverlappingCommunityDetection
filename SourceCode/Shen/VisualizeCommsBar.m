function VisualizeCommsBar(nodeLabels,numComms,makeNewFigure)
% Visualize communities as a bar chart
% Ben Fulcher, 2014-04-17
%-------------------------------------------------------------------------------

barHeight = 1;
barWidth = 1;
numNodes = length(nodeLabels);

if makeNewFigure
    figure('color','w');
    box('on'); hold('on')
end

numOverlapping = sum(cellfun(@(x)length(x)>1,nodeLabels));

% Set colors
if numComms+1 <= 9
    Colors = BF_GetColorMap('set1',numComms+1,0);
elseif numComms+1 <= 21
    Colors = [BF_GetColorMap('set1',numComms+1,0);BF_GetColorMap('set3',numComms+1,0)];
elseif numComms+1 <= 64
    Colors = jet(numComms+1);
else
    error('%u is too many communities for me to color :(',numComms);
end

% Plot community memberships
for j = 1:numComms
    for k = 1:numNodes
        if any(nodeLabels{k}==j);
            if length(nodeLabels{k})>1
                rectangle('Position',[k-barWidth/2,j-barHeight/2,barWidth,barHeight],'FaceColor', ...
                                Colors(j,:)*0.8,'EdgeColor',Colors(j,:)*0.8,'LineWidth',0.02)
            else
                rectangle('Position',[k-barWidth/2,j-barHeight/2,barWidth,barHeight],'FaceColor', ...
                                Colors(j,:),'EdgeColor',Colors(j,:),'LineWidth',0.01)
            end
        end
    end
end

% Add summary barcode
for k = 1:numNodes
    if length(nodeLabels{k})==1;
        theComm = nodeLabels{k};
        rectangle('Position',[k-barWidth/2,numComms+2-barHeight/2,barWidth,barHeight],'FaceColor', ...
                        Colors(theComm,:),'EdgeColor',Colors(theComm,:),'LineWidth',0.01)
    else
        rectangle('Position',[k-barWidth/2,numComms+2-barHeight/2,barWidth,barHeight],'FaceColor', ...
                        Colors(numComms+1,:),'EdgeColor',Colors(numComms+1,:),'LineWidth',0.01)
    end
    
end

plot([1,numNodes],(numComms+1)*ones(2,1),'k','LineWidth',2)

xlim([0.5,numNodes+0.5])
if numComms > 0
    ylim([0.5,numComms+2+0.5])
end

end
