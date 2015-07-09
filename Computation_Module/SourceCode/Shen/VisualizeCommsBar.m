function VisualizeCommsBar(NodeLabels,NumComms,MakeNewFigure)
% Ben Fulcher, 2014-04-17

BarHeight = 1;
BarWidth = 1;
NumNodes = length(NodeLabels);

if MakeNewFigure
    figure('color','w'); box('on'); hold on
end


NumOverlapping = sum(cellfun(@(x)length(x)>1,NodeLabels));

% Set colors
if NumComms+1 <= 9
    Colors = BF_getcmap('set1',NumComms+1,0);
elseif NumComms+1 <= 21
    Colors = [BF_getcmap('set1',NumComms+1,0);BF_getcmap('set3',NumComms+1,0)];
elseif NumComms+1 <= 64
    Colors = jet(NumComms+1);
else
    error('%u is too many communities for me to color :(',NumComms);
end

% Plot community memberships
for j = 1:NumComms
    for k = 1:NumNodes
        if any(NodeLabels{k}==j);
            if length(NodeLabels{k})>1
                rectangle('Position',[k-BarWidth/2,j-BarHeight/2,BarWidth,BarHeight],'FaceColor', ...
                                Colors(j,:)*0.8,'EdgeColor',Colors(j,:)*0.8,'LineWidth',0.02)
            else
                rectangle('Position',[k-BarWidth/2,j-BarHeight/2,BarWidth,BarHeight],'FaceColor', ...
                                Colors(j,:),'EdgeColor',Colors(j,:),'LineWidth',0.01)
            end
        end
    end
end

% Add summary barcode
for k = 1:NumNodes
    if length(NodeLabels{k})==1;
        TheComm = NodeLabels{k};
        rectangle('Position',[k-BarWidth/2,NumComms+2-BarHeight/2,BarWidth,BarHeight],'FaceColor', ...
                        Colors(TheComm,:),'EdgeColor',Colors(TheComm,:),'LineWidth',0.01)
    else
        rectangle('Position',[k-BarWidth/2,NumComms+2-BarHeight/2,BarWidth,BarHeight],'FaceColor', ...
                        Colors(NumComms+1,:),'EdgeColor',Colors(NumComms+1,:),'LineWidth',0.01)
    end
    
end

plot([1,NumNodes],(NumComms+1)*ones(2,1),'k','LineWidth',2)

xlim([0.5,NumNodes+0.5])
if NumComms>0
    ylim([0.5,NumComms+2+0.5])
end

end