function [FigureHandles] = Draw_Lines_on_Plot(FigureHandles,DockingData,FusionData,UniversalData)

set(0,'CurrentFigure',FigureHandles.TraceWindow)
hold on
LineToPlot = ylim;

if strcmp(FusionData.Designation,'2 Fuse')
    XToPlot = [FusionData.FuseFrameNumbers(1), FusionData.FuseFrameNumbers(1)];
    plot(XToPlot,LineToPlot,'g--')
    XToPlot = [FusionData.FuseFrameNumbers(2), FusionData.FuseFrameNumbers(2)];
    plot(XToPlot,LineToPlot,'k--')

    Title = strcat('pH = ',num2str(UniversalData.pHDropFrameNumber),...
        '; Dock = ',num2str(DockingData.StopFrameNum),...
        '; 1fuse = ', num2str(FusionData.FuseFrameNumbers(1)),...
        '; 2fuse = ', num2str(FusionData.FuseFrameNumbers(2)),...
        '; pHtoF = ', num2str(FusionData.pHtoFusionTime(1)));
    title(Title);
    
elseif strcmp(FusionData.Designation,'1 Fuse')
    XToPlot = [FusionData.FuseFrameNumbers(1), FusionData.FuseFrameNumbers(1)];
    plot(XToPlot,LineToPlot,'g--')
    Title = strcat('pH = ',num2str(UniversalData.pHDropFrameNumber),...
        '; Dock = ',num2str(DockingData.StopFrameNum),...
        '; 1fuse = ', num2str(FusionData.FuseFrameNumbers(1)),...
        '; pHtoF = ', num2str(FusionData.pHtoFusionTime(1)));
    title(Title);
    
elseif strcmp(FusionData.Designation,'No Fusion')
    Title = strcat('pH = ',num2str(UniversalData.pHDropFrameNumber),...
        '; Dock = ',num2str(DockingData.StopFrameNum));
    title(Title);
end

if strcmp(DockingData.IsMobile ,'y')
    XToPlot = [DockingData.StopFrameNum, DockingData.StopFrameNum];
    plot(XToPlot,LineToPlot,'b--')
end
     
hold off
drawnow
end