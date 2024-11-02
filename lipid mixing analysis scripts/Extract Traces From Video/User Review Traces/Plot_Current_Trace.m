function [FigureHandles] = Plot_Current_Trace(FigureHandles,CurrentVirusData,UniversalData,...
    CurrentTraceBackSub,PlotCounter,CurrentTraceNumber)

FusionData = CurrentVirusData.FusionData;
DockingData = CurrentVirusData.DockingData;
TraceGradData = CurrentVirusData.TraceGradData;

% Plot the current trace to the current subplot axis
set(FigureHandles.MasterWindow,'CurrentAxes',FigureHandles.SubHandles(PlotCounter));
    cla
    
    hold on

if strcmp(FusionData.Designation,'2 Fuse')
    
    plot(CurrentTraceBackSub,'r-')
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
    
    XToPlot = [FusionData.FuseFrameNumbers(1), FusionData.FuseFrameNumbers(1)];
    plot(XToPlot,LineToPlot,'g--')
    XToPlot = [FusionData.FuseFrameNumbers(2), FusionData.FuseFrameNumbers(2)];
    plot(XToPlot,LineToPlot,'m--')

elseif strcmp(FusionData.Designation,'1 Fuse')
    
    plot(CurrentTraceBackSub,'b-')
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
    
    XToPlot = [FusionData.FuseFrameNumbers(1), FusionData.FuseFrameNumbers(1)];
    plot(XToPlot,LineToPlot,'g--')
    
elseif strcmp(FusionData.Designation,'No Fusion')
    
    plot(CurrentTraceBackSub,'k-')
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
    
elseif strcmp(FusionData.Designation,'Slow')
    
    plot(CurrentTraceBackSub,'r-')
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
    
    TraceRunMedian = TraceGradData.TraceRunMedian;
    DiffPosClusterData = TraceGradData.DiffPosClusterData;
    SlowFusePosFrameNumbers = TraceRunMedian.FrameNumbers(DiffPosClusterData.ClusterStartIndices);
    for d = 1:length(SlowFusePosFrameNumbers)
        XToPlot = [SlowFusePosFrameNumbers(d), SlowFusePosFrameNumbers(d)];
        plot(XToPlot,LineToPlot,'c--')
    end
    
else
    
    plot(CurrentTraceBackSub,'k.')
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
end

    Title = strcat(num2str(PlotCounter),'-',FusionData.Designation,'-',num2str(CurrentTraceNumber));
    title(Title);
         
    XToPlot = [UniversalData.pHDropFrameNumber, UniversalData.pHDropFrameNumber];
    plot(XToPlot,LineToPlot,'k--')
    

if strcmp(DockingData.IsMobile ,'y')
    XToPlot = [DockingData.StopFrameNum, DockingData.StopFrameNum];
    plot(XToPlot,LineToPlot,'b--')
end
     
hold off
drawnow
end