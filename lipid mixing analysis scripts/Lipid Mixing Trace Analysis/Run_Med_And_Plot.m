function [TraceRunMedian,FigureHandles] = ...
    Run_Med_And_Plot(CurrTraceCropped,FigureHandles,UniversalData,Options)
    
%Calc the running median, plot running median and where the pH drop
%occurred
    if strcmp(Options.TypeofFusionData, 'TetheredVesicle') 
        RunMedHalfLength = 1; 
            %Num of points on either side, prev value = 5
    elseif strcmp(Options.TypeofFusionData, 'SLBSelfQuench') 
        RunMedHalfLength = 0; 
            %Num of points on either side, prev value = 5
    else
        disp(' Type of fusion data not specified correctly');
    end        

        StartIdx = RunMedHalfLength + 1;
        EndIdx = length(CurrTraceCropped.Trace)-RunMedHalfLength;
    TraceRunMedian.Trace = zeros(length(StartIdx:EndIdx),1);
        %We also set up a vector with the actual frame numbers
        %corresponding to each index position.
        TraceRunMedian.FrameNumbers = CurrTraceCropped.FrameNumbers(StartIdx:EndIdx);

    for n = StartIdx:EndIdx
        TraceRunMedian.Trace(n-RunMedHalfLength) = median(CurrTraceCropped.Trace(n-RunMedHalfLength:n+RunMedHalfLength));
    end


%Plot running median
    set(0,'CurrentFigure',FigureHandles.TraceWindow)
    cla
    plot(CurrTraceCropped.FrameNumbers,CurrTraceCropped.Trace,'b-')
    hold on
    plot(TraceRunMedian.FrameNumbers,TraceRunMedian.Trace,'r-');
    %Plot a line where the pH drop event occurred
        set(0,'CurrentFigure',FigureHandles.TraceWindow)
        hold on
            LineToPlot = ylim;
            XToPlot = [UniversalData.pHDropFrameNumber, UniversalData.pHDropFrameNumber];
        plot(XToPlot,LineToPlot,'m--')
    hold off
end