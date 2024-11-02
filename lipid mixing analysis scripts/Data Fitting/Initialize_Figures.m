function [FigureHandles] = Initialize_Figures(Options)

    FigureHandles.ResidualsWindow=figure(3);
    FigureHandles.FitWindow=figure(2);
    FigureHandles.Fuse1Wind=figure(1);
    FigureHandles.EfficiencyWindow = figure(4);
    
    if strcmp(Options.ShowBeginningIntensity,'y')
        FigureHandles.IntensityWindow = figure(5);
        set(FigureHandles.IntensityWindow,'Position',[500 50 400 400])
    end
    
    if strcmp(Options.RunBootstrap,'y')
        FigureHandles.BootstrapWindow = figure(6);
        set(FigureHandles.BootstrapWindow,'Position',[-5 50 560 420])
    end

%     Fuse2Wind=figure(4);
%     FigureHandles.HistogramWindow = figure(4);

%     set(FigureHandles.FitWindow, 'Position', [562   366   560   420]);
%     set(FigureHandles.ResidualsWindow, 'Position', [1124  368  560  420]);
%     set(FigureHandles.Fuse1Wind, 'Position', [-5   364   560   420]);
%     
    set(FigureHandles.FitWindow, 'Position', [-5   364   560   420]);
    set(FigureHandles.ResidualsWindow, 'Position', [562   366   560   420]);
    set(FigureHandles.EfficiencyWindow, 'Position', [530    57   548   341]);
    set(FigureHandles.Fuse1Wind, 'Position', [1124  368  560  420]);
    

end