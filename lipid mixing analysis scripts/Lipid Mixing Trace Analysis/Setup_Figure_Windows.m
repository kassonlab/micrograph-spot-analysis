function [FigureHandles] = Setup_Figure_Windows(Options)

FigureHandles.GradientWindow = figure(23);
%     set(FigureHandles.GradientWindow, 'Position', [1    79   560   420]);
    set(FigureHandles.GradientWindow, 'Position', [1 52 450 318]);
FigureHandles.DiagnosticWindow = figure(24);
    set(FigureHandles.DiagnosticWindow, 'Position', [430 386 424 310]);
FigureHandles.TraceWindow = figure(1);
    set(FigureHandles.TraceWindow, 'Position', [2 389 447 316]);

end