function [FigureHandles] = Setup_Figures(Options)

        FigureHandles.BinaryImageWindow = figure(3);
        set(FigureHandles.BinaryImageWindow,'Position',[1 -50 450 341]);
        FigureHandles.BackgroundTraceWindow = figure(4);
        set(FigureHandles.BackgroundTraceWindow,'Position',[452 -130 450 341]);
        FigureHandles.ImageWindow = figure(1);
        set(FigureHandles.ImageWindow,'Position',[6   479   451   338]);
        FigureHandles.CurrentTraceWindow = figure(2);
        set(FigureHandles.CurrentTraceWindow,'Position',[472   476   450   341]);
        
        if Options.UserGrabExampleTrace == 'y'
            FigureHandles.UserExampleTrace = figure(6);
            set(FigureHandles.UserExampleTrace,'Position',[250 -50 700 550]);
        end
end