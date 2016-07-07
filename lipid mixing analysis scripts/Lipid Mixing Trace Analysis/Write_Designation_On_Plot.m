function Write_Designation_On_Plot(TraceWindow, TextToDisplay)
    set(0,'CurrentFigure',TraceWindow)
    hold on
    X = xlim;
    Y = ylim;
        XMidpoint = (X(2)-X(1))/2 + X(1);
        if strcmp(TextToDisplay,'Slow')
            YMidpoint = (Y(2)-Y(1))/4 + Y(1);
        else 
            YMidpoint = (Y(2)-Y(1))/2 + Y(1);
        end
    text(XMidpoint,YMidpoint,TextToDisplay,'FontSize',20,'EdgeColor','k',...
        'HorizontalAlignment','center');
    drawnow
    hold off
end