function [FigureHandles] = Create_Master_Window(Options)

    FigureHandles.MasterWindow = figure(1);
    set(FigureHandles.MasterWindow, 'Position', [2 53 1278 652]);
    set(0, 'DefaultAxesFontSize',11)

    NumPlotsX = Options.NumPlotsX;
    NumPlotsY = Options.NumPlotsY;
    

    Gap = [.04,.01];
    MarginsHeight = [.04,.04];
    MarginsWidth = [.03,.02];


    FigureHandles.SubHandles = tight_subplot(NumPlotsY, NumPlotsX, Gap, MarginsHeight, MarginsWidth);
%     
% for b = 1:NumPlots
%     axes(SubHandles(b));
%     plot(randn(10,b));
%     title(num2str(b));
% end

    

end