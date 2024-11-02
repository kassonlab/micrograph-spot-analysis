
close all
FigureHandles.MasterWindow = figure(1);
set(FigureHandles.MasterWindow, 'Position', [2 53 1278 652]);
set(FigureHandles.MasterWindow, 'Border','tight');
set(0, 'DefaultAxesFontSize',11)

NumPlotsX = 5;
NumPlotsY = 3;
NumPlots = NumPlotsX*NumPlotsY;

for b = 1:NumPlots
    SubHandles(b).Axis = subplot(NumPlotsY,NumPlotsX,b);
    title(num2str(b));
    axis tight
end
