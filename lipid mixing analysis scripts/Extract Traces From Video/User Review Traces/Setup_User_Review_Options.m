 function [DataCounters,Options] = Setup_User_Review_Options()
    Options.Label = '-Revd';
    Options.NumPlotsX = 6;
    Options.NumPlotsY = 3;
    Options.TotalNumPlots = Options.NumPlotsX*Options.NumPlotsY;
    Options.StartingTraceNumber = 1;
    
    Options.SaveAtEachStep = 'y';
    Options.QuickModeNoCorrection = 'n';
    
    DataCounters.CurrentTraceNumber = Options.StartingTraceNumber;
    DataCounters.CurrentErrorRate = 0;
    
% To enter at prompt: PlotNumber.DesignationCode
    % DesignationCode as follows:
    % .0 = No Fusion
    % .1 = 1 Fuse
    % .12 = 1 Fuse designation is already correct, but wait time is wrong
    % .2 = Abnormal (e.g. slow) fusion
    % .9 = Hard To Classify, Ignore This One
    
    Options.AddPresetOptions = 'n';
    if strcmp(Options.AddPresetOptions, 'y')
        PresetOptionsDir = '/Users/bobrawle/Matlab/Virus Fusion Scripts/Preset Options/User Review Traces';
        [PresetOptionsFile, PresetOptionsDir] = uigetfile('*.m','Select pre-set options .m file',...
            char(PresetOptionsDir),'Multiselect', 'off');
        run(strcat(PresetOptionsDir,PresetOptionsFile));
    end
end