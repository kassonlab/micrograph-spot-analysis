function [DefaultPathname] = Start_Trace_Analysis_Program(varargin)
% - - - - - - - - - - - - - - - - - - - - -

% Input:
% Start_Trace_Analysis_Program(), in this case the user navigates to the
%       .mat output file from the Extract Traces From Video program.
%   OR
% Start_Trace_Analysis_Program(DefaultPath), where DefaultPath is the directory to which 
%       the user will be automatically directed to find the 
%       .mat output file from the Extract Traces From Video program.

% Output:
% A .mat file is created which saves all of the variables in the current 
% workspace. This file is saved in a new folder created in the parent 
% directory where the input file came from. The information about the 
% waiting time for each lipid mixing event, as well as the designation 
% of each trace, is contained in the DataToSave.CombinedAnalyzedTraceData 
% structure, as defined in the Compile_Analyzed_Trace_Data.m file.

% Note: This program has been designed to process many sets of data sequentially,
% but it has been tested with individual sets, so keep that 
% in mind if you choose to process many sets at once.

% By Bob Rawle, Kasson Lab, University of Virginia, 2016
% Published online in conjunction with:
% Rawle et al., Disentangling Viral Membrane Fusion from Receptor Binding 
% Using Synthetic DNA-Lipid Conjugates, Biophysical Journal (2016) 
% http://dx.doi.org/10.1016/j.bpj.2016.05.048
% - - - - - - - - - - - - - - - - - - - - -

    %Debugging
    dbstop in Start_Trace_Analysis_Program at 88
%     dbstop in Parse_Events_And_Calc_pHtoF_Vesicle at 292
%     dbstop in Setup_Figure_Windows at 3
    dbstop in Analyze_Current_Data_Set at 4

    [DataFilenames,DefaultPathname] = Load_Data(varargin);
%         Nested function

%     Define the options
    [Options] = Setup_Options(DefaultPathname);
   
        RestartCount = []; CombinedAnalyzedTraceData = [];
        disp(' '); disp(' '); disp (' ');
        
    if iscell(DataFilenames) %This lets us know if there is more than one file
        NumberOfFiles = length(DataFilenames);
    else
        NumberOfFiles = 1;
    end

    for i = 1:NumberOfFiles
        if NumberOfFiles > 1
           CurrDataFileName = DataFilenames{1,i};
        else
           CurrDataFileName = DataFilenames;
        end
        
        CurrDataFilePath = strcat(DefaultPathname,CurrDataFileName);

        % Call the analysis function to analyze the data from the current
        % set
        [AnalyzedTraceData,OtherDataToSave,StatsOfFailures,StatsOfDesignations] =...
        Analyze_Current_Data_Set(CurrDataFilePath,Options);

        Save_Data_At_Each_Step(AnalyzedTraceData,OtherDataToSave,DefaultPathname,Options.Label,Options)

        
        %To combine the data, we have to deal with empty structures,
        %which can create problems.  So we deal with it and then
        %combine the current data with the previous iterations
        [CombinedAnalyzedTraceData,RestartCount]= ...
            Deal_With_Empty_Recorded_Data(i,AnalyzedTraceData,CombinedAnalyzedTraceData,RestartCount);

        disp(strcat('-----------------File_', num2str(i),'_of_',num2str(NumberOfFiles),'-----------------'))
        disp(' ')
        disp(strcat('Filename: ', CurrDataFileName))
        StatsOfFailures
        StatsOfDesignations
        disp('---------------------------------------------')
        disp(' ')
        disp(' ')

    end
    
    %Save combined data in big structure
        disp('Analysis completed for all files.');
        Stop_ThisWillCauseErrortoPreventMatlabfromCrashing

        
disp('Thank you.  Come again.')

        
end

function Save_Data_At_Each_Step(AnalyzedTraceData,OtherDataToSave,DefaultPathname,Label,Options)

DataToSave.OtherDataToSave = OtherDataToSave;

SaveDataFolder = DefaultPathname;
SaveDataFolder = strcat(SaveDataFolder,'/Analysis/');
if exist(SaveDataFolder,'dir') == 0
    mkdir(SaveDataFolder);
end

    if ~isempty(AnalyzedTraceData)
        DataToSave.CombinedAnalyzedTraceData = AnalyzedTraceData;
        save(strcat(SaveDataFolder,Label,'.mat'),'DataToSave');
    end
end


function [CombinedAnalyzedTraceData,RestartCount]= ...
                Deal_With_Empty_Recorded_Data(i,AnalyzedTraceData,CombinedAnalyzedTraceData,RestartCount)

    %Compile the data which will be saved (there are lots of
    %complicated if statements here just to deal with the times
    %that there doesn't happen to be any events recorded in a given file).
    if i == 1
        if ~isempty(AnalyzedTraceData)
            CombinedAnalyzedTraceData = AnalyzedTraceData; %This is a structure
            RestartCount = 'n';
        else
            RestartCount = 'y';
        end
        
    else
        if RestartCount == 'y'
            if ~isempty(AnalyzedTraceData)
                CombinedAnalyzedTraceData = AnalyzedTraceData; %This is a structure
                RestartCount = 'n';
            else
                RestartCount = 'y';
            end
        elseif ~isempty(AnalyzedTraceData)
            StartIdx = length(CombinedAnalyzedTraceData) + 1;
            EndIdx = StartIdx + length(AnalyzedTraceData)-1;
            CombinedAnalyzedTraceData(StartIdx:EndIdx) = AnalyzedTraceData;
        end
    end
end

function [DataFilenames,DefaultPathname] = Load_Data(varargin)
%First, we load the .mat data files.
        if length(varargin) == 1
            [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed',...
                char(varargin{1}),'Multiselect', 'on');
        else
            [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed', 'Multiselect', 'on');
        end
end