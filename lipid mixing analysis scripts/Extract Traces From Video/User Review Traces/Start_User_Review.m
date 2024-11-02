function Start_User_Review(varargin)

close all

% Identify the path to the data you wish to analyze
    [DataFilename,DefaultPathname] = Load_Data(varargin);
        % Nested function

    DataFilePath = strcat(DefaultPathname,DataFilename);

% Setup options
    [DataCounters,Options] = Setup_User_Review_Options();
    
% Extract data that we will need
    InputData = open(DataFilePath);
    DataToSave = InputData.DataToSave;
    UniversalData = InputData.DataToSave.OtherDataToSave.UniversalData;
    PreviousAnalysisData = InputData.DataToSave.CombinedAnalyzedTraceData;
    NumTraces = length(PreviousAnalysisData);
    
% Define variables we will need
    CorrectedAnalysisData = PreviousAnalysisData;
    NumTracesToReview = NumTraces - Options.StartingTraceNumber - 1;
    NumReviewRounds = ceil(NumTracesToReview/Options.TotalNumPlots);
    TraceCounter = Options.StartingTraceNumber;
    ErrorCounter = 0;

% Create master window with subplots
    [FigureHandles] = Create_Master_Window(Options);

% Review plots round by round
    for b = 1:NumReviewRounds

        disp(strcat('Round-', num2str(b),'-of-', num2str(NumReviewRounds)))
        
        if b ~= NumReviewRounds
            CurrentTraceRange = TraceCounter:TraceCounter + Options.TotalNumPlots - 1;
            TraceCounter = max(CurrentTraceRange) +1;
        else
            CurrentTraceRange = TraceCounter:NumTraces;
            for d = length(CurrentTraceRange) + 1: Options.TotalNumPlots
                set(FigureHandles.MasterWindow,'CurrentAxes',FigureHandles.SubHandles(d));
                cla
            end
        end
        
        % Load and plot the current round
        PlotCounter = 1;
        for i = CurrentTraceRange
            
            CurrentTraceNumber = i;
            CurrentVirusData = PreviousAnalysisData(i);
            CurrentTraceBackSub = CurrentVirusData.Trace_BackSub;
            
            % Correct focus problems
            [CurrentTraceBackSub] = Correct_Focus_Problems(CurrentTraceBackSub,UniversalData);
            
            
            
            [FigureHandles] = Plot_Current_Trace(FigureHandles,CurrentVirusData,UniversalData,...
                CurrentTraceBackSub,PlotCounter,CurrentTraceNumber);
            
            CorrectedAnalysisData(i).ChangedByUser = 'Reviewed By User';
            PlotCounter = PlotCounter +1;
        end

        % Ask user if we need to change any of the designations on the current round of plots
        RerunThisRound = 'y';
        while RerunThisRound =='y'
            
            Prompts = {strcat(num2str(b),'/', num2str(NumReviewRounds),'; List IncorrectNumber.DesigCode')};
            DefaultInputs = {'No Correction Needed'};
            Heading = 'Type q to quit';
            UserAnswer = inputdlg(Prompts,Heading, 1, DefaultInputs, 'on');

            if isempty(UserAnswer)
                % There has been an error, re-run the last round to avoid crash
                RerunThisRound = 'y';
                
            elseif strcmp(UserAnswer{1,1},'q')
                disp('You Chose To Quit')
                ThisWillCauseError
            
            elseif strcmp(UserAnswer{1,1},'No Correction Needed')
                % Everything is correct, move to next round
                RerunThisRound = 'n';
                
            else
                
                % Extract User Inputs
                IncorrectPlotIndices = str2num(UserAnswer{1,1}); 
                
                if isvector(IncorrectPlotIndices)
                    % User has indicated that we need to correct some designations                
                    
                    [RerunThisRound, CorrectedAnalysisData, ErrorCounter] = Correct_Designations(IncorrectPlotIndices,...
                        PreviousAnalysisData,CurrentTraceRange,CorrectedAnalysisData,ErrorCounter,Options);
                    
                else
                    % There has been an error, re-run the last round to avoid crash
                    RerunThisRound = 'y';
                end
            end
        end
        
        DataCounters.CurrentErrorCount = ErrorCounter;
        DataCounters.CurrentTraceNumber = CurrentTraceNumber;
        DataCounters.CurrentErrorRate = ErrorCounter/CurrentTraceNumber;
        
        DataCounters
        
        if  strcmp(Options.SaveAtEachStep,'y')
            DataToSave.DataCounters = DataCounters;
            Save_Data_At_Each_Step(DataFilename,DefaultPathname,DataToSave,CorrectedAnalysisData,Options)
        end
    end
      
    if ~strcmp(Options.SaveAtEachStep,'y')
        DataToSave.DataCounters = DataCounters;
        Save_Data_At_Each_Step(DataFilename,DefaultPathname,DataToSave,CorrectedAnalysisData,Options)
    end
    
	disp('Thank You. Come Again.')
    
    ThisWillCauseFunctiontoThrowError

end

function [DataFilenames,DefaultPathname] = Load_Data(varargin)
    if length(varargin) == 1
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed',...
            char(varargin{1}),'Multiselect', 'on');
    else
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed', 'Multiselect', 'on');
    end
end

function Save_Data_At_Each_Step(DataFilename,DefaultPathname,DataToSave,CorrectedAnalysisData,Options)

SaveDataFolder = DefaultPathname;

SaveDataFolder = strcat(SaveDataFolder,'/AnalysisReviewed/');
if exist(SaveDataFolder,'dir') == 0
    mkdir(SaveDataFolder);
end

DataFilenameWOExt = DataFilename(1:end-4);
    if ~isempty(CorrectedAnalysisData)
        DataToSave.ReviewOptions = Options;
        DataToSave.CombinedAnalyzedTraceData = CorrectedAnalysisData;
        save(strcat(SaveDataFolder,DataFilenameWOExt,Options.Label,'.mat'),'DataToSave');
    end
end