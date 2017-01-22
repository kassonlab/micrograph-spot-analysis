function [AnalyzedTraceData,OtherDataToSave,StatsOfFailures,StatsOfDesignations] =...
        Analyze_Current_Data_Set(CurrDataFilePath,Options)

% ------------------- Set up before analysis -------------------
    [FigureHandles] = Setup_Figure_Windows(Options);

    % Set up stats of failures and designations (for tracking purposes)
    [StatsOfFailures,StatsOfDesignations] = Set_Up_Stats();
        % Nested function
        
    % Load current data set from the file path
    InputData = open(CurrDataFilePath);
    if isfield(InputData,'VirusDataToSave')
        InputDataType = 'All Traces Saved';
        InputTraceData = InputData.VirusDataToSave;
        OtherImportedData = InputData.OtherDataToSave;
    end
    
    % Set up output analysis structure variable
    AnalyzedTraceData = [];

    %Determine the frame numbers in which the pH drop occurred, when focusing issues
    %happened, and when all of the virus particles stopped moving. Some/all of
    %these values may have been predefined before this point.
    [FrameAllVirusStoppedBy,PHdropFrameNum,focusframenumbers,focusproblems] =...
        Determine_pH_Focus_Stop_FrameNumbers(OtherImportedData,InputTraceData,Options);

        UniversalData.FrameAllVirusStoppedBy = FrameAllVirusStoppedBy;
        UniversalData.pHDropFrameNumber = PHdropFrameNum;
        UniversalData.FocusFrameNumbers = focusframenumbers;
        UniversalData.FocusProblems = focusproblems;
        
    UniversalData.NumTraces = length(InputTraceData);

    % ------------------- Start analysis for each trace -------------------
    for i = Options.StartingTraceNumber:UniversalData.NumTraces
        UniversalData.TraceNumber = i;
        CurrentVirusData = InputTraceData(i);
        
        % Only the traces of 'good' viruses are analyzed and included. The 
        % rest are excluded from analysis.
        if strcmp(CurrentVirusData.IsVirusGood,'y') || strcmp(InputDataType,'Only Good Saved')
        
            CurrTrace = InputTraceData(i).Trace_BackSub;
%             CurrTrace = InputTraceData(i).Trace;

            % ------------------- Pre-processing of the current trace -------------------
                % Overwrite the focus frame numbers if the user specified
                % different ones than were prerecorded in the Extract Traces From Video scripts.
                CurrentVirusData.focusframenumbers = UniversalData.FocusFrameNumbers;
                
                % Deal with legacy data format
                if strcmp(InputDataType,'Only Good Saved')
                    BoxAroundVirus = InputTraceData(i).BoxAroundSUV;
                elseif strcmp(InputDataType,'All Traces Saved')
                    BoxAroundVirus = InputTraceData(i).BoxAroundVirus;
                end
            
                CurrentVirusData.BoxCoords = [BoxAroundVirus.Right, BoxAroundVirus.Bottom;
                    BoxAroundVirus.Left, BoxAroundVirus.Bottom;
                    BoxAroundVirus.Left, BoxAroundVirus.Top;
                    BoxAroundVirus.Right, BoxAroundVirus.Top;
                    BoxAroundVirus.Right, BoxAroundVirus.Bottom];

                % Correct focus problems
                [CurrTrace] = Correct_Focus_Problems(CurrTrace,UniversalData);

                % Can limit frames to analyze if we need to
                if isnan(Options.FrametoEndAnalysis)
                    Options.FrametoEndAnalysis = length(CurrTrace);
                end
                CurrTraceCropped.Trace = CurrTrace(Options.FrameToStartAnalysis:Options.FrametoEndAnalysis);
                CurrTraceCropped.FrameNumbers = Options.FrameToStartAnalysis:length(CurrTraceCropped.Trace)+Options.FrameToStartAnalysis-1;

                % Calculate the running median of the trace (used for much of the analysis) and plot the trace
                [TraceRunMedian,FigureHandles] = ...
                    Run_Med_And_Plot(CurrTraceCropped,FigureHandles,UniversalData,Options);
                
            %----------------Determine/Apply Filters for the Various Tests --------------------
            
            % Determine and apply the threshold filters
             [TraceGradData,DockingData] =...
            Det_and_Apply_Threshold_Filters(FigureHandles,UniversalData,Options,TraceRunMedian,CurrTraceCropped);

        %         set(0,'CurrentFigure',FigureHandles.DiagnosticWindow)
        % %         hold on
        %         DiagnosticData =TraceRunMedian.Trace(UniversalData.FrameAllVirusStoppedBy:end);
        %         [ Counts, Bins] = hist(DiagnosticData,20);
        %         barh(Bins,Counts)
        % %         hold off

            % --------------Identify fusion events and calc pHtoF Time---------------

            % We use the information from the various tests to identify fusion 
            % events and then calculate the waiting time between pH drop and 
            % fusion (lipid-mixing). This information is then saved for the 
            % current trace to the combined data structure (AnalyzedTraceData).
            if strcmp(Options.TypeofFusionData, 'TetheredVesicle')
                [StatsOfFailures,DockingData,FusionData,...
                StatsOfDesignations,AnalyzedTraceData] =...
                Identify_Events_And_Calc_pHtoF_Vesicle(StatsOfFailures,TraceRunMedian,...
                StatsOfDesignations,CurrentVirusData,FigureHandles,Options,UniversalData,...
                AnalyzedTraceData,DockingData,TraceGradData);
            end

            %------------------------------------------------------------------
            
        % Traces from viruses identified as 'bad' are ignored, and the number 
        % of bad viruses is tabulated in the stats structure.
        elseif strcmp(CurrentVirusData.IsVirusGood,'n')
            StatsOfFailures.BadVirusRegion = StatsOfFailures.BadVirusRegion +1;
        end
    end
    
    % Make sure we record other useful data not associated with each
    % individual trace
    OtherDataToSave.Options = Options;
    OtherDataToSave.UniversalData = UniversalData;
end

function [StatsOfFailures,StatsOfDesignations]= Set_Up_Stats()
    StatsOfFailures.TooManyDock = 0;
    StatsOfFailures.TooManyFuseEvent = 0;
%     StatsOfFailures.FastFuseBeg = 0;
%     StatsOfFailures.PerDecreaseTooSmall = 0;
%     StatsOfFailures.FastFuseEnd = 0;
    StatsOfFailures.pHtoFLess0 = 0;
    StatsOfFailures.StrangeNoFuseEvent = 0;
%     StatsOfFailures.NoDockEvent = 0;
    StatsOfFailures.UserRejected = 0;
%     StatsOfFailures.FutureDock = 0;
    StatsOfFailures.NoFuseEvent = 0;
    StatsOfFailures.WonkyFusionEvent = 0;
    StatsOfFailures.BadVirusRegion = 0;

    StatsOfDesignations.Other = 0;
    StatsOfDesignations.NoFuse = 0;
    StatsOfDesignations.Fuse1 = 0;
    StatsOfDesignations.Fuse2 = 0;
end