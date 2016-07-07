function [AnalyzedTraceData,OtherDataToSave,StatsOfFailures,StatsOfDesignations] =...
        Analyze_Current_Data_Set(CurrDataFilePath,Options)

[FigureHandles] = Setup_Figure_Windows(Options);

%Set up stats of failures and designations
    [StatsOfFailures,StatsOfDesignations] = Set_Up_Stats();

    InputData = open(CurrDataFilePath);
    if isfield(InputData,'VirusDataToSave')
        InputDataType = 'All Traces Saved';
        InputTraceData = InputData.VirusDataToSave;
        OtherImportedData = InputData.OtherDataToSave;
    end
    
    
    AnalyzedTraceData = [];

    %Determine the frame numbers in which the pH drop occurred, when focusing issues
    %happened, and when all of the virus particles stopped moving. Some of
    %these values may have been predefined before this point.
    [FrameAllVirusStoppedBy,PHdropFrameNum,focusframenumbers,focusproblems] =...
        Determine_pH_Focus_Stop_FrameNumbers(OtherImportedData,InputTraceData,Options);

        UniversalData.FrameAllVirusStoppedBy = FrameAllVirusStoppedBy;
        UniversalData.pHDropFrameNumber = PHdropFrameNum;
        UniversalData.FocusFrameNumbers = [focusframenumbers' Options.AdditionalFocusFrameNumbers];
        UniversalData.FocusProblems = focusproblems;
        
    UniversalData.NumTraces = length(InputTraceData);

    for i = Options.StartingTraceNumber:UniversalData.NumTraces
        UniversalData.TraceNumber = i;
        CurrentVirusData = InputTraceData(i);
        
        if strcmp(CurrentVirusData.IsVirusGood,'y') || strcmp(InputDataType,'Only Good Saved')
        
            CurrTrace = InputTraceData(i).Trace_BackSub;
            CurrentVirusData.focusframenumbers = UniversalData.FocusFrameNumbers;
                % This overwrites the focus frame numbers if the user specified
                % different ones than were prerecorded
                
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

            %Deal with focus problems.
                [CurrTrace] = Correct_Focus_Problems(CurrTrace,UniversalData);

            %Can limit frames if we need to
                if isnan(Options.FrametoEndAnalysis)
                    Options.FrametoEndAnalysis = length(CurrTrace);
                end

                CurrTraceCropped.Trace = CurrTrace(Options.FrameToStartAnalysis:Options.FrametoEndAnalysis);
                CurrTraceCropped.FrameNumbers = Options.FrameToStartAnalysis:length(CurrTraceCropped.Trace)+Options.FrameToStartAnalysis-1;

            %----------------Define/Apply Gradients Filters--------------------
            [TraceRunMedian,FigureHandles] = ...
                Run_Med_And_Plot(CurrTraceCropped,FigureHandles,UniversalData,Options);

            %Diagnostic
        %     if strcmp(UniversalData.FocusProblems,'y')
        %         set(0,'CurrentFigure',FigureHandles.TraceWindow)
        %         hold on
        %         plot(CurrTraceCropped.FrameNumbers,OldCurrTrace,'g-')
        %         hold off
        %     end

            %Define and apply the gradient filters
             [TraceGradData,DockingData] =...
            Define_and_Apply_Gradient_Filters(FigureHandles,UniversalData,Options,TraceRunMedian,CurrTraceCropped);

        %         set(0,'CurrentFigure',FigureHandles.DiagnosticWindow)
        % %         hold on
        %         DiagnosticData =TraceRunMedian.Trace(UniversalData.FrameAllVirusStoppedBy:end);
        %         [ Counts, Bins] = hist(DiagnosticData,20);
        %         barh(Bins,Counts)
        % %         hold off
            %------------------------------------------------------------------

            %-------------------------Calc pHtoF Time---------------------------

            %We parse up the current trace into number of 
            %of fusion events and calculate a rough estimate of the pHtoF time 
            %when appropriate.
            if strcmp(Options.TypeofFusionData, 'TetheredVesicle')
                [StatsOfFailures,DockingData,FusionData,...
                StatsOfDesignations,AnalyzedTraceData] =...
                Parse_Events_And_Calc_pHtoF_Vesicle(StatsOfFailures,TraceRunMedian,...
                StatsOfDesignations,CurrentVirusData,FigureHandles,Options,UniversalData,...
                AnalyzedTraceData,DockingData,TraceGradData);
            end

            %------------------------------------------------------------------
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