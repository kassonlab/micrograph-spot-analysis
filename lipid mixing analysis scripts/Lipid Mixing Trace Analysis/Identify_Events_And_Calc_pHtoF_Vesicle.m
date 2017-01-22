function [StatsOfFailures,DockingData,FusionData,...
    StatsOfDesignations,AnalyzedTraceData] =...
    Identify_Events_And_Calc_pHtoF_Vesicle(StatsOfFailures,TraceRunMedian,...
    StatsOfDesignations,CurrentVirusData,FigureHandles,Options,UniversalData,...
    AnalyzedTraceData,DockingData,TraceGradData)

% Pre-set some default values, these will be changed as we go on
    FusionData.Designation='failed';  
    FusionData.FuseFrameNumbers = [];
    FusionData.pHtoFusionNumFrames = [];
    FusionData.pHtoFusionTime = [];
    % The spike data is used later on if we are going to cross correlate each trace with
    % nearby traces to remove spurious fusion events - only relevant for
    % SLB fusion data
    SpikeData.FrameNumbers = [];
    SpikeData.IntensityHeight = [];
    SpikeData.CorrelationData(1) = NaN; 
    SpikeData.CorrelationData(2) = NaN;

    % Unpack the data from the various threshold tests which will be used 
    % to identify fusion events below
    [RangeToFilterPositive,RangeToFilterNegative,GradTraceRunMed,NegFilteredGradTrace,PosFilteredGradTrace,...
    FilteredSpikeTrace,RangeToFilterSpike,SpikeFrameNumbers]...
    = Unpack_Thresh_Filter_Data(TraceGradData);

% Now we determine if the virus was mobile, and if so, we determine
% the frame number at which it stopped
    [DockingData] = Det_If_Mobile(DockingData,UniversalData.pHDropFrameNumber,Options);
        % Note: in the case of fusion to tethered vesicle data, this will 
        % always be identified as 'no docking'
    
    % If the virus was determined to be mobile, then we make sure that spurious 
    % fusion events aren't identified before the virus docks
    % Note: Again, this is not relevant to fusion to tethered vesicle data
    if strcmp(DockingData.IsMobile,'y')
        StopFrameIndex = find(TraceRunMedian.FrameNumbers== DockingData.StopFrameNum);
        PosFilteredGradTrace(1:StopFrameIndex) = 0;
        NegFilteredGradTrace(1:StopFrameIndex) = 0;
        TraceGradData.FilteredDiffTracePos(1:StopFrameIndex) = 0;
        FilteredSpikeTrace(1:StopFrameIndex) = 0;
        TraceGradData.FilteredDiffTraceNeg(1:StopFrameIndex) = 0;
    end
     
    %Now we identify the frame numbers of potential fusion events
        FuseUpFrameNumbers = TraceRunMedian.FrameNumbers(PosFilteredGradTrace);
        FuseDownFrameNumbers = TraceRunMedian.FrameNumbers(NegFilteredGradTrace);
        FuseSpikeFrameNumbers = SpikeFrameNumbers(FilteredSpikeTrace);
        NumFuseUpEvents = sum(PosFilteredGradTrace);
        NumFuseDownEvents = sum(NegFilteredGradTrace);
        NumFuseSpikeEvents = sum(FilteredSpikeTrace);
        NumFuseTotalEvents = NumFuseUpEvents;

    % Filter out spike events which might be erroneously identified as 
    % fusion events (this can happen especially if your running median 
    % length is set to zero)
    if NumFuseSpikeEvents > 0 && NumFuseUpEvents > 0
        CutoffDistance = 2;
            % CutoffDistance = How far apart (in frame numbers) must the spike and 
            % fusion events be in order to not be considered the same event
        for q = 1:NumFuseSpikeEvents
            HowFarfromSpike = abs(FuseUpFrameNumbers - FuseSpikeFrameNumbers(q));
            FrameNumbersNottheSame = HowFarfromSpike > CutoffDistance;
            FuseUpFrameNumbers = FuseUpFrameNumbers(FrameNumbersNottheSame);
        end
        NumFuseUpEvents = length(FuseUpFrameNumbers);
        NumFuseTotalEvents = NumFuseUpEvents;
    end
    
    % Identify num of fusion events
    if NumFuseTotalEvents == 0
        FusionData.Designation = 'No Fusion';
    elseif NumFuseTotalEvents == 1 && NumFuseUpEvents == 1
        FusionData.Designation = '1 Fuse';
        FusionData.FuseFrameNumbers = FuseUpFrameNumbers(1);

    elseif NumFuseTotalEvents == 2 && NumFuseUpEvents == 2
        FusionData.Designation = '2 Fuse';
        FusionData.FuseFrameNumbers = FuseUpFrameNumbers(1:2);
         
    else 
        %then there are too many fuse events (i.e. >2)
            Reason_Failed = 'Too Many Fuse Event';
            Cross_Out_Plot(FigureHandles.TraceWindow,Reason_Failed)
                StatsOfFailures.TooManyFuseEvent = StatsOfFailures.TooManyFuseEvent + 1;
    end
    
    % Now we determine if there is a slow fusion event (where the intensity
    % changes considerably, but not sharply as we would expect for a normal
    % fusion event)
%     DetectionOption = 'Usual Trace Analysis';
    DetectionOption = 'Cluster Analysis';
    [TraceGradData,FusionData] = Is_Slow_Fusion(TraceGradData,FusionData,...
        DockingData,DetectionOption,Options);
    
    % Calculate the number of frames and time between the pH drop and any
    % fusion events which were recorded.
        FusionData.pHtoFusionNumFrames = FusionData.FuseFrameNumbers - UniversalData.pHDropFrameNumber;
        FusionData.pHtoFusionTime = Options.TimeInterval*(FusionData.pHtoFusionNumFrames);
        FusionData.SpikeData = SpikeData;
    
        % Draw lines on the trace window to indicate where fusion, docking,
        % etc. occurred, and write designation on plot
        [FigureHandles] = Draw_Lines_on_Plot(FigureHandles,DockingData,FusionData,UniversalData);
        Write_Designation_On_Plot(FigureHandles.TraceWindow,FusionData.Designation);
    
        % Now we compile all the data for this analyzed trace into the
        % analyzed trace data structure (this is what will ultimately be
        % saved)

            % Debug
            if strcmp(FusionData.Designation,'1 Fuse')
                   debugstop = 1;
            end
            
            ChangedByUser = 'Not analyzed';
            SaveOption = 'DontWritetoDiskyet';
            [AnalyzedTraceData,StatsOfDesignations] = Compile_Analyzed_Trace_Data(UniversalData,FusionData,...
            CurrentVirusData,StatsOfDesignations,DockingData,AnalyzedTraceData,ChangedByUser,Options,...
            TraceGradData,SaveOption);
end

function [RangeToFilterPositive,RangeToFilterNegative,GradTraceRunMed,NegFilteredGradTrace,PosFilteredGradTrace,...
    FilteredSpikeTrace,RangeToFilterSpike,SpikeFrameNumbers]...
    = Unpack_Thresh_Filter_Data(TraceGradData)

    RangeToFilterSpike = TraceGradData.RangeToFilterSpike;
    FilteredSpikeTrace = TraceGradData.FilteredSpikeTrace;
    NegFilteredGradTrace = TraceGradData.NegFilteredGradTrace;
    PosFilteredGradTrace = TraceGradData.PosFilteredGradTrace;
    GradTraceRunMed = TraceGradData.GradTraceRunMed;
    RangeToFilterPositive = TraceGradData.RangeToFilterPositive;
    RangeToFilterNegative = TraceGradData.RangeToFilterNegative;
    SpikeFrameNumbers = TraceGradData.SpikeFrameNumbers;
end