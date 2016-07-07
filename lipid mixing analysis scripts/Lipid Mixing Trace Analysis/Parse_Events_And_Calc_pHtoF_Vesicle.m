function [StatsOfFailures,DockingData,FusionData,...
    StatsOfDesignations,AnalyzedTraceData] =...
    Parse_Events_And_Calc_pHtoF_Vesicle(StatsOfFailures,TraceRunMedian,...
    StatsOfDesignations,CurrentVirusData,FigureHandles,Options,UniversalData,...
    AnalyzedTraceData,DockingData,TraceGradData)

FusionData.Designation='failed';  % by default we say that it has failed
FusionData.FuseFrameNumbers = [];
FusionData.pHtoFusionNumFrames = [];
FusionData.pHtoFusionTime = [];
% The spike data is used later on if we are going to cross correlate each trace with
% nearby traces to remove spurious fusion events 
SpikeData.FrameNumbers = [];
SpikeData.IntensityHeight = [];
SpikeData.CorrelationData(1) = NaN; 
SpikeData.CorrelationData(2) = NaN;

[RangeToFilterPositive,RangeToFilterNegative,GradTraceRunMed,NegFilteredGradTrace,PosFilteredGradTrace,...
    FilteredSpikeTrace,RangeToFilterSpike,SpikeFrameNumbers]...
    = Unpack_Grad_Filter_Data(TraceGradData);

%First we determine if the virus was mobile, and if so, we determine
%the frame number at which it stopped
    [DockingData] = Det_If_Mobile(DockingData,UniversalData.pHDropFrameNumber,Options);
    % The docking data will always be set to no docking because we are
    % analyzing fusion to a tethered vesicle
    
    if strcmp(DockingData.IsMobile,'y')
        StopFrameIndex = find(TraceRunMedian.FrameNumbers== DockingData.StopFrameNum);
        PosFilteredGradTrace(1:StopFrameIndex) = 0;
        NegFilteredGradTrace(1:StopFrameIndex) = 0;
        TraceGradData.FilteredDiffTracePos(1:StopFrameIndex) = 0;
        FilteredSpikeTrace(1:StopFrameIndex) = 0;
        TraceGradData.FilteredDiffTraceNeg(1:StopFrameIndex) = 0;
    end
     
    %Grab the frame numbers of potential fusion events
        FuseUpFrameNumbers = TraceRunMedian.FrameNumbers(PosFilteredGradTrace);
        FuseDownFrameNumbers = TraceRunMedian.FrameNumbers(NegFilteredGradTrace);
        FuseSpikeFrameNumbers = SpikeFrameNumbers(FilteredSpikeTrace);
        NumFuseUpEvents = sum(PosFilteredGradTrace);
        NumFuseDownEvents = sum(NegFilteredGradTrace);
        NumFuseSpikeEvents = sum(FilteredSpikeTrace);
        NumFuseTotalEvents = NumFuseUpEvents;
        
    %Parse Events and Calc pHtoF Time when appropriate
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
    
    % Now we calculate the number of frames and time between the pH drop and any
    % fusion events which were recorded.
        FusionData.pHtoFusionNumFrames = FusionData.FuseFrameNumbers - UniversalData.pHDropFrameNumber;
        FusionData.pHtoFusionTime = Options.TimeInterval*(FusionData.pHtoFusionNumFrames);
        FusionData.SpikeData = SpikeData;
    
        %Now we  draw lines on the Trace window to indicate where fusion, docking,
        %etc. occurred
        [FigureHandles] = Draw_Lines_on_Plot(FigureHandles,DockingData,FusionData,UniversalData);
    
        % Now we compile all the data for this analyzed trace into the
        % analyzed trace data structure (this is what will ultimately be
        % saved)
            ChangedByUser = 'Not analyzed';
            Write_Designation_On_Plot(FigureHandles.TraceWindow,FusionData.Designation);

            SaveOption = 'DontWritetoDiskyet';
            [AnalyzedTraceData,StatsOfDesignations] = Compile_Analyzed_Trace_Data(UniversalData,FusionData,...
            CurrentVirusData,StatsOfDesignations,DockingData,AnalyzedTraceData,ChangedByUser,Options,...
            TraceGradData,SaveOption);
end

function [RangeToFilterPositive,RangeToFilterNegative,GradTraceRunMed,NegFilteredGradTrace,PosFilteredGradTrace,...
    FilteredSpikeTrace,RangeToFilterSpike,SpikeFrameNumbers]...
    = Unpack_Grad_Filter_Data(TraceGradData)

    RangeToFilterSpike = TraceGradData.RangeToFilterSpike;
    FilteredSpikeTrace = TraceGradData.FilteredSpikeTrace;
    NegFilteredGradTrace = TraceGradData.NegFilteredGradTrace;
    PosFilteredGradTrace = TraceGradData.PosFilteredGradTrace;
    GradTraceRunMed = TraceGradData.GradTraceRunMed;
    RangeToFilterPositive = TraceGradData.RangeToFilterPositive;
    RangeToFilterNegative = TraceGradData.RangeToFilterNegative;
    SpikeFrameNumbers = TraceGradData.SpikeFrameNumbers;
end