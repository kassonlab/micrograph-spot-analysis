function [TraceGradData,DockingData] =...
    Define_and_Apply_Gradient_Filters(FigureHandles,UniversalData,Options,TraceRunMedian,CurrTraceCropped)

    ProbableDockingEvent = 'n';
        % The trigger for a probable docking event is set to y only if the type
        % of fusion data is correct and if the event meets certain criteria
    MaxIntensity = max(TraceRunMedian.Trace);
    
    %Calc the gradient of the running median
        GradTraceRunMed = gradient(TraceRunMedian.Trace);
        
    %Set all grad trace before pH drop event to zero (don't see events
    %before pH drop).
        PHdropFrameNumIdx = find(TraceRunMedian.FrameNumbers==UniversalData.pHDropFrameNumber);
        GradTraceRunMed(1:PHdropFrameNumIdx) = 0;
        MaxGrad = max(GradTraceRunMed);    

%Now we define gradient filters

    %%%%Positive Filter%%%%%%%%
    if strcmp(Options.TypeofFusionData, 'TetheredVesicle')
        STDFilterFactor = 5;
        STDFilterFactorDifference = 14;
        PositiveGradFilter = STDFilterFactor*std(GradTraceRunMed);
%         PositiveGradFilter = min([250, Max_Grad-max(.3*Max_Grad,50)]); %prev value = 150
        RangeToFilterPositive = 10;
        RangeToFilterNegative = 10;
%         VirusStopbyIndex = find(TraceRunMedian.FrameNumbers==UniversalData.FrameAllVirusStoppedBy);
        VirusStopbyIndex = find(TraceRunMedian.FrameNumbers==UniversalData.pHDropFrameNumber + 5);
    elseif strcmp(Options.TypeofFusionData, 'SLBSelfQuench')
        % This is a special case in which we need to decide whether there
        % is a docking event at the beginning. If so, then we only look for
        % fusion events in the rest of the trace. To detect whether there
        % is a docking event, we will look for a very large intensity jump
        % near the beginning of the trace which remains high for the
        % remainder
        STDFilterFactor = 4.5;
        STDFilterFactorDifference = 9;
        VirusStopbyIndex = find(TraceRunMedian.FrameNumbers==UniversalData.FrameAllVirusStoppedBy);
        TruncatedGradTraceDock = GradTraceRunMed(PHdropFrameNumIdx:VirusStopbyIndex);
        FrameNumbersDockTrace =TraceRunMedian.FrameNumbers(PHdropFrameNumIdx:VirusStopbyIndex);
        MaxGradientBeforeStop = max(TruncatedGradTraceDock);
        EarlyMedianIntensity = median(TraceRunMedian.Trace(1:UniversalData.pHDropFrameNumber));
        RangeToFilterDock = 5;
        DockGradientFilter = MaxGradientBeforeStop - 0.5*MaxGradientBeforeStop;
        DockFilteredRunMed =  TruncatedGradTraceDock > DockGradientFilter;
        
        if MaxGradientBeforeStop > 0.15*MaxIntensity && EarlyMedianIntensity < 1500
            % Probably There Is a Docking Event
            ProbableDockingEvent = 'y';    
            PositiveGradFilter = STDFilterFactor*std(GradTraceRunMed(VirusStopbyIndex:end));
        else
            PositiveGradFilter = STDFilterFactor*std(GradTraceRunMed(PHdropFrameNumIdx:end));
        end
        RangeToFilterPositive = 5;
        RangeToFilterNegative = 5;
        
    end
%         if PositiveGradFilter < 0 %To account for times when pos gradient is very small.
%             PositiveGradFilter = 0;
%         end

    %%%%Negative Filter%%%%%%%%%
        
        NegativeGradFilter = -STDFilterFactor*std(GradTraceRunMed(VirusStopbyIndex:end));
    %Apply gradient filters
        NegFilteredGradTrace = GradTraceRunMed < NegativeGradFilter;
        PosFilteredGradTrace = GradTraceRunMed > PositiveGradFilter;
            
    % Filter adjacent gradient results that are higher than the cutoff (a
    % single event can sometimes be above the cutoff for several frames on
    % either side)
        FilterOption='First Only';
        [NegFilteredGradTrace] = Filter_Adjacent_Gradient_Results(RangeToFilterNegative,NegFilteredGradTrace,FilterOption,GradTraceRunMed);
        [PosFilteredGradTrace] = Filter_Adjacent_Gradient_Results(RangeToFilterPositive,PosFilteredGradTrace,FilterOption,GradTraceRunMed);
        if strcmp(ProbableDockingEvent, 'y')
            FilterOption='Max Then First';
            [DockFilteredRunMed] = Filter_Adjacent_Gradient_Results(RangeToFilterDock,DockFilteredRunMed,FilterOption,TruncatedGradTraceDock);
        end
        
    % Calculate difference trace (difference between current point and the
    % nth previous point). This will be used to identify either a slow
    % increase in fluorescence or a slow decrease in fluorescence which may
    % indicate an unusual fusion event.
        NumberPointsBack = 25;
        DifferenceTrace = zeros(length(TraceRunMedian.Trace),1);
        DifferenceTrace(VirusStopbyIndex+NumberPointsBack+1:end) =  TraceRunMedian.Trace(VirusStopbyIndex+NumberPointsBack+1:end) -...
            TraceRunMedian.Trace(VirusStopbyIndex+1:end-NumberPointsBack);
    % Filter the difference trace in both the positive and negative
    % directions.
        DifferenceFilterPos = STDFilterFactorDifference* std(GradTraceRunMed(VirusStopbyIndex:end));
        DifferenceFilterNeg = -STDFilterFactorDifference* std(GradTraceRunMed(VirusStopbyIndex:end));

%         DifferenceFilterPos = STDFilterFactorDifference* std(DifferenceTrace(VirusStopbyIndex:end));
%         DifferenceFilterNeg = -STDFilterFactorDifference* std(DifferenceTrace(VirusStopbyIndex:end));
%         DifferenceFilterPos = MaxIntensity/5;
%         DifferenceFilterNeg = -MaxIntensity/5;
        FilteredDiffTracePos = DifferenceTrace > DifferenceFilterPos;
        FilteredDiffTraceNeg = DifferenceTrace < DifferenceFilterNeg;
        RangeToFilterDifference = NumberPointsBack;
%         FilterOption='First Only';
%             [FilteredDiffTracePos] = Filter_Adjacent_Gradient_Results(...
%                 RangeToFilterDifference,FilteredDiffTracePos,FilterOption,[]);
%             [FilteredDiffTraceNeg] = Filter_Adjacent_Gradient_Results(...
%                 RangeToFilterDifference,FilteredDiffTraceNeg,FilterOption,[]);
            [DiffNegClusterData] = Analyze_Gradient_Clusters(FilteredDiffTraceNeg,5);
            [DiffPosClusterData] = Analyze_Gradient_Clusters(FilteredDiffTracePos,5);
    
    % Calculate difference between raw data values and a running median
    % that is of a wider range, this is to identify transient spikes in the
    % data which might indicate a  self quenched fusion event that returns to a similar
    % intensity value
        RunMedHalfLength = 5;
        StartIdx = RunMedHalfLength + 1;
        EndIdx = length(TraceRunMedian.Trace)-RunMedHalfLength;
        TraceRunMedianWiderRange = zeros(length(StartIdx:EndIdx),1);
        %We also set up a vector with the actual frame numbers
        %corresponding to each index position.
        SpikeFrameNumbers = TraceRunMedian.FrameNumbers(StartIdx:EndIdx);

        for n = StartIdx:EndIdx
            TraceRunMedianWiderRange(n-RunMedHalfLength) = median(TraceRunMedian.Trace(n-RunMedHalfLength:n+RunMedHalfLength));
        end
        SpikeTrace = TraceRunMedian.Trace(SpikeFrameNumbers) - TraceRunMedianWiderRange; 
        STDFilterFactorSpike = 5;
        SpikeFilter = STDFilterFactorSpike* std(SpikeTrace(VirusStopbyIndex:end));
        FilteredSpikeTrace = SpikeTrace > SpikeFilter;
        RangeToFilterSpike = 12;
        FilterOption='First Only';
        [FilteredSpikeTrace] = Filter_Adjacent_Gradient_Results(...
            RangeToFilterSpike,FilteredSpikeTrace,FilterOption,[]);
    
        
    %Plot gradient & filters
        set(0,'CurrentFigure',FigureHandles.GradientWindow)
        cla
        plot(TraceRunMedian.FrameNumbers,GradTraceRunMed,'r-');
        hold on
%         plot(TraceRunMedian.FrameNumbers,DifferenceTrace ,'g--');
%         plot(SpikeFrameNumbers,SpikeTrace,'m--');
        plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*NegativeGradFilter,'g--');
        plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*PositiveGradFilter,'g--');
        if strcmp(ProbableDockingEvent, 'y')
            plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*DockGradientFilter,'b--');
        end
        TitleInfo = strcat('Event =',num2str(UniversalData.TraceNumber),'/',num2str(UniversalData.NumTraces),...
            '; STD=',num2str(std(GradTraceRunMed(VirusStopbyIndex:end))));
        title(TitleInfo)
        hold off
        Draw_Events_On_Plot(NegFilteredGradTrace,TraceRunMedian.FrameNumbers,FigureHandles.GradientWindow,'r-')
        Draw_Events_On_Plot(PosFilteredGradTrace,TraceRunMedian.FrameNumbers,FigureHandles.GradientWindow,'k-')
        
        set(0,'CurrentFigure',FigureHandles.DiagnosticWindow)
        cla
        hold on
        plot(TraceRunMedian.FrameNumbers,DifferenceTrace ,'g-');
        plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*DifferenceFilterNeg,'g--');
        plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*DifferenceFilterPos,'g--');
            
        plot(SpikeFrameNumbers,SpikeTrace,'m-');
        plot(SpikeFrameNumbers,ones(1,length(SpikeFrameNumbers))*SpikeFilter,'m--');
        hold off
        
        Draw_Events_On_Plot(FilteredDiffTraceNeg,TraceRunMedian.FrameNumbers,FigureHandles.DiagnosticWindow,'k:')
        Draw_Events_On_Plot(FilteredDiffTracePos,TraceRunMedian.FrameNumbers,FigureHandles.DiagnosticWindow,'b:')
        Draw_Events_On_Plot(FilteredSpikeTrace,SpikeFrameNumbers,FigureHandles.DiagnosticWindow,'r:')
 
        if strcmp(Options.TypeofFusionData,'SLBSelfQuench')
            [DockingData]=Compile_Docking_Data(ProbableDockingEvent,DockFilteredRunMed,VirusStopbyIndex,...
                UniversalData.FrameAllVirusStoppedBy,FrameNumbersDockTrace,RangeToFilterDock);
        else 
            DockingData = [];
        end
        
        
    % Compile all of the trace and gradient data into one structure
        TraceGradData.CroppedRawTrace = CurrTraceCropped;
        TraceGradData.TraceRunMedian = TraceRunMedian;
        TraceGradData.GradTraceRunMed = GradTraceRunMed;
        TraceGradData.PositiveGradFilter = PositiveGradFilter;
        TraceGradData.NegativeGradFilter = NegativeGradFilter;
        TraceGradData.NegFilteredGradTrace = NegFilteredGradTrace;
        TraceGradData.PosFilteredGradTrace= PosFilteredGradTrace;
        TraceGradData.RangeToFilterPositive = RangeToFilterPositive;
        TraceGradData.RangeToFilterNegative = RangeToFilterNegative;

        TraceGradData.SpikeTrace = SpikeTrace;
        TraceGradData.SpikeFilter = SpikeFilter;
        TraceGradData.RangeToFilterSpike = RangeToFilterSpike;
        TraceGradData.FilteredSpikeTrace= FilteredSpikeTrace;
        TraceGradData.SpikeFrameNumbers = SpikeFrameNumbers;

        TraceGradData.DifferenceTrace = DifferenceTrace;
        TraceGradData.DifferenceFilterNeg = DifferenceFilterNeg;
        TraceGradData.DifferenceFilterPos = DifferenceFilterPos;
        TraceGradData.DiffPosClusterData = DiffPosClusterData;
        TraceGradData.DiffNegClusterData = DiffNegClusterData;
        TraceGradData.DiffTraceFrameNumbers = TraceRunMedian.FrameNumbers;
        TraceGradData.FilteredDiffTracePos = FilteredDiffTracePos;
        TraceGradData.FilteredDiffTraceNeg = FilteredDiffTraceNeg;
        TraceGradData.RangeToFilterDifference= RangeToFilterDifference;
end

function [DockingData]=Compile_Docking_Data(ProbableDockingEvent,DockFilteredRunMed,VirusStopbyIndex,VirusStopbyFrameNumber,...
    FrameNumbersDockTrace,RangeToFilterDock)
DockingData.VirusStopbyIndex= VirusStopbyIndex;
DockingData.VirusStopbyFrameNumber= VirusStopbyFrameNumber;
DockingData.ProbableDockingEvent = ProbableDockingEvent; 
DockingData.DockFilteredRunMed =DockFilteredRunMed;
DockingData.FrameNumbersDockTrace=FrameNumbersDockTrace;
DockingData.RangeToFilterDock= RangeToFilterDock;
end

function Draw_Events_On_Plot(FilteredTrace,FrameNumbers,FigureHandle,LineType)
    set(0,'CurrentFigure',FigureHandle)
    hold on
    NumberofEvents = sum(FilteredTrace);
    EventFrameNumbers = FrameNumbers(FilteredTrace);
    LineToPlot = ylim;
    for j= 1:NumberofEvents
        XToPlot = [EventFrameNumbers(j), EventFrameNumbers(j)];
        plot(XToPlot,LineToPlot,LineType)
    end
    hold off
end