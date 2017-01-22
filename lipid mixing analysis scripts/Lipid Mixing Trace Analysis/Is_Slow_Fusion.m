function [TraceGradData,FusionData] = Is_Slow_Fusion(TraceGradData,FusionData,...
        DockingData,DetectionOption,Options)

% We use the difference trace test to identify slow fusion events, which 
% appear to occur over many frames, rather than as a sharp event. Because 
% the difference trace will also identify fast events, we have to do a 
% cross comparison to make sure that we don't identify a fast event as 
% a slow event. There are two different methods of analyzing the difference 
% trace test data - Usual and Cluster Analysis. Right now the Cluster 
% Analysis is more reliable (Bob, December 2016).

% Pre-define some variables
FilteredDiffTraceNeg = TraceGradData.FilteredDiffTraceNeg;
FilteredDiffTracePos = TraceGradData.FilteredDiffTracePos;
RangeToFilterDifference = Options.NumberFramesBackToSubtract;
TraceRunMedian = TraceGradData.TraceRunMedian;

if strcmp(DetectionOption,'Usual Trace Analysis')
    NumSlowFusionEventsPosDetected = sum(FilteredDiffTracePos);
    NumSlowFusionEventsNegDetected = sum(FilteredDiffTraceNeg);
    SlowFusePosFrameNumbers = TraceRunMedian.FrameNumbers(FilteredDiffTracePos);
    SlowFuseNegFrameNumbers = TraceRunMedian.FrameNumbers(FilteredDiffTraceNeg);

    if NumSlowFusionEventsNegDetected == 0 && NumSlowFusionEventsPosDetected == 0
    else
        if NumSlowFusionEventsNegDetected > 0
            if isempty(FusionData.FuseFrameNumbers)
                FusionData.Designation = 'Slow';
            else
                for j= 1:NumSlowFusionEventsNegDetected 
                    % If the potential slow fusion event is nearby a fusion
                    % event that was already detected, assume that everything
                    % is fine and move on. Otherwise, classify this trace as a
                    % slow fusion event.
                    DistanceToFastFusionEvents = abs(FusionData.FuseFrameNumbers - SlowFuseNegFrameNumbers(j));
                    EventsCloseBy = DistanceToFastFusionEvents <= RangeToFilterDifference;
                    if sum(EventsCloseBy) == 0
                        FusionData.Designation = 'Slow';
                        break
                    end
                end
            end

        end

        if NumSlowFusionEventsPosDetected > 0
            if isempty(FusionData.FuseFrameNumbers) && isempty(DockingData.StopFrameNum)
                FusionData.Designation = 'Slow';
            else
                for j= 1:NumSlowFusionEventsPosDetected
                    % If the potential slow fusion event is nearby a fusion or
                    % docking event that was already detected, assume that everything
                    % is fine and move on. Otherwise, classify this trace as a
                    % slow fusion event.
                    OtherEventFrameNumbers = [DockingData.StopFrameNum FusionData.FuseFrameNumbers];
                    DistanceToOtherEvents = abs(OtherEventFrameNumbers - SlowFusePosFrameNumbers(j));
                    EventsCloseBy = DistanceToOtherEvents <= RangeToFilterDifference;
                    if sum(EventsCloseBy) == 0
                        FusionData.Designation = 'Slow';
                        break
                    end
                end
            end

        end
    end

elseif strcmp(DetectionOption,'Cluster Analysis')
    
    ClusterRange = Options.NumFramesBetweenDifferentClusters;
    [TraceGradData.DiffPosClusterData] = Analyze_Trace_Clusters(FilteredDiffTracePos,ClusterRange);
    [TraceGradData.DiffNegClusterData] = Analyze_Trace_Clusters(FilteredDiffTraceNeg,ClusterRange);
    
    DiffNegClusterData = TraceGradData.DiffNegClusterData;
    DiffPosClusterData = TraceGradData.DiffPosClusterData;
    
    NumSlowFusionEventsPosDetected = DiffPosClusterData.NumberOfClusters;
    NumSlowFusionEventsNegDetected = DiffNegClusterData.NumberOfClusters;
    SlowFusePosFrameNumbers = TraceRunMedian.FrameNumbers(DiffPosClusterData.ClusterStartIndices);
    SlowFuseNegFrameNumbers = TraceRunMedian.FrameNumbers(DiffNegClusterData.ClusterStartIndices);
    ClusterSizesPos = DiffPosClusterData.ClusterSizes;
    ClusterSizesNeg = DiffNegClusterData.ClusterSizes;
    
    NumberFramesFastFusionCuttoffPos = Options.NumberFramesBackToSubtract + Options.ClusterSizePosConsideredFastFusion;
    NumberFramesFastFusionCuttoffNeg = Options.NumberFramesBackToSubtract + Options.ClusterSizeNegConsideredFastFusion;

    if NumSlowFusionEventsNegDetected == 0 && NumSlowFusionEventsPosDetected == 0
    else
        if NumSlowFusionEventsNegDetected > 0
            if isempty(FusionData.FuseFrameNumbers)
                if NumSlowFusionEventsNegDetected == 1 && ClusterSizesNeg(1) < ...
                        NumberFramesFastFusionCuttoffNeg && NumSlowFusionEventsPosDetected == 0
                    
                    if strcmp(Options.TypeofFusionData,'TetheredVesicle')
                        FusionData.Designation = 'Slow';
                    else
                        FusionData.Designation = '1 Fuse';
                        FusionData.FuseFrameNumbers = SlowFuseNegFrameNumbers(1);
                    end
                else 
                    FusionData.Designation = 'Slow';
                end
            else
                for j= 1:NumSlowFusionEventsNegDetected 
                    % If the potential slow fusion event is nearby a fusion
                    % event already detected, and the cluster size is small
                    % enough, then assume it is fine. Otherwise, classify
                    % as a slow fusion event.
                    DistanceToFastFusionEvents = abs(FusionData.FuseFrameNumbers - SlowFuseNegFrameNumbers(j));
                    EventsCloseBy = DistanceToFastFusionEvents <= RangeToFilterDifference;
                    NumberBigClusters = sum(ClusterSizesNeg >= NumberFramesFastFusionCuttoffNeg);
                    if sum(EventsCloseBy) == 0 || NumberBigClusters > 0
                        FusionData.Designation = 'Slow';
                        break
                    end
                end
            end

        end

        if NumSlowFusionEventsPosDetected > 0
            if isempty(FusionData.FuseFrameNumbers) && isnan(DockingData.StopFrameNum)
                if NumSlowFusionEventsPosDetected == 1 && ClusterSizesPos(1) < NumberFramesFastFusionCuttoffPos &&...
                        NumSlowFusionEventsNegDetected == 0
                    
                    FusionData.Designation = '1 Fuse';
                    FusionData.FuseFrameNumbers = SlowFusePosFrameNumbers(1);
                else 
                    FusionData.Designation = 'Slow';
                end
            else
                for j= 1:NumSlowFusionEventsPosDetected
                    % If the potential slow fusion event is nearby a fusion or
                    % docking event that was already detected, (And the cluster sizes are small enough)
                    % assume that everything
                    % is fine and move on. Otherwise, classify this trace as a
                    % slow fusion event.
                    if isnan(DockingData.StopFrameNum)
                        OtherEventFrameNumbers = FusionData.FuseFrameNumbers;
                    else 
                        OtherEventFrameNumbers = [DockingData.StopFrameNum FusionData.FuseFrameNumbers];
                    end
                    
                    DistanceToOtherEvents = abs(OtherEventFrameNumbers - SlowFusePosFrameNumbers(j));
                    EventsCloseBy = DistanceToOtherEvents <= RangeToFilterDifference;
                    NumberBigClusters = sum(ClusterSizesPos >= NumberFramesFastFusionCuttoffPos);
                    if sum(EventsCloseBy) == 0 || NumberBigClusters > 0
                        FusionData.Designation = 'Slow';
                        break
                    end
                end
            end

        end
    end
end

end