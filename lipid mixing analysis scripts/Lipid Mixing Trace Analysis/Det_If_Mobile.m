function [DockingData] = Det_If_Mobile(DockingData,PHdropFrameNum,Options)

    % define default values
    pHtoStopNumFrames = NaN;
    pHtoStopTime = NaN;
    StopFrameNum = NaN;
    IsMobile = 'n';    
    
% The docking data will be empty if the type of fusion data is such that we
% don't search for a docking event (i.e. fusion to a tethered vesicle) or
% if this event was obviously not a docking event as determined in the
% define and apply gradient filter function.
if ~isempty(DockingData)
    if DockingData.ProbableDockingEvent == 'y'
        FrameToStopBy = DockingData.VirusStopbyFrameNumber;
        DockFilteredRunMed = DockingData.DockFilteredRunMed;
        FrameNumbersDockTrace = DockingData.FrameNumbersDockTrace;
        NumDockEvents = sum(DockFilteredRunMed);
        DockEventFrameNumbers = FrameNumbersDockTrace(DockFilteredRunMed);

        if NumDockEvents == 1 && DockEventFrameNumbers(1) < FrameToStopBy
            %clean event    
            IsMobile = 'y';
            StopFrameNum = DockEventFrameNumbers(1);
            pHtoStopNumFrames = DockEventFrameNumbers(1)-PHdropFrameNum; %can be negative!!
            pHtoStopTime = pHtoStopNumFrames*Options.TimeInterval;
        elseif NumDockEvents > 1 && DockEventFrameNumbers(end) < FrameToStopBy
            %pretty clean event, bounced around before stopping, choose
            %last docking event as the real one
            IsMobile = 'y';
            StopFrameNum = DockEventFrameNumbers(end);
            pHtoStopNumFrames = DockEventFrameNumbers(end)-PHdropFrameNum; %can be negative!!
            pHtoStopTime = pHtoStopNumFrames*Options.TimeInterval;
        elseif NumDockEvents == 0
            IsMobile = 'n';
    %             Reason_Failed = 'No Docking Event Detected';
    %             Cross_Out_Plot(TraceWindow,Reason_Failed)
    %                 Stats_Of_Failures.NoDockEvent = Stats_Of_Failures.NoDockEvent + 1;
    %             ShouldISkipToNextSUV = 'y';
        else
            IsMobile = 'n';
        end
    else
    end
end

DockingData.StopFrameNum = StopFrameNum;
DockingData.IsMobile = IsMobile;
DockingData.pHtoStopNumFrames= pHtoStopNumFrames;
DockingData.pHtoStopTime = pHtoStopTime;
    
end