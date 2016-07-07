function [AnalyzedTraceData,StatsOfDesignations] = Compile_Analyzed_Trace_Data(UniversalData,FusionData,...
    CurrentVirusData,StatsOfDesignations,DockingData,AnalyzedTraceData,ChangedByUser,Options,...
    TraceGradData,SaveOption)

% First we set up the default values, then change them as needed depending
% on the designation
CurrentVirusData.ChangedByUser = ChangedByUser;
CurrentVirusData.Designation = FusionData.Designation;

CurrentVirusData.TraceGradData = TraceGradData;
CurrentVirusData.DockingData = DockingData;
CurrentVirusData.FusionData = FusionData;

% CurrentVirusData.StopFrameNum = DockingData.StopFrameNum;
% CurrentVirusData.IsMobile = DockingData.IsMobile;
% CurrentVirusData.pHtoStopNumFrames= DockingData.pHtoStopNumFrames;
% CurrentVirusData.pHtoStopTime = DockingData.pHtoStopTime;
   
% CurrentVirusData.Fuse1FrameNum = NaN;
% CurrentVirusData.Fuse2FrameNum = NaN;
% CurrentVirusData.pHtoFuse1Time = NaN;
% CurrentVirusData.pHtoFuse2Time = NaN;
% CurrentVirusData.pHtoFuse1NumFrames = NaN;
% CurrentVirusData.pHtoFuse2NumFrames = NaN;
CurrentVirusData.FusionData.pHDropFrameNumber = UniversalData.pHDropFrameNumber;
CurrentVirusData.FusionData.TimeInterval = Options.TimeInterval;

CurrentVirusData.IntensityJumpFuse1 = NaN;
CurrentVirusData.IntensityBeforeFuse1 = NaN;
CurrentVirusData.IntensityAfterFuse1 = NaN;
CurrentVirusData.IntensityJumpFuse2 = NaN;
CurrentVirusData.IntensityBeforeFuse2 = NaN;
CurrentVirusData.IntensityAfterFuse2 = NaN;

if strcmp(FusionData.Designation,'1 Fuse')
    
%     CurrentVirusData.Fuse1FrameNum = FusionData.FuseFrameNumbers(1);
%     CurrentVirusData.pHtoFuse1Time = FusionData.pHtoFusionTime(1);
%     CurrentVirusData.pHtoFuse1NumFrames = FusionData.pHtoFusionNumFrames(1);
%     
        StatsOfDesignations.Fuse1 = StatsOfDesignations.Fuse1 + 1;

elseif strcmp(FusionData.Designation,'2 Fuse')

%Save the relevant data for the fusion event
   
%         CurrentVirusData.Fuse1FrameNum = FusionData.FuseFrameNumbers(1);
%         CurrentVirusData.Fuse2FrameNum = FusionData.FuseFrameNumbers(2);
%         CurrentVirusData.pHtoFuse1Time = FusionData.pHtoFusionTime(1);
%         CurrentVirusData.pHtoFuse2Time = FusionData.pHtoFusionTime(2);
%         CurrentVirusData.pHtoFuse1NumFrames = FusionData.pHtoFusionNumFrames(1);
%         CurrentVirusData.pHtoFuse2NumFrames = FusionData.pHtoFusionNumFrames(2);                 

        StatsOfDesignations.Fuse2 = StatsOfDesignations.Fuse2 + 1;

elseif strcmp(FusionData.Designation,'No Fusion')                    

    StatsOfDesignations.NoFuse = StatsOfDesignations.NoFuse + 1;

else 
    StatsOfDesignations.Other = StatsOfDesignations.Other + 1;
end

% Tally up the total number of events that we have recorded so far
    NumNoFuseRecorded = StatsOfDesignations.NoFuse;
    NumFuse2Recorded = StatsOfDesignations.Fuse2;
    NumFuse1Recorded = StatsOfDesignations.Fuse1;
TotalEventsRecorded = NumFuse1Recorded + NumFuse2Recorded + NumNoFuseRecorded + ...
    StatsOfDesignations.Other;

if TotalEventsRecorded == 1
    AnalyzedTraceData = CurrentVirusData;
else
    AnalyzedTraceData(TotalEventsRecorded) = CurrentVirusData;
end

% OtherDataToSave.Options = Options;
% OtherDataToSave.UniversalData = UniversalData;

if rem(TotalEventsRecorded,5) == 0 && strcmp(SaveOption,'WritetoDisk')
    Save_Data_At_Each_Step(AnalyzedTraceData,Options,OtherDataToSave);
end

end


function Save_Data_At_Each_Step(AnalyzedTraceData,Options,OtherDataToSave)

DataToSave.LastTraceNumberSaved = OtherDataToSave.UniversalData.TraceNumber;
DataToSave.OtherDataToSave = OtherDataToSave;

    if ~isempty(AnalyzedTraceData)
        DataToSave.CombinedAnalyzedTraceData = AnalyzedTraceData;
    end

    save(strcat(Options.IncompleteSavePath,'/',Options.Label,'.mat'),'DataToSave');
end