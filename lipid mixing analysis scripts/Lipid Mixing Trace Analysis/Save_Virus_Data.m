function [Fuse1Data,Fuse2Data,NoFuseData,StatsOfDesignations] = Save_Virus_Data(i,Designation,CurrentVirusData,StatsOfDesignations,...
    Dock_FrameNums,IsMobile,pHtoStopNumFrames,pHtoStopTime,StopFrameNum,Fuse_FrameNums,pHtoF_Time,...
    pHtoF_NumFrames,PHdropFrameNum,TimeInterval,Fuse1Data,Fuse2Data,NoFuseData,ChangedByUser,IncompleteSavePath,Label)

CurrentVirusData.ChangedByUser = ChangedByUser;

if strcmp(Designation,'1 Fuse')
                    
            %Save the relevant data for the fusion event
                    CurrentVirusData.Dock_FrameNum = Dock_FrameNums(1);
                    CurrentVirusData.IsMobile = IsMobile;
                    CurrentVirusData.pHtoStopNumFrames= pHtoStopNumFrames;
                    CurrentVirusData.pHtoStopTime = pHtoStopTime;
                    CurrentVirusData.StopFrameNum = StopFrameNum;
                    CurrentVirusData.Fuse1FrameNum = Fuse_FrameNums(1);
                    CurrentVirusData.pHtoFuse1Time = pHtoF_Time(1);
                    CurrentVirusData.pHtoFuse1NumFrames = pHtoF_NumFrames(1);
                    CurrentVirusData.PHdropFrameNum = PHdropFrameNum;
                    CurrentVirusData.TimeInterval = TimeInterval;
                    
                    CurrentVirusData.IntensityJumpFuse1 = NaN;
                    CurrentVirusData.IntensityBeforeFuse1 = NaN;
                    CurrentVirusData.IntensityAfterFuse1 = NaN;
                    
                    StatsOfDesignations.Fuse1 = StatsOfDesignations.Fuse1 + 1;
                    NumFuse1Recorded = StatsOfDesignations.Fuse1;

                if NumFuse1Recorded == 1
                    Fuse1Data = CurrentVirusData;
                else
                    Fuse1Data(NumFuse1Recorded) = CurrentVirusData;
                end                    
elseif strcmp(Designation,'2 Fuse')

            %Save the relevant data for the fusion event
                    CurrentVirusData.Dock_FrameNum = Dock_FrameNums(1);
                    CurrentVirusData.IsMobile = IsMobile;
                    CurrentVirusData.pHtoStopNumFrames= pHtoStopNumFrames;
                    CurrentVirusData.pHtoStopTime = pHtoStopTime;
                    CurrentVirusData.StopFrameNum = StopFrameNum;
                    CurrentVirusData.Fuse1FrameNum = Fuse_FrameNums(1);
                    CurrentVirusData.Fuse2FrameNum = Fuse_FrameNums(2);
                    CurrentVirusData.pHtoFuse1Time = pHtoF_Time(1);
                    CurrentVirusData.pHtoFuse2Time = pHtoF_Time(2);
                    CurrentVirusData.pHtoFuse1NumFrames = pHtoF_NumFrames(1);
                    CurrentVirusData.pHtoFuse2NumFrames = pHtoF_NumFrames(2);
                    CurrentVirusData.PHdropFrameNum = PHdropFrameNum;
                    CurrentVirusData.TimeInterval = TimeInterval;
                    
                    CurrentVirusData.IntensityJumpFuse1 = NaN;
                    CurrentVirusData.IntensityBeforeFuse1 = NaN;
                    CurrentVirusData.IntensityAfterFuse1 = NaN;
                    CurrentVirusData.IntensityJumpFuse2 = NaN;
                    CurrentVirusData.IntensityBeforeFuse2 = NaN;
                    CurrentVirusData.IntensityAfterFuse2 = NaN;

                    
                    StatsOfDesignations.Fuse2 = StatsOfDesignations.Fuse2 + 1;
                    NumFuse2Recorded = StatsOfDesignations.Fuse2;

                if NumFuse2Recorded == 1
                    Fuse2Data = CurrentVirusData;
                else
                    Fuse2Data(NumFuse2Recorded) = CurrentVirusData;
                end                    
elseif strcmp(Designation,'No Fusion')                    
     %Save the relevant data for the no fusion event
                    CurrentVirusData.Dock_FrameNum = Dock_FrameNums(1);
                    CurrentVirusData.IsMobile = IsMobile;
                    CurrentVirusData.pHtoStopNumFrames= pHtoStopNumFrames;
                    CurrentVirusData.pHtoStopTime = pHtoStopTime;
                    CurrentVirusData.StopFrameNum = StopFrameNum;
                    CurrentVirusData.PHdropFrameNum = PHdropFrameNum;
                    CurrentVirusData.TimeInterval = TimeInterval;

                    StatsOfDesignations.NoFuse = StatsOfDesignations.NoFuse + 1;
                    NumNoFuseRecorded = StatsOfDesignations.NoFuse;

                if NumNoFuseRecorded == 1
                    NoFuseData = CurrentVirusData;
                else
                    NoFuseData(NumNoFuseRecorded) = CurrentVirusData;
                end
                
end

Save_Data_At_Each_Step(Fuse1Data,Fuse2Data,NoFuseData,IncompleteSavePath,Label,i);

end


function Save_Data_At_Each_Step(Fuse1_Data,Fuse2_Data,NoFus_Data,IncompleteSavePath,Label,i)

Useful_Data_To_Save.LastTraceNumberSaved = i;

    if ~isempty(Fuse1_Data)
        Useful_Data_To_Save.Combined_Fuse1_Data = Fuse1_Data;
    end

    if ~isempty(NoFus_Data)
        Useful_Data_To_Save.Combined_NoFuse_Data = NoFus_Data;
    end

    if ~isempty(Fuse2_Data)
        Useful_Data_To_Save.Combined_Fuse2_Data = Fuse2_Data;
    end

    save(strcat(IncompleteSavePath,'/',Label,'.mat'),'Useful_Data_To_Save');
end