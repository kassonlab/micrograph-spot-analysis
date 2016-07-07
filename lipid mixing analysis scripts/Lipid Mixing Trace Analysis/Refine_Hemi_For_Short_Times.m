function [ShouldISkipToNextSUV,Stats_Of_Failures,...
    Fuse_FrameNums,...
    DtoF_NumFrames,...
    DtoF_Time] =...
        Refine_Hemi_For_Short_Times(DtoF_Time,Type_Of_Fusion,Fuse_FrameNums,...
            TrToAnalyze_FrameNums,TraceWindow,Stats_Of_Failures,CurrTrace,TimeInterval,Dock_FrameNums)

ShouldISkipToNextSUV = 'n';

FastDtoF_Cutoff = 2; %in sec
if DtoF_Time(1) < FastDtoF_Cutoff && strcmp(Type_Of_Fusion,'Hemi or Hemi-Then-Full')
    %Select range to look at from raw data (not running median)
        HalfRange = 30;
        OldHFuseFrameNum = Fuse_FrameNums(1);
        %Make sure that the selection range doesn't go outside the
        %data
            while ((OldHFuseFrameNum + HalfRange) > max(TrToAnalyze_FrameNums))
                HalfRange = HalfRange - 1;
            end

    if HalfRange < 5
        Reason_Failed = 'Fast fuse event at end of stream' ;
        Cross_Out_Plot(TraceWindow,Reason_Failed)
            Stats_Of_Failures.FastFuseEnd = Stats_Of_Failures.FastFuseEnd + 1;
        ShouldISkipToNextSUV = 'y';
    end
    %Calc int after hemi-fuse
            %HF_Int = median(CurrTrace(OldHFuseFrameNum+1:OldHFuseFrameNum+HalfRange));

        %Collapsed hemifusion test.  The first frame closer to the HF
        %intensity than to the docking intensity (following
        %docking) is considered the hemifusion event.
            %FramesToSelect_ForHF = Dock_FrameNums(1):OldHFuseFrameNum+HalfRange;
            %TraceSelected_HF = CurrTrace(FramesToSelect_ForHF);
            %    DistToDocked = (TraceSelected_HF - Docked_Int).^2;
            %    DistToHF = (TraceSelected_HF - HF_Int).^2;
            %Collapsed_Trace_HF = DistToHF < DistToDocked;
            %Collapsed_Trace_HF_FrameNums = FramesToSelect_ForHF(Collapsed_Trace_HF);
                %If the newly calculated HF frame num is > the
                %old one, then use the old one.  We assume that the
                %earliest calculated value is probably correct.
            %NewHFFrameNum = min(Collapsed_Trace_HF_FrameNums(1),OldHFuseFrameNum);

        %Hemifusion test--gradient on the raw data.
            FramesToSelect_ForHF = Dock_FrameNums(1):OldHFuseFrameNum+HalfRange;
            TraceSelected_HF = CurrTrace(FramesToSelect_ForHF);
            Grad_HF = gradient(TraceSelected_HF);
                Neg_Grad_Filter = Grad_HF < (min(Grad_HF) - .2*min(Grad_HF));
                HF_Filtered_FrameNums = FramesToSelect_ForHF(Neg_Grad_Filter);
            NewHFFrameNum = min(HF_Filtered_FrameNums(1), OldHFuseFrameNum);

            %Make sure we don't change the fusion frame too much by
            %accident.  Assume the old value is more accurate if
            %the new estimate is off by more than 4 from the old.
            if abs(NewHFFrameNum - OldHFuseFrameNum) > 4
                NewHFFrameNum = OldHFuseFrameNum;
            end

    %Overwrite old DtoF values with the newly calculated ones
        Fuse_FrameNums(1) = NewHFFrameNum;
        DtoF_NumFrames(1) = Fuse_FrameNums(1) - Dock_FrameNums(1);
        DtoF_Time(1) = TimeInterval*(DtoF_NumFrames(1));
end

%Bob, Oct 2011: By analyzing ~40 points by eye, I noticed that the hemifusion
%frame (especially at fast timescales) was consistently 1 frame too
%slow.  So I subtract one frame from the HF frame num.  This
%should reduce the total amount of error in the data quantification.
    Fuse_FrameNums(1) = Fuse_FrameNums(1) - 1;
    DtoF_NumFrames(1) = Fuse_FrameNums(1) - Dock_FrameNums(1);
    DtoF_Time(1) = TimeInterval*(DtoF_NumFrames(1));

    %This if statement is also added in to account for times when
    %the subtraction happened to make the DtoF time impractically
    %small.
    if DtoF_Time(1) <= 0
        Reason_Failed = 'DtoF Time <= 0' ;
        Cross_Out_Plot(TraceWindow,Reason_Failed)
            Stats_Of_Failures.DtoFLess0 = Stats_Of_Failures.DtoFLess0 + 1;
        ShouldISkipToNextSUV = 'y';            
    end