function [Neg_Filtered_RunMed,Pos_Filtered_RunMed, ExitFlag2] = Re_Filter_Gradient_Results(Neg_Filtered_RunMed,...
    Pos_Filtered_RunMed,Grad_Trace_RunMed, Max_Grad_Near_Dock,TypeofFusionData)

%By default, no exit flag
ExitFlag2 = 0;

%Filter out gradients that occur next to each other, if needed
if sum(Neg_Filtered_RunMed) > 1
    IdxToCheck = find(Neg_Filtered_RunMed>0);

    if strcmp(TypeofFusionData, 'Normal')
        RangeToFilter = 10;
    elseif strcmp(TypeofFusionData, 'StrobeLight,NormalDye')
        RangeToFilter =  5;
    elseif strcmp(TypeofFusionData, 'ContinuousLightVesicle')
        RangeToFilter =  10;
    elseif strcmp(TypeofFusionData, 'StrobeLight,HiDye')
       RangeToFilter =  5;
    elseif strcmp(TypeofFusionData, 'SLBSelfQuench')
       RangeToFilter =  5;
    end

    OldIdx = -1*RangeToFilter; %Guarantees that we always see the first one.
    for b = 1:length(IdxToCheck)
        NewIdx = IdxToCheck(b);

        if NewIdx - OldIdx <= RangeToFilter
            IdxToFilterOut = NewIdx;
        else
            IdxToFilterOut = IdxToCheck(b)+1:IdxToCheck(b)+1+RangeToFilter;
            OldIdx = NewIdx;
        end

        Neg_Filtered_RunMed(IdxToFilterOut) = 0;
    end
end

if strcmp(TypeofFusionData, 'Normal') || strcmp(TypeofFusionData, 'StrobeLight,NormalDye')
    if sum(Pos_Filtered_RunMed) > 1
        %First select out the maximum positive gradient and filter
        %out any points around it that also passed the filter test
        %This means that the maximum gradient will be recognized as
        %the docking event, rather than the first frame above the
        %filter test.
        IdxMaxGrad = find(Grad_Trace_RunMed == Max_Grad_Near_Dock);
            IdxMaxGrad = IdxMaxGrad(1); %This is to help if there were two frames with the exact same max grad
        HalfLengthToFilter = 15;
        if (IdxMaxGrad-HalfLengthToFilter) > 0 
            IdxBelowMax = IdxMaxGrad-HalfLengthToFilter:IdxMaxGrad-1;
            Pos_Filtered_RunMed(IdxBelowMax) = 0;
        elseif (IdxMaxGrad-HalfLengthToFilter) < 0 && IdxMaxGrad > 1
            IdxBelowMax = 1:IdxMaxGrad-1;
            Pos_Filtered_RunMed(IdxBelowMax) = 0;
        else
            %This means the maximum gradient is the first index (too close to
            %beginning)
%             ExitFlag2 = 1;
        end


        if (IdxMaxGrad + HalfLengthToFilter) <= length(Grad_Trace_RunMed)
            IdxAboveMax = IdxMaxGrad + 1: IdxMaxGrad+HalfLengthToFilter;
            Pos_Filtered_RunMed(IdxAboveMax) = 0;
        else
            %This means the maximum gradient is within the last few frames (too
            %close to end)
%             ExitFlag2 = 1;
        end
    end
   
elseif strcmp(TypeofFusionData, 'StrobeLight,HiDye')|| strcmp(TypeofFusionData, 'ContinuousLightVesicle')...
        || strcmp(TypeofFusionData, 'SLBSelfQuench')
    if strcmp(TypeofFusionData, 'ContinuousLightVesicle')
       RangeToFilter =  10;
    elseif strcmp(TypeofFusionData, 'StrobeLight,HiDye')
       RangeToFilter =  3;
    elseif strcmp(TypeofFusionData, 'SLBSelfQuench')
       RangeToFilter =  5;
    end
    
    %Filter out gradients that occur next to each other, if needed
    if sum(Pos_Filtered_RunMed) > 1
        IdxToCheck = find(Pos_Filtered_RunMed>0);
        OldIdx = -1*RangeToFilter; %Guarantees that we always see the first one.
        for b = 1:length(IdxToCheck)
            NewIdx = IdxToCheck(b);

            if NewIdx - OldIdx <= RangeToFilter
                IdxToFilterOut = NewIdx;
            else
                IdxToFilterOut = IdxToCheck(b)+1:IdxToCheck(b)+1+RangeToFilter;
                OldIdx = NewIdx;
            end

            Pos_Filtered_RunMed(IdxToFilterOut) = 0;
        end
    end


end