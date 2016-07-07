function [FilteredTrace] = Filter_Adjacent_Gradient_Results(RangeToFilter,TraceToFilter,FilterOption,Grad_Trace_RunMed)

   FilteredTrace=TraceToFilter; 

if strcmp(FilterOption,'First Only') 
    % TheFilterIsAppliedattheEnd
elseif strcmp(FilterOption,'Max Then First') 
    if sum(FilteredTrace) > 1
        %First select out the maximum positive gradient and filter
        %out any points around it that also passed the filter test
        %This means that the maximum gradient will be recognized as
        %the docking event, rather than the first frame above the
        %filter test.
        IndexMaxGrad = find(Grad_Trace_RunMed == max(Grad_Trace_RunMed));
            IndexMaxGrad = IndexMaxGrad(1); %This is to help if there were two frames with the exact same max grad
        if IndexMaxGrad == 1
            % ThenWeOnlyDeleteGoingForward 
            if (IndexMaxGrad + RangeToFilter) <= length(Grad_Trace_RunMed)
                IndexAboveMax = IndexMaxGrad + 1: IndexMaxGrad+RangeToFilter;
            else
                IndexAboveMax = IndexMaxGrad + 1: length(Grad_Trace_RunMed); 
            end
            FilteredTrace(IndexAboveMax) = 0;
            
        elseif IndexMaxGrad == length(Grad_Trace_RunMed)  
            % ThenWeOnlyDeleteGoingBackward
            if (IndexMaxGrad-RangeToFilter) > 0 
                IndexBelowMax = IndexMaxGrad-RangeToFilter:IndexMaxGrad-1;
            elseif (IndexMaxGrad-RangeToFilter) <= 0
                IndexBelowMax = 1:IndexMaxGrad-1;   
            end
            FilteredTrace(IndexBelowMax) = 0; 
            
        else
            % ThenWeDeleteonEitherSideoftheMax
            if (IndexMaxGrad-RangeToFilter) > 0 
                IndexBelowMax = IndexMaxGrad-RangeToFilter:IndexMaxGrad-1;
            elseif (IndexMaxGrad-RangeToFilter) <= 0
                IndexBelowMax = 1:IndexMaxGrad-1;   
            end
            FilteredTrace(IndexBelowMax) = 0; 
            
            if (IndexMaxGrad + RangeToFilter) <= length(Grad_Trace_RunMed)
                IndexAboveMax = IndexMaxGrad + 1: IndexMaxGrad+RangeToFilter;
            else
                IndexAboveMax = IndexMaxGrad + 1: length(Grad_Trace_RunMed); 
            end
            FilteredTrace(IndexAboveMax) = 0;
            
        end
    end
end

%Filter out gradients that occur next to each other, if needed
if sum(FilteredTrace) > 1
    IdxToCheck = find(FilteredTrace>0);

    OldIdx = -1*RangeToFilter; %Guarantees that we always see the first one.
    for b = 1:length(IdxToCheck)
        NewIdx = IdxToCheck(b);

        if NewIdx - OldIdx <= RangeToFilter
            IdxToFilterOut = NewIdx;
        else
            IdxToFilterOut = IdxToCheck(b)+1:IdxToCheck(b)+RangeToFilter-2;
            OldIdx = NewIdx;
        end

        FilteredTrace(IdxToFilterOut) = 0;
    end
end

end