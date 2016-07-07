function [ClusterData] = Analyze_Gradient_Clusters(FilteredTrace,ClusterRange)
NumberOfClusters = 0;
ClusterStartIndices = [];
ClusterSizes = [];

if sum(FilteredTrace) > 1
    IdxToCheck = find(FilteredTrace>0);
    NumberOfClusters = 1;
    ClusterStartIndices = IdxToCheck(1);
    ClusterSizes = 1;
    
    if length(IdxToCheck)<2
        % Move On
    else
        OldIndex = IdxToCheck(1);
        
        for b = 2:length(IdxToCheck)
            NewIndex = IdxToCheck(b);

            if NewIndex - OldIndex <= ClusterRange
                % Same cluster
                ClusterSizes(NumberOfClusters) = ClusterSizes(NumberOfClusters) + (NewIndex - OldIndex);
            else
                % New cluster
                NumberOfClusters = NumberOfClusters + 1;
                ClusterStartIndices(NumberOfClusters) = NewIndex;
                ClusterSizes(NumberOfClusters) = 1;
            end
            OldIndex = NewIndex;
        end
    end
    
end
ClusterData.NumberOfClusters = NumberOfClusters;
ClusterData.ClusterStartIndices = ClusterStartIndices;
ClusterData.ClusterSizes = ClusterSizes;
ClusterData.ClusterRange = ClusterRange;
end