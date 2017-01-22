function [ClusterData] = Analyze_Trace_Clusters(FilteredTrace,ClusterRange)
    % This function takes a filtered trace (logical vector) and finds 
    % clusters of 1's that are within a certain distance from each other, 
    % specified by the ClusterRange. It then returns a structure ClusterData 
    % which contains various metrics about the indices and size etc. of 
    % each cluster.

    % Pre-define some values
    NumberOfClusters = 0;
    ClusterStartIndices = [];
    ClusterSizes = [];

    % Cluster analysis
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