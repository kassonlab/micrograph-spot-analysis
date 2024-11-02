function [AllResults] = Calculate_Randomness_Parameter(FileNumber,AllResults)

SortedpHtoFList = AllResults(FileNumber).CDFData.SortedpHtoFList;

MeanWaitingTimeSquared = (mean(SortedpHtoFList))^2;
MeanOfSquared = mean(SortedpHtoFList.^2);

RandomnessParameter = (MeanOfSquared - MeanWaitingTimeSquared)/MeanWaitingTimeSquared;
Nmin = 1/RandomnessParameter;

AllResults(FileNumber).RandomnessParameter = RandomnessParameter;
AllResults(FileNumber).Nmin = Nmin;

end