function [RerunThisRound, CorrectedAnalysisData, ErrorCounter] = Correct_Designations(IncorrectPlotIndices,...
    PreviousAnalysisData,CurrentTraceRange,CorrectedAnalysisData, ErrorCounter,Options)

    RerunThisRound = 'n';
    NumberIncorrect = length(IncorrectPlotIndices);

    for j = 1:NumberIncorrect
        CurrentIndex = floor(IncorrectPlotIndices(j));

        TraceNumberIndex = CurrentTraceRange(CurrentIndex);

        if  strcmp(Options.QuickModeNoCorrection,'y')
            CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Not Changed';
            
        elseif strcmp(Options.QuickModeNoCorrection,'n')
            
            CurrentDesignationCode = rem(IncorrectPlotIndices(j), 1);
            if CurrentDesignationCode == 0.0
                if strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, 'No Fusion')
                    % There has been a mistake, re-run the last round
                    RerunThisRound = 'y';
                    break
                else

                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = 'No Fusion';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = 'No Fusion';
                end
            elseif CurrentDesignationCode > 0.09 && CurrentDesignationCode < 0.11
                if strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, '1 Fuse')
                    % There has been a mistake, re-run the last round
                    RerunThisRound = 'y';
                    break
                else
                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = '1 Fuse';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = '1 Fuse';
                end
            elseif CurrentDesignationCode > 0.11 && CurrentDesignationCode < 0.13
                if ~strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, '1 Fuse')
                    % There has been a mistake, re-run the last round
                    RerunThisRound = 'y';
                    break
                else
                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Correct Designation, Incorrect Wait Time';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = '1 Fuse';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = '1 Fuse';
                end
            elseif CurrentDesignationCode > 0.19 && CurrentDesignationCode < 0.21
                if strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, '2 Fuse') ...
                        || strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, 'Slow')
                    % There has been a mistake, re-run the last round
                    RerunThisRound = 'y';
                    break
                else
                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = 'Slow';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = 'Slow';
                end
            elseif CurrentDesignationCode > 0.89 && CurrentDesignationCode < 0.91
                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = 'Strange-Ignore';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = 'Strange-Ignore';
            else
                % There has been an error, re-run the last round to avoid crash
                RerunThisRound = 'y';
                break
            end
        end
    end
    
    if  strcmp( RerunThisRound,'n')
        ErrorCounter = ErrorCounter + NumberIncorrect;
    end
end