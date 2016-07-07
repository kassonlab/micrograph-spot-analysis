function [Options] = Setup_Options(DefaultPathname)

    Options.Label = 'Test';
        %This is label for the output file
%     Options.DataLabelSuffix = '-1';
        %This suffix will be added to the end of every datafile
        
    Options.FrameAllVirusStoppedBy = NaN;
        % This value will not be used if it was already predefined in the
        % raw trace data. It will also not be used if your viruses are not mobile.
    
    Options.TypeofFusionData  = 'TetheredVesicle';
   
    Options.TimeInterval = .288; %time per frame, in sec
    
    Options.FrameToStartAnalysis = 1; 
%         If you want to skip over some frames at the beginning
    Options.FrametoEndAnalysis = NaN;
        %NaN to indicate length of current trace as frame to end analysis
    
    Options.StartingTraceNumber = 1;
        % If you want to skip to a specific trace number (i.e. for debugging)
    Options.AdditionalFocusFrameNumbers = [];
        % e.g. [361, 369, 370, 678, 679, 681];
        %[] if you don't want to add any additional focus frame numbers
        
    Options.MinImageShow = 90;
    Options.MaxImageShow = 200;
    
end