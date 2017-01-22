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
        % e.g. [842:847,901];
        %[] if you don't want to add any additional focus frame numbers
    Options.ChangepHNumber = [];
        % This will change the pH drop number (normally specified during the 
        % Extract Traces From Video program, but added here in case that 
        % number is incorrect and/or you want to play with it).
        % [] If you don't want to change.
        
    Options.MinImageShow = 90;
    Options.MaxImageShow = 200;
   
% --------Analysis Parameters--------
% These parameters determine how the traces will be analyzed, setting 
% thresholds to decide what is considered a fusion event. They should 
% be optimized for your particular application; the values listed below 
% are only suggestions that work well with my data. As different parameters 
% may be needed for fusion to tethered vesicles or to supported bilayers,
% they are listed separately as needed.

    Options.NoFusionAllowedBeforepHDrop = 'y';
        % This will not allow fusion events to be found by the gradient test 
        % before the pH drop (most likely a good idea, as these are probably 
        % spurious events anyway).
    Options.ClipFramesAtEnd = 'y';
        % This will not allow fusion events to be found by the gradient test within a certain number 
        % of frames at the end of the video. Empirically, I have noticed that 
        % very occasionally events are erroneously found within the last couple of 
        % frames, so this is a way to correct that. 'y' OR 'n'
    Options.NumFramesToClip = 2;
        % This is the number of frames from the end of the video in which fusion 
        % events will not be found by the gradient test.
        
        
    if strcmp(Options.TypeofFusionData, 'TetheredVesicle') 
        
        % --------Gradient Test Parameters--------
        Options.RunMedHalfLength = 0; 
            % This is the number of frames on either side which will be used to 
            % calculate the running median trace, which will then be used to calculate 
            % the gradient and other values.
        Options.GradientThreshold = 5;
            % This is the number of standard deviations above which the 
            % gradient must reach (either positive or negative) in order 
            % to be considered a possible fusion event. In the case of tethered 
            % vesicles, the negative value is not used.
        Options.NumFramesBetweenGradientEvents = 10;
            % This is the number of frames apart below which possible fusion events 
            % identified by the gradient test will be considered the same event. 
        
        % --------Spike Trace Test Parameters--------
        Options.RunMedHalfLengthSpike = 5;
            % This is the number of frames on either side which will be used to 
            % calculate a separate running median trace, which will then be used to 
            % calculate the spike trace. This value should be large relative to 
            % the duration of the spikes which you are trying to detect.
        Options.SpikeThreshold = 5;
            % This is the number of standard deviations above which the spike trace 
            % must reach in order to be considered a true spike.
        Options.NumFramesBetweenSpikeEvents = 12;
            % This is the number of frames apart below which spike events will be 
            % considered the same event.
            
        % --------Difference Trace Test Parameters--------
        Options.NumberFramesBackToSubtract = 25;
            % This is the number of frames to subtract backward in order 
            % to create the difference trace.
        Options.DifferenceTraceThreshold = 14;
            % This is the number of standard deviations above which the 
            % difference trace must reach in order to be considered a 
            % possible fusion event.
        Options.ClusterSizePosConsideredFastFusion = 15;
        Options.ClusterSizeNegConsideredFastFusion = 15;
            % This is the number of adjacent frames (clusters) flagged as a fusion event by 
            % the difference trace test below which we will consider them to be a 
            % fast (i.e. normal) fusion event. Separate values are defined for 
            % difference trace values that are positive (i.e. an increase in 
            % fluorescence) or negative (i.e. a decrease in fluorescence). 
            % Note that these values are added to the NumberFramesBackToSubtract 
            % value, since a fluorescence intensity sharp jump will likely have 
            % at least that many adjacent frames which will pass the difference 
            % trace test.
        Options.NumFramesBetweenDifferentClusters = 10;
            % This is the number of frames apart below which possible fusion events 
            % identified by the difference trace test (clusters) will be considered the same event.
        
    elseif strcmp(Options.TypeofFusionData, 'SLBSelfQuench') 
        % WARNING: This option is not fully set up
        % Note that if this option is chosen, there are additional parameters 
        % that are defined in order to identify the frame in which the virus 
        % docks (i.e. stops moving). These are listed in the 
        % Define_and_Apply_Gradient_Filters.m script. If you are doing this 
        % regularly, you may wish to move them here.
        
        Options.RunMedHalfLength = 0; 
            % Defined above
        Options.GradientThreshold = 4.5;
            % Defined above
        Options.NumFramesBetweenFusionEvents_Gradient = 5;
            % Defined above
        Options.RunMedHalfLengthSpike = 5; 
            % Defined above
        Options.SpikeThreshold = 5;
            % Defined above
        Options.NumFramesBetweenSpikeEvents = 12;
            % Defined above
        
        Options.NumberFramesBackToSubtract = 25;
            % Defined above
        Options.DifferenceTraceThreshold = 9;
            % Defined above        
        Options.ClusterSizePosConsideredFastFusion = 15;
        Options.ClusterSizeNegConsideredFastFusion = 5;
        Options.NumFramesBetweenDifferentClusters = 10;
    else
        disp(' Type of fusion data not specified correctly');
    end        
    
end