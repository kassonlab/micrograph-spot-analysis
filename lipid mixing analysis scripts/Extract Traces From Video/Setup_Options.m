function [Options] = Setup_Options()

    Options.Threshold = 80;
        % This is the number of counts above background which will be used
        % to detect virus particles. You will need to optimize this number 
        % for each set of imaging conditions and/or each data set. An optimal 
        % threshold value has been reached when you are able to see that the
        % program is accurately finding the particles in each image. To avoid
        % bias introduced by the optimization process, you should make sure 
        % that changing the optimal threshold value somewhat doesn't 
        % significantly affect your data. Assuming similar particle 
        % densities/intensities between data sets, and that the same imaging 
        % conditions are used, you shouldn't need to change the threshold 
        % value much if at all between data sets.

% ---------Inputs specific to each data set---------
    Options.PHdropFrameNum = 30;
        % This is the frame number in which the pH drop occurs to trigger fusion, 
        % defined as the frame number in which the intensity of the pH indicator 
        % drops below one half its value. Will be used to determine the pH 
        % drop to lipid mixing waiting time in the trace analysis program.
        % NOTE: This number will be ignored if you choose to extract inputs from filename.
    Options.FrameNumToFindParticles = 25;
        % This is the frame number that will be used as the finding image, 
        % to find viral particles. 
        % NOTE: This number will be ignored if you choose to extract inputs from filename.
    Options.FindFramesToAverage = 2;
        % Number of frames before/after to average for the finding image       
    Options.FocusFrameNumbers = NaN;
        % These are frame numbers where a focus event occurred and needs to be 
        % corrected in the intensity traces for each viral particle. This 
        % correction will be performed in the Trace Analysis program. 
        % 'NaN' for no events. [FrameNumber1,FrameNumber2, etc.] for multiple 
        % frames. 
        % NOTE: This number will be ignored if you choose to extract inputs from filename.

    Options.ExtractInputsFromFilename = 'n';
        % 'y' OR 'n'
        % Choose 'y' to automatically extract the following inputs from 
        % the filename: PHdropFrameNum, FrameNumToFindParticles,
        % FocusFrameNumbers. In that case, the numbers put in for those 
        % values above will be ignored. Also, you should examine the notes 
        % in Extract_Analysis_Inputs.m for the proper formatting of your filename.

% ---------Parameters Used To Find Particles/Assess Particle 'Goodness'---------
    Options.MinParticleSize = 4;
        % This is the minimum particle size (defined as the number of connected pixels 
        % above the threshold) in order for a particle to be found. 
        % Particles smaller than this size will not be found (i.e. the program 
        % will assume that they are noise).
    Options.MaxParticleSize = 100; 
        % This is the maximum particle size (defined as the number of connected pixels 
        % above the threshold) in order for a particle to be considered "good". 
        % Particles larger than this size will be designated as "bad".
    Options.MaxEccentricity = 0.8;
        % This is the maximum eccentricity that a particle can have in order to still 
        % be considered "good". 0 = perfect circle, 1 = straight line. If the 
        % eccentricity is too high, that may indicate that the particle being 
        % analyzed is actually two diffraction limited particles close together.

    Options.MinROISize = 5; 
    Options.MaxROISize = NaN;
        % These determine the minimum and maximum size allowed for the region 
        % of interest around each particle. In between those values, the ROI 
        % scales with the size of the particle itself (i.e. bigger particles 
        % have larger regions of interest). The ROI is a square, and the 
        % values indicate the number of pixels along one side of the square 
        % (so 5 means a 5x5 pixel ROI). The minimum size must always be 
        % specified. ?NaN? can be used to indicate no maximum size.

% ---------Image Visualization Options---------
    
    Options.MinImageShow = 90;
    Options.MaxImageShow = 300;
        % These determine the minimum and maximum intensity counts that will 
        % be used to set the contrast for the grayscale images that are displayed.
        % The minimum value will be displayed as black and the maximum value 
        % will be displayed as white.
    
% ---------Types Of Analysis Options---------
    Options.UseGaussianIntensity = 'n';
        % 'y' OR 'n'
        % Choose 'y' if you want to use a Gaussian fit to determine the local 
        % background around a viral particle and use that to determine the 
        % background subtracted intensity of the particle. Choose 'n' if 
        % you want to determine the background as an average across the entire 
        % image (i.e. the background value will be the same for all viral 
        % particles within a given frame). Default is 'n'.
    Options.GrabTotalIntensityOnly = 'n';
        % 'y' OR 'n'
        % Choose 'y' if you don't want to analyze each particle individually and 
        % only want to grab the intensity trace of the entire video. Not typical.
        
% ---------Workflow Options---------
    Options.AutoCreateLabels = 'n';
        % 'y' OR 'n'
        % Choose 'y' to automatically use information in the pathname and/or 
        % filename to create labels for the output data file and/or save folder. 
        % If so, you should modify Create_Save_Folder_And_Grab_Data_Labels.m so it 
        % extracts the information that you want.
    
    Options.UserGrabExampleTrace = 'n';
        % 'y' OR 'n'
        % Choose 'y' if you want to plot an example trace in a separate window. 
        % If you choose 'y', you should modify User_Grab_Example_Trace.m with 
        % the specific details of the trace that you want to plot and how you 
        % want to plot it.
        
    Options.FrameNumberLimit = NaN; 
        % Determines the number of frames which will be loaded and analyzed. 
        % Use 'NaN' to indicate no limit (i.e. all frames will be included).        
        % You would likely only use this option if you wanted to quickly 
        % assess how the program is running without having 
        % to wait for an entire video to load. Alternatively, you could 
        % use it to exclude frames at the end of the video from your analysis.
        

end