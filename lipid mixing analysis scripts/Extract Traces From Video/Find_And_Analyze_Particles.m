function [Results,VirusDataToSave, OtherDataToSave,Options] =...
    Find_And_Analyze_Particles(VideoFilePath,VideoFilename, ...
    FileNumber, DefaultPathname,Options)
    
    %The first image in the video is read and displayed. A 3D array
    %(VideoMatrix) is created which will contain the data for all images in
    %the video.  VideoMatrix is pre-allocated with zeros to make the 
    %for loop faster.
    
    VideoInfo = imfinfo(VideoFilePath);
    if isnan(Options.FrameNumberLimit)
        NumFrames = length(VideoInfo);
    else
        NumFrames = Options.FrameNumberLimit;
    end
    
        ImageWidth = VideoInfo.Width; %in pixels
        ImageHeight = VideoInfo.Height; %in pixels
        BitDepth = VideoInfo.BitDepth;
        VideoMatrix = zeros(ImageHeight, ImageWidth, NumFrames, 'uint16');
        
        %Create a logical matrix the same size as the video matrix.
        BWVideoMatrix = VideoMatrix > 0; 
        
        %Preallocate various vectors as well
        ThresholdToFindParticles = zeros(NumFrames,1);
        TotalVideoIntensity = zeros(NumFrames,1);
        AverageVideoIntensity = zeros(NumFrames,1);
        RoughBackground = zeros(NumFrames,1);
        
        %Set up figures
        [FigureHandles] = Setup_Figures(Options);
        
    %This for loop populates the VideoMatrix with the data from each image
    %in the video.  The 1st two dimensions are the x,y of the image plane and the 3rd 
    %dimension is the frame number.
    for b = 1:NumFrames
        CurrentFrameImage = imread(VideoFilePath,b);
        VideoMatrix(:,:,b) = CurrentFrameImage;
        
        if b == Options.FrameNumToFindParticles + Options.FindFramesToAverage + 1
            %Display the finding image, including logical image. This is done 
            %before loading the rest of the frames, just in case the threshold 
            %needs to be changed and you don't want to have to wait for it 
            %to load up all the frames every time.
                FindingImage = mean(VideoMatrix(:,:,Options.FrameNumToFindParticles - Options.FindFramesToAverage...
                    :Options.FrameNumToFindParticles + Options.FindFramesToAverage),3);
                FindingImage = uint16(FindingImage);
                CurrentRoughBackground = mean(min(FindingImage));
                FindingThreshold = (CurrentRoughBackground + Options.Threshold)/2^BitDepth;
                BinaryFindingImage = im2bw(FindingImage, FindingThreshold);
                BinaryFindingImage = bwareaopen(BinaryFindingImage , Options.MinParticleSize, 8);

            %Plot the images
                set(0,'CurrentFigure',FigureHandles.ImageWindow);
                hold off
                imshow(FindingImage, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
                hold on

                set(0,'CurrentFigure',FigureHandles.BinaryImageWindow);
                imshow(BinaryFindingImage, 'InitialMagnification', 'fit','Border','tight');
                drawnow
        end
            
        
        % For each frame, the background intensity, average intensity, and 
        % integrated intensity are calculated. The threshold for each image 
        % is also calculated (this would be used if particles were being found 
        % or tracked in each image, which is currently not being done).
        RoughBackground(b) = mean(median(CurrentFrameImage));
        TotalVideoIntensity(b) = sum(sum(CurrentFrameImage));
        AverageVideoIntensity(b) = mean(mean(CurrentFrameImage));
        ThresholdToFindParticles(b) = (RoughBackground(b) + Options.Threshold)/2^BitDepth;
        
        %We apply the threshold to create a big logical matrix
        CurrThresh = ThresholdToFindParticles(b);
        BWVideoMatrix(:,:,b) = im2bw(CurrentFrameImage, CurrThresh);
        BWVideoMatrix(:,:,b) = bwareaopen(BWVideoMatrix(:,:,b), Options.MinParticleSize, 8);
        
        %Display the progress of loading the frames
        if rem(b,20)==0
            set(0,'CurrentFigure',FigureHandles.CurrentTraceWindow);
            title(strcat('Loading Frame :', num2str(b),'/', num2str(NumFrames)));
            drawnow
        end
        
    end
    
if Options.GrabTotalIntensityOnly ~= 'y'
    % Plot a trace of the background intensity
    set(0,'CurrentFigure', FigureHandles.BackgroundTraceWindow);
    hold on
    plot(RoughBackground,'r-')
    title('Calculated Background Intensity Versus Frame Number')
        
    %Set up counters
    NumGoodParticles = 0;
    NumBadParticles = 0;
        
    % Now we find all of the particles in the finding image, decide whether 
    % they are "good" or "bad" particles, and grab the integrated intensity 
    % trace within the region of interest around the particle for the entire length 
    % of the video. All of this data is then saved.
    % NOTE: There is a for loop here that only runs for one iteration. The reason 
    % that there is a for loop is because this offers flexibility if you wish to 
    % analyze multiple frames (such as in the case if you are tracking the particles 
    % frame by frame). Obviously, in such a case you would need to considerably modify the script.
    for CurrFrameNum = Options.FrameNumToFindParticles:Options.FrameNumToFindParticles
        
        %Re-define the finding image (CurrImage). This image will be an average to boost signal-to-noise.
        %Also re-define the logical image of this averaged finding image.
        CurrImage = mean(VideoMatrix(:,:,CurrFrameNum - Options.FindFramesToAverage...
            :CurrFrameNum + Options.FindFramesToAverage),3);
        CurrImage = uint16(CurrImage);
        CurrentRoughBackground = mean(min(CurrImage));
        CurrentThreshold = (CurrentRoughBackground + Options.Threshold)/2^BitDepth;
        BinaryCurrImage = im2bw(CurrImage, CurrentThreshold);
        BinaryCurrImage = bwareaopen(BinaryCurrImage , Options.MinParticleSize, 8);

        %All of the isolated regions are "particles" and will
        %be analyzed.
            ParticleComponentArray = bwconncomp(BinaryCurrImage,8);
            ParticleProperties = regionprops(ParticleComponentArray, CurrImage, 'Centroid',...
                'Eccentricity', 'PixelValues', 'Area','PixelIdxList');
            NumberOfParticlesFound = length(ParticleProperties);
            
        %Re-plot the finding image
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
            hold off
            imshow(CurrImage, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
            title('Finding Image');
            hold on
            
            set(0,'CurrentFigure',FigureHandles.BinaryImageWindow);
            imshow(BinaryCurrImage, 'InitialMagnification', 'fit','Border','tight');
            title('Finding Image, Thresholded');
            drawnow
            
            
        %Analyze each particle region
        for n = 1:NumberOfParticlesFound
            CurrentParticleProperties = ParticleProperties(n);
            CurrVesX = round(ParticleProperties(n).Centroid(1)); 
            CurrVesY = round(ParticleProperties(n).Centroid(2));
                CurrentParticleProperties.Centroid = [CurrVesX, CurrVesY];
            
            %Apply many tests to see if Particle is good
            [IsParticleGood,  ~, CurrentParticleBox, ~, ~, ~, ~, ~,...
                ~, ~, ReasonParticleFailed] =...
                Simplified_Test_Goodness(CurrImage,CurrentParticleProperties,BitDepth,...
                CurrThresh, Options.MinParticleSize, Options.MaxEccentricity,ImageWidth,...
                ImageHeight,Options.MaxParticleSize,BinaryCurrImage,Options);
            
            if strcmp(IsParticleGood,'y')
                LineColor = 'g-';
                NumGoodParticles = NumGoodParticles + 1;
                
            elseif strcmp(IsParticleGood,'n')
                LineColor = 'r-';
                NumBadParticles = NumBadParticles + 1;
                disp(ReasonParticleFailed)    
            end
                    
            %Plot a box around the Particle. Green particles are "good" and 
            %red particles are "bad".
                CVB = CurrentParticleBox;
                BoxToPlot = [CVB.Bottom,CVB.Left;CVB.Bottom,CVB.Right;CVB.Top,CVB.Right;CVB.Top,CVB.Left;CVB.Bottom,CVB.Left];

                set(0,'CurrentFigure',FigureHandles.ImageWindow);
                plot(BoxToPlot(:,2),BoxToPlot(:,1),LineColor)
                hold on
                drawnow
                         
            % The integrated intensity trace within the region of interest around 
            % the particle is then calculated. The background subtracted intensity 
            % trace is also calculated, where the background is calculated as an 
            % average across the entire image.
                CurrentParticleCroppedVideo = VideoMatrix(...
                    CurrentParticleBox.Top:CurrentParticleBox.Bottom,...
                    CurrentParticleBox.Left:CurrentParticleBox.Right,...
                    1:NumFrames);

                CurrentTraceSumArray = sum(sum((CurrentParticleCroppedVideo)));

                %Because CurrentTraceSumArray is a 3D array, the summed data is
                %transferred to row vector (CurrentTraceIntensity), so it can be
                %plotted.
                CurrentTraceIntensity = shiftdim(CurrentTraceSumArray(1,1,:),1);

                CurrentTraceIntensityBackSub = CurrentTraceIntensity -...
                    (RoughBackground(1:NumFrames).*(CurrentParticleBox.Bottom - CurrentParticleBox.Top + 1)^2)';                 

            % ALTERNATE: Use gaussian fit to determine local background intensity
            % as an alternative method for background subtraction
            if strcmp(Options.UseGaussianIntensity,'y')
                [GaussianQuantResults]  = ...   
                    Gaussian_Quantification(NumFrames,CurrentParticleCroppedVideo,Options.FrameNumToFindParticles,... 
                    CurrentParticleBox,VideoMatrix,FigureHandles,Options);
            end

            %The FigureHandles.CurrentTraceWindow is set as the current figure.  The axes are 
            %cleared (in case there was a previous trace) and then the new trace is plotted.
                set(0,'CurrentFigure',FigureHandles.CurrentTraceWindow)
                cla
                plot(CurrentTraceIntensityBackSub,LineColor)
                hold on
                if strcmp(Options.UseGaussianIntensity,'y')
                    plot(GaussianQuantResults.TraceBackSub,LineColor)
                    plot(GaussianQuantResults.TraceGauss,'b-')
                end
                title(strcat('Particle :', num2str(n),'/', num2str(NumberOfParticlesFound)));
                drawnow

            
            %Copy the data for each particle in a structure
                VirusDataToSave(n).Trace = CurrentTraceIntensity;
                VirusDataToSave(n).Trace_BackSub = CurrentTraceIntensityBackSub;
                VirusDataToSave(n).FrameNumFound = CurrFrameNum;
                VirusDataToSave(n).Coordinates = ParticleProperties(n).Centroid;
                VirusDataToSave(n).Eccentricity = ParticleProperties(n).Eccentricity;
                VirusDataToSave(n).FullFilePath = VideoFilePath;
                VirusDataToSave(n).StreamFilename = VideoFilename;
                VirusDataToSave(n).BoxAroundVirus = CurrentParticleBox;
                VirusDataToSave(n).PHdropFrameNum = Options.PHdropFrameNum;
                VirusDataToSave(n).FocusFrameNumbers  = Options.FocusFrameNumbers;
                VirusDataToSave(n).IsVirusGood = IsParticleGood;
                VirusDataToSave(n).ReasonVirusFailed = ReasonParticleFailed;
                    
        end
       
    end
    
    % If chosen, plot a specific trace in a separate window, together 
    % with the intensity trace of a pH indicator ROI. If you decide to use 
    % this, you should modify User_Grab_Example_Trace.m as it is currently set up 
    % for the last trace that I grabbed.
    if Options.UserGrabExampleTrace == 'y'
        User_Grab_Example_Trace(FigureHandles,VirusDataToSave,VideoMatrix,RoughBackground,NumFrames);
    end
    
elseif Options.GrabTotalIntensityOnly == 'y'
    VirusDataToSave=[];
end

    % Record any other data that we wish to save
    OtherDataToSave.ThresholdsUsed = ThresholdToFindParticles;
    OtherDataToSave.RoughBackground = RoughBackground;
    OtherDataToSave.TotalVideoIntensity = TotalVideoIntensity; 
    OtherDataToSave.AverageVideoIntensity = AverageVideoIntensity; 
    OtherDataToSave.PHdropFrameNum =Options.PHdropFrameNum;    
    OtherDataToSave.FocusFrameNumbers = Options.FocusFrameNumbers; 
    OtherDataToSave.Options = Options;

    % Plot the total video intensity
    set(0,'CurrentFigure',FigureHandles.CurrentTraceWindow)
    cla
    plot(TotalVideoIntensity)
    title('Total video intensity versus frame number')
    
    % Record any results that we wish to display in the command prompt window
    Results.Filename = VideoFilename;
    Results.NumBadParticles = NumBadParticles;
    Results.NumGoodParticles = NumGoodParticles;
    
end