function Start_Analysis_Photobleach(DefaultPathname,SaveFolderDir)

% - - - - - - - - - - - - - - - - - - - - -

% Input:

% Start_Analysis_Photobleach(DefaultPathname,SaveFolderDir), where 
%       DefaultPath is the directory to which the user will be automatically
%       directed to find the image video stacks. 
%       SavePath is the parent folder where the output analysis files will be saved

% Output:
% A .mat file is created which saves all of the variables in the current 
% workspace. This file will be the input for a trace analysis program. Within 
% this file, the intensity traces for each single molecule, 
% together with additional data for each virus, will be in the 
% ParticleData structure, as defined below.

% By Bob Rawle, Kasson Lab, University of Virginia, 2016
% Published online in conjunction with:
% Rawle et al., Disentangling Viral Membrane Fusion from Receptor Binding 
% Using Synthetic DNA-Lipid Conjugates, Biophysical Journal (2016) 
% http://dx.doi.org/10.1016/j.bpj.2016.05.048
% - - - - - - - - - - - - - - - - - - - - -

    [Filename,DefaultPathname] = Load_Data(DefaultPathname);
    Back_Sub_Method = 'Manual'; %Choose 'Gauss Fit' or 'Manual'--only used for back sub of bleaching stream
    
    ImageDisplayLow = 90;
    ImageDisplayHigh = 300;
    NumParticlesToChoose = 30;
    NumBackToChoose = 15;
    
    Datalabel=  DefaultPathname;
            IndexofSlash = find(Datalabel=='/');
            DataLabelForSaveFolder = Datalabel(IndexofSlash(end- 2) : IndexofSlash(end- 1));
    
            FileFolderInfo = Datalabel(IndexofSlash(end-1):end);
            IndexofSemi = find(FileFolderInfo == ';');
            DataFileLabel = FileFolderInfo(3:IndexofSemi(1)-1);
            
            IndexofDash = find(Filename == '-');
            DataFileLabel = strcat(Filename(IndexofDash(1)+1:IndexofDash(2)-1),DataFileLabel,'-');
            
    DataFolderName = strcat(DataLabelForSaveFolder,'PhotoBleachTraces/');
    SaveDataPathname = strcat(SaveFolderDir,DataFolderName);
    mkdir(SaveDataPathname);

    RadiusOfAnalysis = 6; 
    GenericBox = [-1, -1; -1, 1; 1,1; 1,-1; -1,-1];
    Box_Bleach = GenericBox.*RadiusOfAnalysis;
    
    
    %-----------------------Open 1st Frame Of Bleaching Stream-------------    
    %Only the first stream is opened b/c it takes a while to load the whole
    %stream, so we'll do that later after we've manually picked out the
    %points we want.

    BleachingStreamFilePath = strcat(DefaultPathname,Filename);
    StreamInfo = imfinfo(BleachingStreamFilePath);
    NumFrames = length(StreamInfo);
    
        Width_Bl = StreamInfo(1,1).Width;
        Height_Bl = StreamInfo(1,1).Height;
        
    BleachingStream = zeros(Height_Bl,Width_Bl,NumFrames,'uint16');
    
    for i = 1:1
        BleachingStream(:,:,i) = imread(BleachingStreamFilePath,i);
    end
    
        %Set up a figure with the first image
        StreamFirstImageWindow = figure(4);
        set(StreamFirstImageWindow,'Position', [300 300 450 450]);
        
        imshow(BleachingStream(:,:,1), [ImageDisplayLow, ImageDisplayHigh],...
            'InitialMagnification', 'fit');
        hold on
        drawnow
    
    %---------------------------------------------------------------------- 
    
    
    
    %-----------------------Finding Vesicle Regions------------------------
        %The vesicles are selected manually in the bleaching stream--this
        %is to look at the dim vesicles that might be more reliable in
        %terms of single-step counting
        
        for i = 1:NumParticlesToChoose
                set(0,'CurrentFigure',StreamFirstImageWindow)
                ImageTitle = strcat('Pick Particles--num left to choose = ',num2str(NumParticlesToChoose-i+1));
                title(ImageTitle)
                [ginx,giny] = ginput(1);

                
                Centroid_Bleach(i,1:2) = [round(ginx),round(giny)];
                    plot(Box_Bleach(:,1) + round(ginx), Box_Bleach(:,2) + round(giny), 'r-')
                    hold on

                ParticleData(i).Centroid = Centroid_Bleach(i,1:2);
         end
        
    %----------------------------------------------------------------------

    
    %-----------------------Manually Choose Background Regions-------------
    if strcmp(Back_Sub_Method,'Manual')
        
        for i = 1:NumBackToChoose
            set(0,'CurrentFigure',StreamFirstImageWindow)
            ImageTitle = strcat('Now pick background regions--regions left to choose = ',num2str(NumBackToChoose-i+1));
            title(ImageTitle)
            [ginx,giny] = ginput(1);

            BackCoordinates(i,1:2) = [round(ginx),round(giny)];
                plot(Box_Bleach(:,1) + round(ginx), Box_Bleach(:,2) + round(giny), 'c-')
                hold on
        end
        BackgroundData.BackgroundCoordinates = BackCoordinates(i,1:2);
    end
    %----------------------------------------------------------------------
    
    %------------------Opening the Bleaching Stream------------------------

    %Now we open the rest.
    
    for i = 1:NumFrames
        BleachingStream(:,:,i) = imread(BleachingStreamFilePath,i);
    end

    %----------------------------------------------------------------------
        
    
    BlankTrace(1:NumFrames) = 0;

    %Define the edges of the calibration image
    Left_edge = 1; Right_edge = size(BleachingStream(:,:,1),2);
    Top_edge = 1; Bot_edge = size(BleachingStream(:,:,1),1);
    NumOfVesicleRegions = length(ParticleData);
    
    for n = 1:NumOfVesicleRegions
       CurrVesX = ParticleData(n).Centroid(1);
       CurrVesY = ParticleData(n).Centroid(2);

       if (CurrVesX <= Right_edge - RadiusOfAnalysis+1 && ...
               CurrVesX >= Left_edge + RadiusOfAnalysis + 1 && ...
               CurrVesY <= Bot_edge - RadiusOfAnalysis+1 && ...
               CurrVesY >= Top_edge + RadiusOfAnalysis + 1)

           %The type of quantification depends on the method of background
           %subtraction.
           if strcmp(Back_Sub_Method,'Manual')
           %---------------------------Manual----------------------------------
               %-------Determine which background area is the closest--------
               BackTest = ((BackCoordinates(:,1)-CurrVesX).^2 + (BackCoordinates(:,2)-CurrVesY).^2).^5; 
               IdxClosest = BackTest==min(BackTest);
%                ClosestBackValue_Calib = BackInt_Calib(IdxClosest); 
               ClosestBackCoords_Bleach = BackCoordinates(IdxClosest,1:2);
               %-------------------------------------------------------------

%                 
%                 ParticleData(n).CalibRawIntensity = VesRawIntensity;
%                 ParticleData(n).CalibBackSubInt = VesBackSubInt;


               %Then do the bleaching trace analysis
               [RawCurrTrace, BackSubTrace] = Calc_Trace_Manual_Back(CurrVesX, CurrVesY, BleachingStream,BlankTrace,...
                   StreamFirstImageWindow,Box_Bleach,RadiusOfAnalysis,ClosestBackCoords_Bleach);

%                figure(3)
%                plot(RawCurrTrace,'b-')

               figure(5)
               plot(BackSubTrace,'r-')
               drawnow

               ParticleData(n).RawBleachTrace = RawCurrTrace;
               ParticleData(n).BackSubBleachTrace = BackSubTrace;
               ParticleData(n).CoordsOfBackUsed_Bleach = ClosestBackCoords_Bleach;
               ParticleData(n).Back_Sub_Method = Back_Sub_Method;

           %-------------------------------------------------------------------



           elseif strcmp(Back_Sub_Method,'Gauss Fit')
           %---------------------------Gauss Fit-------------------------------
               set(0,'CurrentFigure',InitImageWindow)
                    plot(BoxCalib(:,1) + CurrVesX, BoxCalib(:,2) + CurrVesY, 'g-')
                    hold on

               VesArea = InitialImage(round(CurrVesY)-RoA_Calib:round(CurrVesY)+RoA_Calib,...
                   round(CurrVesX)-RoA_Calib:round(CurrVesX)+RoA_Calib);

               VesRawIntensity = sum(sum(VesArea));
               %Do Gaussian Fit for calibration image
                   try
                        [OffsetFrom2DFit, Noise] = Vesicle_Gaussian_Fit(VesArea);
                        VesBackSubInt = sum(sum(VesArea - (abs(OffsetFrom2DFit) + Noise)));
                   catch
                        disp('Gauss Fit Failed')
                        VesBackSubInt = 0;
                   end
                   
               ParticleData(n).CalibRawIntensity = VesRawIntensity;
               ParticleData(n).CalibBackSubInt = VesBackSubInt;


               [RawCurrTrace, BackSubTrace] = Calc_Trace_Gauss_Fit(CurrVesX, CurrVesY, BleachingStream,BlankTrace,...
                   StreamFirstImageWindow,Box_Bleach,RadiusOfAnalysis);

%                figure(3)
%                plot(RawCurrTrace,'b-')

               figure(5)
               plot(BackSubTrace,'r-')
               drawnow

               ParticleData(n).RawBleachTrace = RawCurrTrace;
               ParticleData(n).BackSubBleachTrace = BackSubTrace;
               %PhotobleachTraceData(n).CoordsOfBackUsed_Bleach = ClosestBackCoords_Bleach;
               ParticleData(n).Back_Sub_Method = Back_Sub_Method;

               save(strcat(SaveDataPathname,FileLabel,'StructOfVesiclesProps_ManPick_GaussBack.mat'),...
                   'PhotobleachTraceData');
           %-------------------------------------------------------------------
           end
       end 

    end
    
    PhotobleachTraceData.ParticleData = ParticleData;
    PhotobleachTraceData.BackgroundData = BackgroundData;
    save(strcat(SaveDataPathname,DataFileLabel,'Traces.mat'),...
                   'PhotobleachTraceData');
               
ThisWillGiveErrorToStopFunction


end

function [RawCurrTrace, BackSubTrace] = Calc_Trace_Manual_Back(CurrVesX, CurrVesY, BleachingStream, BlankTrace,...
    StreamFirstImageWindow,Box_Bleach,RadiusOfAnalysis,ClosestBackCoords_Bleach)

    RawCurrTrace = BlankTrace;
    BackSubTrace = BlankTrace;
    CurrBackX = ClosestBackCoords_Bleach(1,1);
    CurrBackY = ClosestBackCoords_Bleach(1,2);
           
        CurrVesX = round(CurrVesX);
        CurrVesY = round(CurrVesY);
        
        set(0,'CurrentFigure',StreamFirstImageWindow)
        plot(Box_Bleach(:,1) + CurrVesX, Box_Bleach(:,2) + CurrVesY, 'g-')
        hold on

    
    VesAreaStream = BleachingStream(CurrVesY-RadiusOfAnalysis:CurrVesY+RadiusOfAnalysis,...
        CurrVesX-RadiusOfAnalysis:CurrVesX+RadiusOfAnalysis,:);
    
    BackAreaStream = BleachingStream(CurrBackY-RadiusOfAnalysis:CurrBackY+RadiusOfAnalysis,...
        CurrBackX-RadiusOfAnalysis:CurrBackX+RadiusOfAnalysis,:);
  
    Trace3DArrayRaw = sum(sum(VesAreaStream));
    RawCurrTrace(1:size(Trace3DArrayRaw,3)) = Trace3DArrayRaw(1,1,1:size(Trace3DArrayRaw,3));
    

    %Calculate running median of RawBackTrace
        TraceBackArrayRaw = sum(sum(BackAreaStream));
        RawBackTrace(1:size(TraceBackArrayRaw,3)) = TraceBackArrayRaw(1,1,1:size(TraceBackArrayRaw,3));
        Length_Median = 10;
            l_med = Length_Median;
            
            BackTrace_RunMed = zeros(size(RawBackTrace));
        for i = 1:l_med
            BackTrace_RunMed(i) = median(RawBackTrace(1:i+l_med));
        end
        
        for i = l_med+1:length(RawBackTrace)-l_med
            BackTrace_RunMed(i) = median(RawBackTrace(i-l_med:i+l_med));
        end
        
        for i = length(RawBackTrace)-l_med+1:length(RawBackTrace)
            BackTrace_RunMed(i) = median(RawBackTrace(i-l_med:length(RawBackTrace)));
        end
        
        %Diagnostic
%             figure(7)
%                 plot(RawBackTrace,'g-');
%                 hold on
%                 plot(BackTrace_RunMed,'k-');
%                 drawnow
%                 hold off
                
    %Do background subtraction
        %Trace3DArrayBackSub = Trace3DArrayRaw - sum(sum(BackAreaStream));
        %BackSubTrace(1:size(Trace3DArrayRaw,3)) = Trace3DArrayBackSub(1,1,1:size(Trace3DArrayRaw,3));
        BackSubTrace = RawCurrTrace - BackTrace_RunMed;
end

function [RawCurrTrace, BackSubTrace] = Calc_Trace_Gauss_Fit(CurrVesX, CurrVesY, BleachingStream, BlankTrace,...
    StreamFirstImageWindow,Box_Bleach,RadiusOfAnalysis)


    RawCurrTrace = BlankTrace;
    BackSubTrace = BlankTrace;
               
    %Translate centroid into the 2x2 bin coords
        CurrVesX = round(CurrVesX/2);
        CurrVesY = round(CurrVesY/2);
        
        set(0,'CurrentFigure',StreamFirstImageWindow)
        plot(Box_Bleach(:,1) + CurrVesX, Box_Bleach(:,2) + CurrVesY, 'g-')
        hold on

    
    VesAreaStream = BleachingStream(CurrVesY-RadiusOfAnalysis:CurrVesY+RadiusOfAnalysis,...
        CurrVesX-RadiusOfAnalysis:CurrVesX+RadiusOfAnalysis,:);
    
    Trace3DArrayRaw = sum(sum(VesAreaStream));
    RawCurrTrace(1:size(Trace3DArrayRaw,3)) = Trace3DArrayRaw(1,1,1:size(Trace3DArrayRaw,3));
    
    %VesArea_CurrFrame = zeros(size(VesAreaStream,1),size(VesAreaStream,2),'uint16');
    for i = 1:size(VesAreaStream,3)
        %Fit the SUV to a Gaussian
        VesArea_CurrFrame = VesAreaStream(:,:,i);
        %Do Gaussian Fit
            try
                [OffsetFrom2DFit, Noise] = Vesicle_Gaussian_Fit(VesArea_CurrFrame);
                BackSubTrace(i) = sum(sum(VesArea_CurrFrame - (abs(OffsetFrom2DFit) + Noise)));
            catch
                %disp('Gauss Fit Failed')
                if i == 1
                    break
                else
                    BackSubTrace(i) = BackSubTrace(i-1);
                end
            end
            
            %If the fit fails 10 frames in a row, then don't bother fitting
            %the rest
            if (BackSubTrace(i) == BackSubTrace(i-1) && i > 10)
                if (sum(BackSubTrace(i-10:i))==BackSubTrace(i)*11)
                    break
                end
            end
    end


    
end



function [DataFilenames,DefaultPathname] = Load_Data(varargin)
        if length(varargin) == 1
            [DataFilenames, DefaultPathname] = uigetfile('*.tif','Select .tif files to be analyzed',...
                char(varargin{1}),'Multiselect', 'off');
        end    
end