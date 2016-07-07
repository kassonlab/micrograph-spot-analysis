function Visualize_Find_Data(varargin)

    %Debugging
    dbstop in Visualize_Find_Data at 102
    set(0, 'DefaultAxesFontSize',15)
    
    %First, we load the .mat data files.
    if length(varargin) == 1
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed',...
            varargin{1},'Multiselect', 'on');
    elseif length(varargin) == 2
        DefaultPathname = varargin{1,1}; DataFilenames = varargin{1,2};
    else
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed', 'Multiselect', 'on');
    end
    
    close all
    
    if iscell(DataFilenames)
        NumberofFiles = length(DataFilenames);
    else
        NumberofFiles = 1;
    end
    
    for CurrentFileNumber = 1:NumberofFiles
        if iscell(DataFilenames)
            CurrDataFileName = DataFilenames{1,CurrentFileNumber};
        else
            CurrDataFileName = DataFilenames;
        end
        CurrDataFilePath = strcat(DefaultPathname,CurrDataFileName);

        InputData = open(CurrDataFilePath);
        BindingDataToSave = InputData.BindingDataToSave;
        FileNumbers(CurrentFileNumber) = CurrentFileNumber;

        CurrentListNumberVirus = [];
        CurrentListVirusPerArea = [];
        CurrentListVirusIntensity = [];
        VirusIndex = 0;
        for b = 1:length(BindingDataToSave)
%                 NumberVirusesBound = BindingDataToSave(b).TotalVirusesBound; 
            NumberVirusesBound = BindingDataToSave(b).NumberGoodViruses; 
            NumberMicronsAnalyzed = BindingDataToSave(b).NumberMicronsNotBlackOut; 
            CurrentListNumberVirus(b) = NumberVirusesBound;
            CurrentListVirusPerArea(b) = NumberVirusesBound/NumberMicronsAnalyzed;
            CurrentVirusData = BindingDataToSave(b).VirusData;
            for p = 1: length(CurrentVirusData)
                if CurrentVirusData(p).IsVirusGood == 'y'
                    VirusIndex = VirusIndex +1;
                    CurrentListVirusIntensity(VirusIndex) = CurrentVirusData(p).IntensityBackSub;
                end
%                     VirusIndex = VirusIndex +1;
%                     CurrentListVirusIntensity(VirusIndex) = CurrentVirusData(p).IntensityBackSub;
            end
        end

        DataToPlot = CurrentListVirusPerArea* 6400;

%             figure(CurrentFileNumber)
%             plot(1:length(BindingDataToSave),DataToPlot)
%             hold on
%             FigureTitle = strcat(CurrDataFileName);
%             title(FigureTitle)
%             hold off

        AverageNumberVirusCombined(CurrentFileNumber) = mean(DataToPlot);
        STDCombined(CurrentFileNumber) = std(DataToPlot);

        AverageVirusIntensityCombined(CurrentFileNumber) = mean(CurrentListVirusIntensity);
        STDIntensityCombined(CurrentFileNumber) = std(CurrentListVirusIntensity);
        LowerErrorIntensity(CurrentFileNumber) = prctile(CurrentListVirusIntensity,50) - prctile(CurrentListVirusIntensity,25);
        UpperErrorIntensity(CurrentFileNumber) = prctile(CurrentListVirusIntensity,75) - prctile(CurrentListVirusIntensity,50);
        MedianVirusIntensityCombined(CurrentFileNumber) = median(CurrentListVirusIntensity);

        CharFileName = char(CurrDataFileName);
            IndexofMarker = find(CharFileName=='-');
            FileNameWOExt = CharFileName(IndexofMarker(1)+1:IndexofMarker(2)-1);
            FileNameWOExt = FileNameWOExt(1:min(25,length(FileNameWOExt)));
            XLabels{CurrentFileNumber,1}= strcat(FileNameWOExt);

    end
            figure(CurrentFileNumber +1)
            FigureTitle = strcat('Combined Data');
            title(FigureTitle)
            hold on
            errorbar(FileNumbers,AverageNumberVirusCombined,STDCombined,'bo');
            set(gca,'XTick',1:length(FileNumbers))
            set(gca,'XTickLabel',XLabels)
            Axes = gca;
            Axes.XTickLabelRotation = -45;
%             xlabel('FileName');
            ylabel(' Number of virus bound/6400 um^2');
            set(gcf, 'Position', [696 285 560 420]);
            
            IntensityWindow = figure(47);
            FigureTitle = strcat('Combined Data');
            title(FigureTitle)
            hold on
            errorbar(FileNumbers,MedianVirusIntensityCombined,LowerErrorIntensity,UpperErrorIntensity,'mo');
            set(gca,'XTick',1:length(FileNumbers))
            set(gca,'XTickLabel',XLabels)
            Axes = gca;
            Axes.XTickLabelRotation = -45;
%             xlabel('FileName');
            ylabel('Virus Intensity (25-75 prctile)');
            set(IntensityWindow, 'Position', [33 387 560 420]);
    
    
    
    disp('Thank you.  Come Again.')