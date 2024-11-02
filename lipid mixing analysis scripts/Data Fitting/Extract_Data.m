function [SortedpHtoFList,CumX,CumY,UsefulInfo] = Extract_Data(InputData,TypeOfInputData,UsefulInfo,Options,FigureHandles,FileNumber,CurrentColor)

if strcmp(TypeOfInputData,'Normal CDF-Improved Analysis')
    AnalyzedTraceData = InputData.DataToSave.CombinedAnalyzedTraceData;
    NumberofTraces = length(AnalyzedTraceData);
    NumberFuse1 = 0;
    NumberFuse1ToPlot = 0;
    NumberFuse2 = 0;
    NumberNoFuse = 0;
    NumberSlow = 0;
    NumberOther = 0;
    PHtoFusionList = [];
    
    for k = 1:NumberofTraces
        CurrentFusionData = AnalyzedTraceData(k).FusionData;
        if strcmp(AnalyzedTraceData(k).ChangedByUser,'Incorrect Designation-Not Changed')
            % The designation is wrong, but has not been corrected, so we will skip it.
            disp('Some designations have been previously flagged as incorrect, but have not yet been changed, and so have been skipped.')
            
        elseif strcmp(CurrentFusionData.Designation,'Strange-Ignore')
            % This trace was flagged as hard to classify, so we will ignore it
            NumberOther = NumberOther + 1;
            
        else
            if strcmp(CurrentFusionData.Designation,'No Fusion')
                NumberNoFuse = NumberNoFuse + 1;
            elseif strcmp(CurrentFusionData.Designation,'1 Fuse')
                if strcmp(AnalyzedTraceData(k).ChangedByUser,'Reviewed By User') ||...
                        strcmp(AnalyzedTraceData(k).ChangedByUser,'Not analyzed')
                    NumberFuse1 = NumberFuse1 + 1;
                    NumberFuse1ToPlot = NumberFuse1ToPlot + 1;
                    PHtoFusionList(NumberFuse1ToPlot) = CurrentFusionData.pHtoFusionTime(1);
                else
                    % We Can't Necessarily Trust The Wait Time
                    NumberFuse1 = NumberFuse1 + 1;
                end

            elseif strcmp(CurrentFusionData.Designation,'2 Fuse')
                NumberFuse2 = NumberFuse2 + 1;
            elseif strcmp(CurrentFusionData.Designation,'Slow')
                NumberSlow = NumberSlow + 1;
            else
                NumberOther = NumberOther + 1;
            end
        end
    end
    
    IdxToUse = PHtoFusionList > Options.TimeCutoffLow; 
%     if UsefulInfo.FileNumber == 6
%        IdxToUse = PHtoFusionList > 4; 
%     end
    if ~isnan(Options.TimeCutoffHigh)
        Index2 = PHtoFusionList<Options.TimeCutoffHigh;
        IdxToUse =IdxToUse  +Index2;
%         IdxToUse =logical(IdxToUse );
        IdxToUse = IdxToUse ==2;
    end
    PHtoFusionList = PHtoFusionList(IdxToUse);
    PHtoFusionList = PHtoFusionList';
    SortedpHtoFList = sort(PHtoFusionList);

    %Change to real cumulative distribution function (not multiple
    %data points at repeated time points)
    if ~isempty(SortedpHtoFList)
        [CumX, CumY] = Generate_Prop_Cum(SortedpHtoFList);
    else
        CumX = 0; 
        CumY = 0;
    end
    
    NumberTotalAnalyzed = NumberFuse1 + NumberFuse2 + NumberNoFuse + NumberSlow;
        
    UsefulInfo.NumberTotalAnalyzed = NumberTotalAnalyzed;
    UsefulInfo.NumberDataPoints = length(SortedpHtoFList);
    UsefulInfo.MeanFusion1Time = mean(SortedpHtoFList);
    UsefulInfo.PercentFuse1 = NumberFuse1/NumberTotalAnalyzed;
    UsefulInfo.PercentAnyFusion = (NumberFuse1 + NumberFuse2 + NumberSlow)/NumberTotalAnalyzed;
    
    if strcmp(Options.ShowBeginningIntensity,'y')
        IntensityValues = zeros(1,NumberofTraces);
        for k = 1:NumberofTraces
            CurrentTrace = AnalyzedTraceData(k).TraceGradData.TraceRunMedian.Trace;
            IntensityValues(k) = mean(CurrentTrace(1:15));
        end
%         figure
%         hist(IntensityValues,20)
        
        set(0,'CurrentFigure',FigureHandles.IntensityWindow)
        hold on
        errorbar(FileNumber,median(IntensityValues),prctile(IntensityValues,25),prctile(IntensityValues,75),CurrentColor.DataPoints)
        ylabel('Virus Intensity')
        hold off
    end
    
elseif strcmp(TypeOfInputData,'Normal CDF')
    CompiledFuse1Data = InputData.Useful_Data_To_Save.Combined_Fuse1_Data; %InputData.CompiledData; %
    if isfield(InputData.Useful_Data_To_Save, 'Combined_Fuse2_Data')
        CompiledFuse2Data = InputData.Useful_Data_To_Save.Combined_Fuse2_Data; %InputData.CompiledData; %
        NumFuse2 = length(CompiledFuse2Data);
    else
        NumFuse2 = 0;
    end
    
    if isfield(InputData.Useful_Data_To_Save, 'Combined_NoFuse_Data')
        CompiledNoFuseData = InputData.Useful_Data_To_Save.Combined_NoFuse_Data;
        NumNoFuse = length(CompiledNoFuseData);
    else
        NumNoFuse = 0;
    end
     
        PHtoFusionList =  zeros(1,length(CompiledFuse1Data)); % InputData./1000;
        NumMobile = 0;
        NumFuse1 = length(CompiledFuse1Data);

        NumTotAnalyzed = NumFuse1+NumFuse2+NumNoFuse;

    for b = 1:length(CompiledFuse1Data)
        PHtoFusionList(b) =  CompiledFuse1Data(b).pHtoFuse1Time;
    end

    IdxToUse = PHtoFusionList > Options.TimeCutoffLow; 
    if ~isnan(Options.TimeCutoffHigh)
        Index2 = PHtoFusionList<Options.TimeCutoffHigh;
        IdxToUse =IdxToUse  +Index2;
        IdxToUse =(IdxToUse == 2 );
    end
        PHtoFusionList = PHtoFusionList(IdxToUse);
            PHtoFusionList = PHtoFusionList';
            SortedpHtoFList = sort(PHtoFusionList);

    %Change to real cumulative distribution function (not multiple
    %data points at repeated time points)
        [CumX, CumY] = Generate_Prop_Cum(SortedpHtoFList);
        
    UsefulInfo.NumberDataPoints = length(SortedpHtoFList);
    UsefulInfo.MeanFusion1Time = mean(SortedpHtoFList);
    UsefulInfo.PercentFuse1 = NumFuse1/(NumTotAnalyzed);
    UsefulInfo.PercentAnyFusion = (NumFuse1 + NumFuse2)/(NumTotAnalyzed);
    
        
elseif strcmp(TypeOfInputData,'Total Video Intensity')
    
    [CumX, CumY,MeanpHtoFuse1,NumberDataPoints] = Extract_Total_Intensity_Data(InputData,FileNumber,FigureHandles,CurrentColor);

    IdxToUse = CumX > Options.TimeCutoffLow; 
    if ~isnan(Options.TimeCutoffHigh)
        Index2 = CumX<Options.TimeCutoffHigh;
        IdxToUse =IdxToUse  +Index2;
        IdxToUse =(IdxToUse == 2 );
    end
    
    CumX = CumX(IdxToUse);
    CumY = CumY(IdxToUse);
    
    UsefulInfo.NumberDataPoints = NumberDataPoints;
    UsefulInfo.MeanFusion1Time = MeanpHtoFuse1;
    UsefulInfo.PercentFuse1 = NaN;
    UsefulInfo.PercentAnyFusion = NaN;
    SortedpHtoFList = [];
end

end