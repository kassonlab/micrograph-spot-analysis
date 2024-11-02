function [AllResults, ResultsReport] = Start_Fit_Time_Lists(Data, DataNames)
%   PK modified to take list of waiting times (or cell array of such lists) rather than data files
     % dbstop in Setup_And_Run_Fit at 195
    close all
    set(0, 'DefaultAxesFontSize',20)
    
    if iscell(Data)
        NumberDataFiles = length(Data);
    else 
        NumberDataFiles = 1;
        Data = {Data};
        DataNames = {DataNames};
    end    
    
    [Options] = Setup_Fit_Options();
    

    %Set up empty result structure and initialize figures

    AllResults = [];
    ResultsReport = [];
    ExtraSourceNumber = 0;
    [FigureHandles] = Initialize_Figures(Options);


    for FileNumber = 1:NumberDataFiles
        [AllResults,ResultsReport,FigureHandles,UsefulInfo] = Setup_And_Run_Fit(Data{FileNumber}, DataNames{FileNumber}, ...
            FileNumber, ResultsReport, FigureHandles, Options, AllResults);
    
    end
    
    %Display results and show figures
    
    NumberFitsToPerform = UsefulInfo.NumberFitsToPerform;
    for b = 1:length(ResultsReport)
        disp(ResultsReport(b))
        
        
        % Add information to legends for plots
        LegendInfoFuse1{1,b} = strcat(ResultsReport(b).Name,'; N=',num2str(ResultsReport(b).NumVirus),'; %=',num2str(ResultsReport(b).PercentFuse1*100,'%.1f'));
        LegendInfoFit{1,(NumberFitsToPerform+1)*b-NumberFitsToPerform} =...
            strcat(ResultsReport(b).Name,'; N=',num2str(ResultsReport(b).NumVirus),'; %=',num2str(ResultsReport(b).PercentFuse1*100,'%.1f'));
        if NumberFitsToPerform == 1
            k = 1;
            LegendInfoFit{1,(NumberFitsToPerform+1)*b-NumberFitsToPerform+k} =...
                strcat(ResultsReport(b).Name,'; fit=',Options.FitTypes);
            
            LegendInfoResid{1,(NumberFitsToPerform)*b-NumberFitsToPerform+k} =...
                strcat(ResultsReport(b).Name,'; fit=',Options.FitTypes);
        else
            for k = 1:NumberFitsToPerform
                LegendInfoFit{1,(NumberFitsToPerform+1)*b-NumberFitsToPerform+k} =...
                    strcat(ResultsReport(b).Name,'; fit=',Options.FitTypes{1,k});

                LegendInfoResid{1,(NumberFitsToPerform)*b-NumberFitsToPerform+k} =...
                    strcat(ResultsReport(b).Name,'; fit=',Options.FitTypes{1,k});
            end
        end
    end
    
    if strcmp(Options.RunKSTest,'y')
        [StatsMatrix] = Run_Stats_CDF(AllResults);
        StatsMatrix
    end
    
    if strcmp(Options.RunBootstrapMedian,'y')
        [BootstrapMedianMatrix] = Run_Bootstrap_Median(AllResults,Options);
        BootstrapMedianMatrix
    end
    
    [FitParamMatrix,GroupedResultsMatrix,NumberMatrices] = Display_Results_Matrix(AllResults,UsefulInfo,Options);
    for k = 1:NumberMatrices
        GroupedResultsMatrix(k).ResultsMatrix
    end
    for k = 1:length(FitParamMatrix)
        FitParamMatrix(k).ResultsMatrix
    end
    
    if  strcmp(Options.XLimitForPlot,'Max')
        for w = 1:NumberDataFiles
            MaxData(w) = max(AllResults(w).CDFData.SortedpHtoFList);
        end
        XLimit = max(MaxData);
    else
        XLimit = Options.XLimitForPlot;
    end
    
    Graph_Efficiencies(FigureHandles,ResultsReport,Options);
    
    set(0,'CurrentFigure',FigureHandles.Fuse1Wind)
        legend(LegendInfoFuse1,'Location','southeast');
        xlabel('Waiting Time (s)');
        ylabel('Prop of Lipid Mixing Events');
        ylim([0 1]);
        xlim([0 XLimit]);

    set(0,'CurrentFigure',FigureHandles.ResidualsWindow)
        legend(LegendInfoResid,'Location','southeast');
        xlabel('Waiting Time (s)');
        ylabel('Residuals');
        xlim([0 XLimit]);
        
    if strcmp(Options.RunBootstrap,'y')
        set(0,'CurrentFigure',FigureHandles.BootstrapWindow)
            xlabel('Waiting Time (s)');
            ylabel('Proportion Fused');
            ylim([0 1]);
            xlim([0 XLimit]);
    end

    figure(FigureHandles.FitWindow)
%     set(0,'CurrentFigure',FigureHandles.FitWindow)
        legend(LegendInfoFit,'Location','southeast');
        xlabel('Waiting Time (s)');
        ylabel('Proportion Fused');
        ylim([0 1]);
        xlim([0 XLimit]);
        
        shg
    
%     uistack(FigureHandles.FitWindow)
%     drawnow
%        ThisWillGiveError
%         Save_Figure(FigureHandles.Fuse1Wind,DefaultPathname)
%close all
end
    
function [AllResults,ResultsReport,FigureHandles,UsefulInfo] = Setup_And_Run_Fit(DataList,...
    CurrentName, FileNumber, ResultsReport,FigureHandles,Options,AllResults)
        
        UsefulInfo.FitTypes =Options.FitTypes;
        UsefulInfo.FitMethods = Options.FitMethods;
        UsefulInfo.Name = CurrentName;
        UsefulInfo.FileNumber = FileNumber;
        if iscell(Options.FitTypes)
            UsefulInfo.NumberFitsToPerform = length(Options.FitTypes);
        else
            UsefulInfo.NumberFitsToPerform = 1;
        end
        
        % Change color as we go along
        [CurrentColor]=Choose_Color(UsefulInfo);
        
    UsefulInfo.TimeCutoff = Options.TimeCutoffLow;

    TypeOfInputData = 'Normal CDF-Improved Analysis';

% [SortedpHtoFList,CumX,CumY,UsefulInfo] = Extract_Data(InputData,TypeOfInputData,UsefulInfo,Options,FigureHandles,FileNumber,CurrentColor);
    % snippets from ExtractData
    SortedpHtoFList = DataList;
    NumberTotalAnalyzed = NaN;
    % [CumX, CumY] = Generate_Prop_Cum(SortedpHtoFList);
    [CumY, CumX] = ecdf(SortedpHtoFList);
    CumY = CumY .* length(SortedpHtoFList);
    UsefulInfo.NumberTotalAnalyzed = NumberTotalAnalyzed;
    UsefulInfo.NumberDataPoints = length(SortedpHtoFList);
    UsefulInfo.MeanFusion1Time = mean(SortedpHtoFList);
    UsefulInfo.PercentFuse1 = length(SortedpHtoFList)/NumberTotalAnalyzed;
    UsefulInfo.PercentAnyFusion = length(SortedpHtoFList)/NumberTotalAnalyzed;


    %Compile useful information to pass along to fitting function
        CumYDecay = max(CumY)-CumY;
        CumYDecayNorm = CumYDecay/max(CumY);
        CumYNorm = CumY/max(CumY);
        
    % Plot the cumulative distribution function
    if ~isnan(CumX)
      set(0,'CurrentFigure',FigureHandles.Fuse1Wind)
      hold on
      plot(CumX,CumYNorm, CurrentColor.DataPoints);
    end

%     % Plot histograms
%     set(0,'CurrentFigure',FigureHandles.HistogramWindow)
%             xbins = 0:5:175;
%             hist(SortedpHtoFList,xbins)
%             xlabel('Waiting Time (s)');
%             ylabel('Num Fusion Events');
%             xlim([-4 200]);
%             title(TextFilenameWODot);

    %Record the CDF data
        AllResults(FileNumber).CDFData.CumX = CumX;
        AllResults(FileNumber).CDFData.CumY = CumY;
        AllResults(FileNumber).CDFData.CumYNorm = CumYNorm;
        AllResults(FileNumber).CDFData.SortedpHtoFList = SortedpHtoFList;
        AllResults(FileNumber).CDFData.Name = UsefulInfo.Name;
        AllResults(FileNumber).Name = UsefulInfo.Name;
        AllResults(FileNumber).TypeOfInputData = TypeOfInputData;
    
    % Call function to fit data using appropriate method and fit
    [AllResults,ResultsReport,FigureHandles] =...
        Fit_Data_Sorter(CumX,CumYDecayNorm,CumYNorm,FigureHandles,CurrentColor,FileNumber,ResultsReport,UsefulInfo,Options,AllResults);

    
        if strcmp(Options.ShowIntensityComparison,'y')
            [FigureHandles] = Analyze_Intensity_Data(FigureHandles,Compiled_Fuse1_Data,UsefulInfo,Options.XLimitForPlot);
        end
        
        if strcmp(TypeOfInputData,'Total Video Intensity')
            ResultsReport(FileNumber).NumVirus=NaN;
        end
    
    
%     %Run statistical tests if desired
%     if strcmp(Options.RunKSTest , 'y')
%         [ResultsReport] = Run_Stats_CDF(CDFData,ResultsReport,FileNumber);
%     end
   
    % Calculate the randomness parameter and Nmin
    [AllResults] = Calculate_Randomness_Parameter(FileNumber,AllResults);

end

function [SortedpHtoFList,CurrentColor]=Generate_Test_Data()
        mu = 30;
        NumPts = 200;
    Test_Data = exprnd(mu,1,NumPts);
    
    IdxToUse = Test_Data > 0;
    
    SortedpHtoFList = sort(Test_Data(IdxToUse));
    
    CurrentColor = 'bo';
end


function [DataFilenames,DefaultPathname] = Load_Data(varargin)
%First, we load the .mat data files.
        if length(varargin) == 1
            [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed',...
                char(varargin{1}),'Multiselect', 'on');
        elseif length(varargin) == 2
            DefaultPathname = varargin{1,1}; DataFilenames = varargin{1,2};
        else
            [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed', 'Multiselect', 'on');
        end    %[SortedpHtoFList,CurrentColor.DataPoints]=Generate_Test_Data();
end
