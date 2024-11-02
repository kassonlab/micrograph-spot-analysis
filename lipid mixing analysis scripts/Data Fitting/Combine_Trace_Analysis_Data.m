  function Combine_Trace_Analysis_Data(varargin)
    Label = 'combined'
    RemoveErroneousData = 'n';
    RemoveErroneousDataFirstSet = 'n';
    %Debugging
%     dbstop in Combine_Trace_Analysis_Data at 6
    
    %First, we load the .mat data files.
    if length(varargin) == 1
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed',...
            varargin{1},'Multiselect', 'on');
    elseif length(varargin) == 2
        DefaultPathname = varargin{1,1}; DataFilenames = varargin{1,2};
    else
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed', 'Multiselect', 'on');
    end
    
    if iscell(DataFilenames) 
        NumberFiles = length(DataFilenames);
    else
        NumberFiles = 1;
    end
    
        for i = 1:NumberFiles
            if iscell(DataFilenames) 
                CurrDataFileName = DataFilenames{1,i};
            else
                CurrDataFileName = DataFilenames;
            end
            CurrDataFilePath = strcat(DefaultPathname,CurrDataFileName);
            
            InputData = open(CurrDataFilePath);
            
            if i == 1
                DataToSave = InputData.DataToSave;

                %User can remove erroneous data
                if strcmp(RemoveErroneousDataFirstSet,'y')
                    NumbertoCheck = length(DataToSave.CombinedAnalyzedTraceData);
                    for k = 1:NumbertoCheck
                        CurrentFusionData = DataToSave.CombinedAnalyzedTraceData(k).FusionData;

                        if strcmp(CurrentFusionData.Designation,'1 Fuse') && ...
                                CurrentFusionData.pHtoFusionTime(1) < 2 && ...
                                CurrentFusionData.pHtoFusionTime(1) > 0

                            DataToSave.CombinedAnalyzedTraceData(k).FusionData.Designation = ...
                                'User Erased';
                        end
                    end
                end
            else
                NewDataToAdd = InputData.DataToSave;
                if isfield(NewDataToAdd,'CombinedAnalyzedTraceData')
                    if isfield(DataToSave,'CombinedAnalyzedTraceData')
                        PreviousNumberTraces = length(DataToSave.CombinedAnalyzedTraceData);
                        NumbertoAdd = length(NewDataToAdd.CombinedAnalyzedTraceData);
                        
                        %User can remove erroneous data
                        if strcmp(RemoveErroneousData,'y')
                            if i == 2
                                for k = 1:NumbertoAdd
                                    CurrentFusionData = NewDataToAdd.CombinedAnalyzedTraceData(k).FusionData;

                                    if strcmp(CurrentFusionData.Designation,'1 Fuse') && ...
                                            CurrentFusionData.pHtoFusionTime(1) < 2 && ...
                                            CurrentFusionData.pHtoFusionTime(1) > 0

                                        NewDataToAdd.CombinedAnalyzedTraceData(k).FusionData.Designation = ...
                                            'User Erased';
                                    end
                                end
                            end
                        end

                        DataToSave.CombinedAnalyzedTraceData(PreviousNumberTraces+1:NumbertoAdd+PreviousNumberTraces)...
                            = NewDataToAdd.CombinedAnalyzedTraceData;
                    else
                        DataToSave.CombinedAnalyzedTraceData = NewDataToAdd.CombinedAnalyzedTraceData;
                    end
                    
                end
            
            end
        end
    
    save(strcat(DefaultPathname,Label,'.mat'),'DataToSave');
    
disp('Thank you.  Come Again.')
Stop_ThisWillCauseMatlabtoThrowError