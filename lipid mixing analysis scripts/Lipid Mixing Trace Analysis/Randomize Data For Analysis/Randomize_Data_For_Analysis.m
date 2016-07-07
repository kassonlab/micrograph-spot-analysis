function Randomize_Data_For_Analysis(DefaultPathname)

%This program will scramble up all the filenames of the data and the traces
%which were pulled from the data so that analyzer bias will be minimized.
%After the data has been analyzed, everything will be unscrambled.

%Debugging
    %clear all
    
    RandFolderDir = strcat(DefaultPathname,'\Coded_Analysis\');
    mkdir(RandFolderDir);
    
    AreThereMoreDataSets = 'y';
    NumSets = 0;
    NumDataFiles = 0;
    
    while AreThereMoreDataSets == 'y'
        NumSets = NumSets + 1;
        
        %Choose trace data files
        [TraceData(NumSets).Filenames, TraceData(NumSets).Path] = uigetfile('*.mat','Select .mat files to be randomized',...
            DefaultPathname,'Multiselect', 'on');
        
        %Choose original data files
        [Data(NumSets).Filenames, Data(NumSets).Path] = uigetfile('*.tif','Select original .tif files to be randomized',...
            DefaultPathname,'Multiselect', 'on');
        
        %Choose patch that goes with original data files
        [Data(NumSets).PatchFilename, Data(NumSets).PatchPath] = uigetfile('*.tif','Select patch file corresponding to video streams',...
            Data(NumSets).Path,'Multiselect', 'off');
        
        NumDataFiles = NumDataFiles + length(TraceData(NumSets).Filenames);
        
        MoreDataQuestion = inputdlg('Are there more data sets? (y/n)', 'More Data Question', 1, {'y'}, 'on');
        AreThereMoreDataSets = char(MoreDataQuestion(1,1));
    end
    
    %Now we generate the randomized code, copy the data and data analysis files to the Random Folder and
    %re-name them to their randomized form.
    
    RandomCode = randperm(NumDataFiles);
    
    NumFilesDone = 0;
    for Set = 1:NumSets
        for FileNum = 1:length(TraceData(Set).Filenames)
            NumFilesDone = NumFilesDone + 1;
            
            %Get true names and paths to files
            TrueTraceDataFilename = TraceData(Set).Filenames(FileNum);
            TrueDataFilename = Data(Set).Filenames(FileNum);
            
            CurrTrueTraceDataFile = strcat(TraceData(Set).Path,TrueTraceDataFilename);
            CurrTrueDataFile = strcat(Data(Set).Path,TrueDataFilename);
            CurrTruePatchFile = strcat(Data(Set).PatchPath,Data(Set).PatchFilename);
            
            %Check that files are matched correctly
            TraceDataIdentifier = TrueTraceDataFilename{1}(1:3);
            DataIdentifier = TrueDataFilename{1}(1:3);
            
            if ~strcmp(TraceDataIdentifier,DataIdentifier)
                disp('Error: Files didnt match up correctly')
            end
            
            CurrRandTraceDataFile = strcat(RandFolderDir,num2str(RandomCode(NumFilesDone)),'-TraceData.mat');
            CurrRandDataFile = strcat(RandFolderDir,num2str(RandomCode(NumFilesDone)),'-RawData.tif');
            CurrRandPatchFile = strcat(RandFolderDir,num2str(RandomCode(NumFilesDone)),'-Patch.tif');
            
            %Copy the files with randomized name
            disp(strcat('Copying file ',num2str(NumFilesDone), ' of ', num2str(NumDataFiles),'....'))
            copyfile(CurrTrueTraceDataFile{1},CurrRandTraceDataFile);
            copyfile(CurrTruePatchFile,CurrRandPatchFile);
            copyfile(CurrTrueDataFile{1},CurrRandDataFile);
            
            
            %Write down the de-coded information
            DecodeInfo(NumFilesDone).TrueTraceDataName = TrueTraceDataFilename{1};
            DecodeInfo(NumFilesDone).TrueRawDataName = TrueDataFilename{1};
            DecodeInfo(NumFilesDone).CodedTraceDataName = strcat(num2str(RandomCode(NumFilesDone)),'-TraceData.mat');
            DecodeInfo(NumFilesDone).CodedDataName = strcat(num2str(RandomCode(NumFilesDone)),'-RawData.tif');
        end
    end
        
    %We also write a .mat file that contains the de-coded randomization
    save(strcat(RandFolderDir,'DecodeInfo.mat'),'DecodeInfo')
    
    disp('Thank you. Come again.')
        
