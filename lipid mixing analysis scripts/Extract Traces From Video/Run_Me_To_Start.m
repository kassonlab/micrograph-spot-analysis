function [] = Run_Me_To_Start(varargin)

% - - - - - - - - - - - - - - - - - - - - -

% Input:
% Run_Me_To_Start(), in this case the user navigates to the image video
%       stacks and also chooses the parent folder where the output analysis files will be saved
%   OR
% Run_Me_To_Start(DefaultPath), where DefaultPath is the directory to which 
%       the user will be automatically directed to find the image video
%       stacks. After choosing the video stacks, the user then chooses the 
%       parent folder where the output analysis files will be saved.
%   OR
% Run_Me_To_Start(DefaultPath,SavePath), where DefaultPath is as above, and 
%       SavePath is the parent folder where the output analysis files will be saved

% Output:
% A .mat file is created which saves all of the variables in the current 
% workspace. This file will be the input for the Lipid Mixing Trace Analysis program. Within 
% this file, the intensity traces for each viral particle which has been 
% found, together with additional data for each virus, will be in the 
% VirusDataToSave structure, as defined in Find_And_Analyze_Particles.m

% Note: This program has been designed to process many videos sequentially,
% but it has been tested with individual video streams, so keep that 
% in mind if you choose to process many videos at once.

% By Bob Rawle, Kasson Lab, University of Virginia, 2016
% Published online in conjunction with:
% Rawle et al., Disentangling Viral Membrane Fusion from Receptor Binding 
% Using Synthetic DNA-Lipid Conjugates, Biophysical Journal (2016) 
% http://dx.doi.org/10.1016/j.bpj.2016.05.048
% - - - - - - - - - - - - - - - - - - - - -

%Define which options will be used
[Options] = Setup_Options();

close all

%Load the image files, chosen by the user
[NumberOfFiles,SaveParentFolder,StackFilenames,DefaultPathname] = Load_Image_Files(varargin);

% If selected, info is automatically grabbed from the data filenames and/or pathnames to make more 
% informative save folder directory and output analysis filenames. The save 
% folder is then created inside the parent directory.
if strcmp(Options.AutoCreateLabels,'y')
    [DataFileLabel,SaveDataPathname] = Create_Save_Folder_And_Grab_Data_Labels(DefaultPathname,...
    SaveParentFolder,Options);
else
    % Otherwise, the label and save folder are defined as below.
    DataFileLabel = 'TestLabel';
    SaveDataPathname = SaveParentFolder;
    mkdir(SaveDataPathname);
end

% Analyze each video stream one by one
for i = 1:NumberOfFiles
    
    if NumberOfFiles > 1
        CurrentFilename = StackFilenames{1,i};
    else
        CurrentFilename = StackFilenames;
    end
    
    CurrStackFilePath = strcat(DefaultPathname,CurrentFilename);

    % Extract focus frame numbers, pH drop frame number, and frame to find
    % the viruses from the data filename if it is there
        if strcmp(Options.ExtractInputsFromFilename,'y')
            [Options] = Extract_Analysis_Inputs(Options,CurrentFilename);
        end
    % Print out options to commandline
        Options     
        
    % Now we call the function to find the virus particles and extract
    % their fluorescence intensity traces
    [Results,VirusDataToSave, OtherDataToSave,Options] =...
        Find_And_Analyze_Particles(CurrStackFilePath,CurrentFilename, ...
            i, DefaultPathname,Options);

    % Analysis output file is saved to the save folder. All variables are saved.
    save(strcat(SaveDataPathname,DataFileLabel,'-Traces','.mat'));
    
    % Results are displayed in the command prompt window
    Results
end

disp('---------------------')
disp ('Thank you.  Come again.')

ThisWillCauseMatlabtoThrowError
end