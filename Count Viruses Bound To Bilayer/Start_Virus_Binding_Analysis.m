function [] = Start_Virus_Binding_Analysis(varargin)
% - - - - - - - - - - - - - - - - - - - - -

% Input:
% Start_Virus_Binding_Analysis(), in this case the user navigates to the image 
%       stacks and also chooses the parent folder where the output analysis files will be saved
%   OR
% Start_Virus_Binding_Analysis(DefaultPath), where DefaultPath is the directory to which 
%       the user will be automatically directed to find the image 
%       stacks. After choosing the stacks, the user then chooses the 
%       parent folder where the output analysis files will be saved.
%   OR
% Start_Virus_Binding_Analysis(DefaultPath,SavePath), where DefaultPath is as above, and 
%       SavePath is the parent folder where the output analysis files will be saved

% Output:
% A .mat file is created which saves all of the variables in the current 
% workspace. The information about the number of viruses bound, together 
% with details about each virus particle, will be in the BindingDataToSave 
% structure, as defined in Find_And_Process_Virus.m.

% Note: This program has been designed to process many sets of images sequentially,
% but it has been tested with individual sets, so keep that 
% in mind if you choose to process many sets at once.

% By Bob Rawle, Kasson Lab, University of Virginia, 2016
% Published online in conjunction with:
% Rawle et al., Disentangling Viral Membrane Fusion from Receptor Binding 
% Using Synthetic DNA-Lipid Conjugates, Biophysical Journal (2016) 
% http://dx.doi.org/10.1016/j.bpj.2016.05.048
% - - - - - - - - - - - - - - - - - - - - -

%Define which options will be used
[Options] = Setup_Options();
Threshold = Options.Threshold;
close all
     
    %First, we load the .tif files.  Should be an image stack.  We'll also
    %set up the save folder.
    if length(varargin) == 1
        [StackFilenames, DefaultPath] = uigetfile('*.tif','Select .tif files to be processed',...
            varargin{1},'Multiselect', 'on');
        SavePath = uigetdir(varargin{1},'Choose the directory where data folder will be saved');
    elseif length(varargin) == 2
        SavePath = varargin{1,2};
        [StackFilenames, DefaultPath] = uigetfile('*.tif','Select .tif files to be processed',...
            varargin{1,1},'Multiselect', 'on');
    else
        [StackFilenames, DefaultPath] = uigetfile('*.tif','Select .tif files to be processed', 'Multiselect', 'on');
        SavePath = uigetdir(DefaultPath,'Choose the directory where data folder will be saved');
    end

% Automatically create folders in the save directory
    SaveDataPathname = strcat(SavePath,'Binding Analysis','/');
    mkdir(SaveDataPathname);

if iscell(StackFilenames) %This lets us know if there is more than one file
    NumberOfFiles = length(StackFilenames);
else
	NumberOfFiles = 1;
end

for i = 1:NumberOfFiles
    if iscell(StackFilenames) 
        CurrentFilename = StackFilenames{1,i};
    else
        CurrentFilename = StackFilenames;
    end
    CurrStackFilePath = strcat(DefaultPath,CurrentFilename);
    CharFileName = char(CurrentFilename);

    % Now we call the function to find the virus particles and extract
    % their fluorescence intensity
    [BindingDataToSave, OtherDataToSave] = ...
        Find_And_Process_Virus(CurrStackFilePath,Threshold,CurrentFilename,...
            i, DefaultPath,Options);

    save(strcat(SaveDataPathname,Options.DataFileLabel,'-BindingData','.mat'));
end

close all
disp('---------------------')
disp ('Thank you.  Come again.')

end