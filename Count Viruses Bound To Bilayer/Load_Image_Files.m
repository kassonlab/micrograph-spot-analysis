function [NumberOfFiles,SaveFolderDir,StackFilenames,DefaultPathname] = Load_Image_Files(varargin)

% Output:
% SaveFolderDir = directory where the analysis files will be saved
% StackFilenames = file names of the image video stacks (should be a stack of .tif)
% DefaultPathname = directory where the image video stacks are located

InputPaths = varargin{1,1};

    %First, we load the .tif files.  Should be an image stack.  We'll also
    %set up the save folder.
    if length(InputPaths) == 1
        [StackFilenames, DefaultPathname] = uigetfile('*.tif','Select .tif files to be processed',...
            InputPaths{1},'Multiselect', 'on');
        SaveFolderDir = uigetdir(InputPaths{1},'Choose the directory where data folder will be saved');
    elseif length(InputPaths) == 2
        SaveFolderDir = InputPaths{1,2};
        [StackFilenames, DefaultPathname] = uigetfile('*.tif','Select .tif files to be processed',...
            InputPaths{1,1},'Multiselect', 'on');
    else
        [StackFilenames, DefaultPathname] = uigetfile('*.tif','Select .tif files to be processed', 'Multiselect', 'on');
        SaveFolderDir = uigetdir(DefaultPathname,'Choose the directory where data folder will be saved');
    end

    %Determine the number of files selected by the user
    if iscell(StackFilenames)
        NumberOfFiles = length(StackFilenames);
    else
        NumberOfFiles = 1;
    end

end