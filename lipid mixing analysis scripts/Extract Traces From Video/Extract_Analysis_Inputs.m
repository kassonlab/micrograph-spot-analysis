function [Options] = Extract_Analysis_Inputs(Options,CurrentFilename)
% The script will automatically extract the frame numbers for the pH drop, 
% the finding image, and any focus events from the filename being analyzed. 
% The filename should be formatted as: 'IdentifyingLabel-pH #1 find #2 foc #3,#4,?##.tif'
% where #1 is the pH drop frame number, #2 is the frame number of the finding 
% image, #3, #4? ## is a comma separated list of all of the frame numbers 
% for the focus events. For example, a typical file name might be 
% ?160414-102-pH 24 find 21 foc 365,366,367,689.tif'. If there are no focus 
% events, then the filename should read '? foc NaN.tif?

FilenameNoComma = strrep(CurrentFilename, ',', ' ');

        Key = 'foc';
        IndexOfKey = strfind(FilenameNoComma, Key);
        if ~isempty(IndexOfKey)
           Options.FocusFrameNumbers = sscanf(FilenameNoComma(IndexOfKey+length(Key)+1:end), '%i');
        end

        Key = 'pH';
        IndexOfKey = strfind(FilenameNoComma, Key);
        if ~isempty(IndexOfKey)
           Options.PHdropFrameNum = sscanf(FilenameNoComma(IndexOfKey+length(Key)+1:end), '%i');
        end

        Key = 'find';
        IndexOfKey = strfind(FilenameNoComma, Key);
        if ~isempty(IndexOfKey)
           Options.FrameNumToFindParticles = sscanf(FilenameNoComma(IndexOfKey+length(Key)+1:end), '%i');
        end

end