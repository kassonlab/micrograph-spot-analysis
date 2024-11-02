function [CurrTrace] = Correct_Focus_Problems(CurrTrace,UniversalData)
% Focus problems are corrected by determining the difference in intensity 
% before and after the focus frame number (taking the median over some 
% number of frames before and after the focus frame number). This 
% difference is then subtracted from the intensity trace for all frame 
% numbers after the focus problem (essentially erasing the focus problem 
% from the trace which will be used to perform the analysis). The actual 
% frame number where the focus problem occurred is set to the median 
% intensity value before the focus event.

if strcmp(UniversalData.FocusProblems,'y')
        for m=1:length(UniversalData.FocusFrameNumbers)
            currentfocusframenumber = UniversalData.FocusFrameNumbers(m);
            widthtoaverage = 3;
            offset = 1;
            numberpointseithersidetoreplace = 1;
            if currentfocusframenumber+ widthtoaverage + offset < length(CurrTrace)
                intafterfocus = median(CurrTrace(currentfocusframenumber+offset:currentfocusframenumber+offset+widthtoaverage));
                intbeforefocus = median(CurrTrace(currentfocusframenumber-offset-widthtoaverage:currentfocusframenumber-offset));
                diffromfocus = intafterfocus - intbeforefocus;

                for b = currentfocusframenumber-numberpointseithersidetoreplace:currentfocusframenumber+numberpointseithersidetoreplace
                    CurrTrace(b) = intbeforefocus;
                end
                CurrTrace(b+1:end) = CurrTrace(b+1:end) - diffromfocus;
            end
        end
end

end