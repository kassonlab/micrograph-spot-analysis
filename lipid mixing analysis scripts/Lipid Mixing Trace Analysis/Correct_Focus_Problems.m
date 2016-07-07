function [CurrTrace] = Correct_Focus_Problems(CurrTrace,UniversalData)

if strcmp(UniversalData.FocusProblems,'y')
        for m=1:length(UniversalData.FocusFrameNumbers)
            currentfocusframenumber = UniversalData.FocusFrameNumbers(m);
            widthtoaverage = 3;
            offset = 1;
            numberpointstoreplace = 1;
            if currentfocusframenumber+ widthtoaverage + offset < length(CurrTrace)
                intafterfocus = median(CurrTrace(currentfocusframenumber+offset:currentfocusframenumber+offset+widthtoaverage));
                intbeforefocus = median(CurrTrace(currentfocusframenumber-offset-widthtoaverage:currentfocusframenumber-offset));
                diffromfocus = intafterfocus - intbeforefocus;

                for b = currentfocusframenumber-numberpointstoreplace:currentfocusframenumber+numberpointstoreplace
                    CurrTrace(b) = intbeforefocus;
                end
%                 CurrTrace(currentfocusframenumber-1:currentfocusframenumber+1) = [intbeforefocus intbeforefocus intbeforefocus];
%                 CurrTrace(currentfocusframenumber+2:end) = CurrTrace(currentfocusframenumber+2:end) - diffromfocus;
                CurrTrace(b+1:end) = CurrTrace(b+1:end) - diffromfocus;
            end
        end
end

end