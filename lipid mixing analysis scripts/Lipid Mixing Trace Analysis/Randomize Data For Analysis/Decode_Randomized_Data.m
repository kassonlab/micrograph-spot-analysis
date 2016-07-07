function Decode_Randomized_Data(PathToCodedFolder)
    PathToCodedFolder = strcat(PathToCodedFolder,'\');
    MatFile = open(strcat(PathToCodedFolder,'DecodeInfo.mat'));
    DecodeInfo = MatFile.DecodeInfo;
    SuffixOfAnalyzedData = '_Anlyzd_UsRej.mat';
    
    for n = 1:length(DecodeInfo)
        %Store the decode information in variables.
            TrueTraceDataName = DecodeInfo(n).TrueTraceDataName;
            TrueRawDataName = DecodeInfo(n).TrueRawDataName;
            CodedTraceDataName = DecodeInfo(n).CodedTraceDataName;
            CodedRawDataName = DecodeInfo(n).CodedDataName;
            CodedAnalyzedDataName = strcat(CodedTraceDataName(1:end-4),SuffixOfAnalyzedData);
            TrueAnalyzedDataName = strcat(TrueRawDataName(1:3),SuffixOfAnalyzedData);
        
        %Rename the files to their true names.  Note that move file is a
        %bit clunky--I think it actually tries to move the file and so the
        %big image stacks take some time.
            disp(strcat('Decoding file ',num2str(n),' of ', num2str(length(DecodeInfo)), '...'))
            movefile(strcat(PathToCodedFolder,CodedTraceDataName),strcat(PathToCodedFolder,TrueTraceDataName));
            movefile(strcat(PathToCodedFolder,CodedRawDataName),strcat(PathToCodedFolder,TrueRawDataName));
            movefile(strcat(PathToCodedFolder,CodedAnalyzedDataName),strcat(PathToCodedFolder,TrueAnalyzedDataName));
        
    end

    disp('Thank you. Come again')

end