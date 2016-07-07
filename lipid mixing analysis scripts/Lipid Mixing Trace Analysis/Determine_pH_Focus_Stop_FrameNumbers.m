function [FrameAllVirusStoppedBy,PHdropFrameNum,focusframenumbers,focusproblems] =...
    Determine_pH_Focus_Stop_FrameNumbers(OtherImportedData,InputTraceData,Options)

%Determine ph drop frame num
    if isfield(InputTraceData(1), 'PHdropFrameNum')
        if isnan(InputTraceData(1).PHdropFrameNum)
            PHdropFrameNum = Auto_Find_pH_FrameNumber(OtherImportedData);
        else
            PHdropFrameNum = InputTraceData(1).PHdropFrameNum;
            disp('user def phdrop')
        end
    else
        disp('error no phdrop')
        %PHdropFrameNum = [];
    end

%Determine frame number all virus particles have stopped moving
    if isfield(InputTraceData(1), 'FrameAllVirusStoppedBy')
        if isnan(InputTraceData(1).FrameAllVirusStoppedBy)
            FrameAllVirusStoppedBy = Options.FrameAllVirusStoppedBy;
        else
            FrameAllVirusStoppedBy = InputTraceData(1).FrameAllVirusStoppedBy;
        end
    else
%         FrameAllVirusStoppedBy = Options.FrameAllVirusStoppedBy;
        FrameAllVirusStoppedBy = [];
    end

%Determine focus frame nums
    if isfield(InputTraceData(1), 'focusframenumbers')
        if isnan(InputTraceData(1).focusframenumbers)
            focusproblems = 'n';
            focusframenumbers = NaN;
        else
            focusproblems = 'y';
            focusframenumbers = InputTraceData(1).focusframenumbers;
%             indextwouse =focusframenumbers  >PHdropFrameNum;
%             focusframenumbers = [ 124, 126, 127, 128, 129, 130, 232, 404]; %focusframenumbers(indextwouse);
            disp(strcat('user def focus problems, fr = ',num2str(focusframenumbers)));
        end
    else
        focusproblems = 'n';
        focusframenumbers = NaN;
    end
end

function PHdropFrameNum = Auto_Find_pH_FrameNumber(OtherImportedData)

ThresholdVector = OtherImportedData.ThresholdsUsed;

%RoughBackVector = ThresholdVector*2^16;
RoughBackVector = OtherImportedData.RoughBack_Med;
    %RoughBackVector = fliplr(RoughBackVector');
    
    RunMedHalfLength = 1; %Num of points on either side, prev value = 5
        StartIdx = RunMedHalfLength + 1;
        EndIdx = length(RoughBackVector)-RunMedHalfLength;
    TraceRunMedian = zeros(length(RoughBackVector),1);
    
    for n = StartIdx:EndIdx
        TraceRunMedian(n) = median(RoughBackVector(n-RunMedHalfLength:n+RunMedHalfLength));
    end
    for n = 1:StartIdx-1
        TraceRunMedian(n) = TraceRunMedian(StartIdx);
    end
    for n = EndIdx+1:length(RoughBackVector)
        TraceRunMedian(n) = TraceRunMedian(EndIdx);
    end
    
 TestWindow = figure(33);
 set(0,'CurrentFigure',TestWindow);
 plot(RoughBackVector);
 hold on
 plot(TraceRunMedian,'r');

GradRoughBack = gradient(TraceRunMedian);

MaxGrad = max(GradRoughBack);
MinGrad = min(GradRoughBack(1:50));

if abs(MinGrad) > 2*std(GradRoughBack);
    PHdropFrameNum = find(gradient(TraceRunMedian)==MinGrad(1));
else
    % Allow the user to choose the pH drop frame number
    Prompts = {' Enter pH drop frame number:',...
        };
    DefaultInputs = {'16',...
        };
    UserDef = inputdlg(Prompts,' Error finding pH drop automatically', 1, DefaultInputs, 'on');

    PHdropFrameNum = str2double(UserDef(1,1)); 
end

if PHdropFrameNum > 1/5 * length(RoughBackVector)
    disp(strcat('PH drop Frame High = ', num2str(PHdropFrameNum)))
end

 title(strcat('SUV Inj Frame = ', num2str(PHdropFrameNum)));
 drawnow
disp(strcat('pH drop Frame = ', num2str(PHdropFrameNum)));

end