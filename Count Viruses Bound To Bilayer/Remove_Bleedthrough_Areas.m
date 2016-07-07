function [CurrentImage,BinaryCurrentImage] = Remove_Bleedthrough_Areas(CurrFrameNum,CurrentImage,BinaryCurrentImage,...
    ImageStackMatrix,Options,FigureHandles,ImageHeight,ImageWidth)

MinimumAreaSize = 8;
BilayerThreshold = 450;
DisplayFactor = 0.3;
MinimumBlackoutArea = 50;

BilayerImage = ImageStackMatrix(:,:,CurrFrameNum+1);
    BilayerMinimumShow = median(min(BilayerImage));
    BilayerMaximumShow = max(max(BilayerImage)) - DisplayFactor*max(max(BilayerImage));
    BilayerIntensityCutoff = mean(median(BilayerImage)) + BilayerThreshold;
BinaryBilayerImage = BilayerImage >BilayerIntensityCutoff;
BinaryBilayerImage = bwareaopen(BinaryBilayerImage,MinimumAreaSize,8);


%Plot the bilayer image
        set(0,'CurrentFigure',FigureHandles.BilayerWindow);
        hold off
        imshow(BilayerImage, [BilayerMinimumShow, BilayerMaximumShow], 'InitialMagnification', 'fit','Border','tight');
        hold on
        drawnow

        BrightAreaComponentArray = bwconncomp(BinaryBilayerImage,8);

    %The properties associated with each BrightArea in the binary image are
    %extracted.
        BrightAreaProperties = regionprops(BrightAreaComponentArray, BilayerImage, 'Centroid',...
            'Eccentricity', 'PixelValues', 'Area','PixelIdxList');
        NumberOfBrightAreaesFound = length(BrightAreaProperties);

    for j = 1:NumberOfBrightAreaesFound
        if BrightAreaProperties(j).Area > MinimumAreaSize
            CurrentBrightAreaArea = max(BrightAreaProperties(j).Area,MinimumBlackoutArea);
            CurrentBrightAreaCoordinates = BrightAreaProperties(j).Centroid;
                CurrentBrightAreaX = CurrentBrightAreaCoordinates(1);
                CurrentBrightAreaY = CurrentBrightAreaCoordinates(2);

            SizeOfSquareAroundCurrBrightArea = (sqrt(CurrentBrightAreaArea)*2);

                %Coordinates are set up to define the area around the BrightArea
                %of interest (i.e. the ROI)  
                CurrentBoxLeft = max([round(CurrentBrightAreaX) - round(SizeOfSquareAroundCurrBrightArea/2),...
                    1]);
                CurrentBoxRight = min([round(CurrentBrightAreaX) + round(SizeOfSquareAroundCurrBrightArea/2),...
                    ImageWidth]);
                CurrentBoxTop = max([round(CurrentBrightAreaY) - round(SizeOfSquareAroundCurrBrightArea/2),...
                    1]);
                CurrentBoxBottom = min([round(CurrentBrightAreaY) + round(SizeOfSquareAroundCurrBrightArea/2),...
                    ImageHeight]);

 
                BilayerImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;

                BinaryBilayerImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;
                
                
                CurrentImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;

                BinaryCurrentImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;
        end

    end

   %Plot the corrected image
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
            hold off
            imshow(CurrentImage, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
            hold on
            
            set(0,'CurrentFigure',FigureHandles.BinaryImageWindow);
            imshow(BinaryCurrentImage, 'InitialMagnification', 'fit','Border','tight');
            drawnow
        
    %Plot the corrected image
            set(0,'CurrentFigure',FigureHandles.BilayerWindow);
            hold off
            imshow(BilayerImage, [BilayerMinimumShow, BilayerMaximumShow], 'InitialMagnification', 'fit','Border','tight');
            hold on
        
end