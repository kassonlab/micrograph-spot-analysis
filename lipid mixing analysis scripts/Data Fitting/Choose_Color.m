function [CurrentColor] = Choose_Color(InputInfo)
    if isfield(InputInfo,'FileNumber')
        FileNumber = InputInfo.FileNumber;
    end
        if FileNumber == 1
            CurrentColor.DataPoints = 'bo';
            CurrentColor.FitLine = 'b--';
        elseif FileNumber == 2
            CurrentColor.DataPoints = 'mo';
            CurrentColor.FitLine = 'm--';
        elseif FileNumber == 3
            CurrentColor.DataPoints = 'go';
            CurrentColor.FitLine = 'g--';
        elseif FileNumber == 4
            CurrentColor.DataPoints = 'ko';
            CurrentColor.FitLine = 'k--';
        elseif FileNumber == 5
            CurrentColor.DataPoints = 'co';
            CurrentColor.FitLine = 'c--';
        elseif FileNumber == 6
            CurrentColor.DataPoints = 'ro';
            CurrentColor.FitLine = 'r--';
        elseif FileNumber == 7
            CurrentColor.DataPoints = 'yo';
            CurrentColor.FitLine = 'y--';
        elseif FileNumber == 8
            CurrentColor.DataPoints = 'bx';
            CurrentColor.FitLine = 'b--';
        elseif FileNumber == 9
            CurrentColor.DataPoints = 'mx';
            CurrentColor.FitLine = 'm--';
        else
            CurrentColor.DataPoints = 'bo';
            CurrentColor.FitLine = 'b-';
        end
        CurrentColor.ResidualPoints = CurrentColor.DataPoints;

    if isfield(InputInfo,'NumberFitsToPerform')
        NumberFitsToPerform = InputInfo.NumberFitsToPerform;
        if NumberFitsToPerform > 1
            if isfield(InputInfo,'FitNumber')
                FitNumber = InputInfo.FitNumber;
                
                if FitNumber == 1
            %             CurrentColor.DataPoints = 'bo';
                        CurrentColor.FitLine = 'k-';
                        CurrentColor.ResidualPoints = 'ko';
                    elseif FitNumber == 2
            %             CurrentColor.DataPoints = 'mo';
                        CurrentColor.ResidualPoints = 'mo';
                        CurrentColor.FitLine = 'm-';
                    elseif FitNumber == 3
            %             CurrentColor.DataPoints = 'go';
                        CurrentColor.FitLine = 'g-';
                        CurrentColor.ResidualPoints = 'go';
                    elseif FitNumber == 4
            %             CurrentColor.DataPoints = 'ko';
                        CurrentColor.FitLine = 'c-';
                        CurrentColor.ResidualPoints = 'co';
                end
            end
        end

    end

end