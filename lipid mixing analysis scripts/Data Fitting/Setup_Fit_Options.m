function [Options] = Setup_Fit_Options()

% Define the limit for the x-axis of all the plots
    Options.XLimitForPlot = 'Max';
%     Options.XLimitForPlot = 600;

% All fusion events with waiting times smaller or higher than the time cut off will
% not be included in the analysis
    Options.TimeCutoffLow = 0;
    Options.TimeCutoffHigh = NaN;
        % NaN means no cut off on the high end
        
% Choose the model for the fit
%     FitTypeRep = 'Gamma W Lag';
%     FitTypeRep = '2 Exp'; 
%     FitTypeRep = 'Gamma';
%      FitTypeRep = '2 Gamma';
    FitTypeRep = 'Gamma W Lag';

     Options.FitTypes = FitTypeRep;
%     Options.FitTypes =  {'Gamma','Gamma','Gamma','Gamma'};
%     Options.FitTypes =  {FitTypeRep,FitTypeRep,FitTypeRep};
%     Options.FitTypes =  {FitTypeRep,FitTypeRep,FitTypeRep,FitTypeRep};
%     Options.FitTypes =  {FitTypeRep,FitTypeRep,FitTypeRep,FitTypeRep,FitTypeRep};
%     Options.FitTypes =  {'Gamma','Gamma W Lag'};
%      Options.FitTypes =  {'1 Exp W Lag','Gamma'};
%      Options.FitTypes =  {'1 Exp W Lag','2 Exp','Gamma W Lag','2 Gamma'};

% Choose the fit method, either nonlinear least squares or maximum
% likelihood estimation,
%     Options.FitMethods = {'Least Squares','Max Like'};
%     Options.FitMethods = {'Max Like','Max Like'};
    Options.FitMethods = 'Max Like';
    
% Fix values in the fits
    Options.FixGammaN = 'n';
%     Options.FixedNValues = [2,3];
    Options.FixedNValues = [1,2,3];
%     Options.FixedNValues = [2,3,4,5];
%     Options.FixedNValues = [1, 2, 3, 4, 5];

    
% This will show a plot of the median intensity of the viruses at the
% beginning of the trace (± some percentile)
    Options.ShowBeginningIntensity = 'n';
    
% Add data from nontraditional source (like manually entered data from
% another paper), and define how many extra data sets will be included
    Options.AddExtraData = 'n';
    Options.NumberExtra = 2;

% Choose y if you want to scale the CDF according to the intensity jump
% upon fusion (this will only work if the intensity jump was calculated,
% and really only makes sense for fusion to tethered vesicles)
    Options.ShowIntensityComparison= 'n';
    
% Choose y if you want to run a two-sided KS test on each CDF, comparing it
% to the next CDF in the for loop
    Options.RunKSTest = 'y';
    
% Choose y if you want to visualize and run the confidence intervals from the
% bootstrap calculation for the CDF
    Options.RunBootstrap = 'y';
    Options.ConfidenceInterval = 95;
    Options.RunBootstrapMedian = 'n';
%     Options.ShowBootstrap = 'y';
    Options.NumberBootstraps = 500;
    
% Choose y if you want to calculate bootstrap errors for the fusion
% efficiencies. Otherwise, errors will be set to zero.
    Options.BootstrapEfficiencies = 'n';
    Options.EfficiencyConfidenceInterval = 95;
    Options.EfficiencyNumberBootstraps = 10000;
    
end
