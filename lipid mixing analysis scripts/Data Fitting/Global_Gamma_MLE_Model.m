function SumLogLike = Global_Gamma_MLE_Model(p,xdata,ydata)
    % fit gamma with lag across multiple datasets
    SumLogLike = 0;
    for i=1:length(xdata)
      x = xdata{i};
      n = length(xdata{i})
      y = n * ydata{i};
      NumberSteps = p(1);
      Tau = p(1+i);
      LagTime = p(2+length(xdata));
    
%     x =x -LagTime;
      Prob = zeros(size(x))+1e-5;
      IndextoUse = x>LagTime;
      %Model
      Prob(IndextoUse) = gamcdf(x(IndextoUse)-LagTime,NumberSteps,Tau);
      Prob = Prob + (Prob == 0).*1e-5 - (Prob==1).*1e-5; %Ensure 0<Prob<1
    
      %Log Likelihood
      LogLike = (-1)*(y.*log(Prob)+(n-y).*log(1-Prob));
      SumLogLike = SumLogLike + sum(LogLike);
  end
end

