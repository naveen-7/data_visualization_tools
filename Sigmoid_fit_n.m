
function Sigmoid_fit_n(xData,yData,init_params)


% init_params = [max NaN NaN]



nonnan = isnan(xData)|isnan(yData);
xData = xData(~nonnan);
yData = yData(~nonnan);


xData(1:length(xData),1)=xData;
yData(1:length(yData),1)=yData;

% % x(1) = C;  % Max value
% % x(2) = A;  % point of inflexion
% % x(3) = B;  % slope parameter


% C/(1+(exp(-B*(xtry-A))));
F = @(x,xdata)x(1)/(1+(exp(-x(3)*(xdata-x(2)))));
x0 = init_params;
x = lsqcurvefit(F,x0,xData,yData);









end