%% function to plot density plots
% written by naveen at VA on 7/6/20

function plot_density_n(vector,col)

h = histcounts(vector);
[yg,xg] = ksdensity(vector);
X = xg;
Y = yg*(max(h)/max(yg));

plot(X,Y,'color',col);
end