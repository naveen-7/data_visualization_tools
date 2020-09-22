
%% function to plot a pair of data points and do a test between them
% written by Naveen at JLG on 9/19/19

function P = comparison_plot_n(MAT,X,COLOR,MS,MARKER)

% MAT*    :  N x 2  : matrix with two columns of data; can contain NaN
% X*      :  1 x 2  : vector of x values for the plot
% COLOR*  :  2 x 3  : matrix of colors, each row is a color
% MS*     :  1 x 1  : markersize

% MARKER = 'o';
% MARKER = 's';

hold on;
for i=1:size(MAT,1)
    plot([X(1) X(2)],[MAT(i,1) MAT(i,2)],'-','color',[0.7 0.7 0.7],'linewidth',0.2);
end
plot(X(1),MAT(:,1),MARKER,'color',COLOR(1,:),'markersize',MS);
plot(X(2),MAT(:,2),MARKER,'color',COLOR(2,:),'markersize',MS);
P = ttest_NN(MAT(:,1),MAT(:,2))
text(X(1)+0.1,nanmax(MAT(:))+nanmax(nanmean(abs(diff(MAT)))),star_n(P));
ylim([nanmin(MAT(:)) nanmax(MAT(:))+nanmax(nanmean(abs(diff(MAT))))])
end