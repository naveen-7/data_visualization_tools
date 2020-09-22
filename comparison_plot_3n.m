
%% function to plot a pair of data points and do a test between them
% written by Naveen at JLG on 9/19/19

function P = comparison_plot_3n(MAT,X,COLOR,MS)

% MAT    :  N x 2  : matrix with two columns of data; can contain NaN
% X      :  1 x 2  : vector of x values for the plot
% COLOR  :  2 x 3  : matrix of colors, each row is a color
% MS     :  1 x 1  : markersize

hold on;

for i=1:size(MAT,1)
    plot([X(1) X(2)],[MAT(i,1) MAT(i,2)],'-','color',[0.7 0.7 0.7],'linewidth',0.2);
    plot([X(2) X(3)],[MAT(i,2) MAT(i,3)],'-','color',[0.7 0.7 0.7],'linewidth',0.2);
end
plot(X(1),MAT(:,1),'o','color',COLOR(1,:),'markersize',MS);
plot(X(2),MAT(:,2),'o','color',COLOR(2,:),'markersize',MS);
plot(X(3),MAT(:,3),'o','color',COLOR(3,:),'markersize',MS);

P(1) = ttest_NN(MAT(:,1),MAT(:,2)); P(2) = ttest_NN(MAT(:,3),MAT(:,2)); P(3) = ttest_NN(MAT(:,1),MAT(:,3));

text(X(1)+0.2,nanmax(MAT(:))+nanmax(nanmean(abs(diff(MAT)))),star_n(P(1)));
text(X(2)+0.2,nanmax(MAT(:))+nanmax(nanmean(abs(diff(MAT)))),star_n(P(2)));
text(X(2)-0.4,nanmax(MAT(:))+1.3*nanmax(nanmean(abs(diff(MAT)))),star_n(P(3)));

xlim([X(1)-2 X(end)+2]);

end