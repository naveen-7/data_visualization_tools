%% function to plot transparent scatter units
% written by naveen at JLG on 7/13/20

function scatter_n(MAT,COL,alpha)

MAT = deNaN_n(MAT,'r'); c = corr(MAT); c=c(2);
p = scatter(MAT(:,1),MAT(:,2),2,'ok');
p.MarkerFaceColor = COL;
p.MarkerEdgeColor = COL;
p.MarkerFaceAlpha = alpha;
p.MarkerEdgeAlpha = 0;
%title(strcat('corr = ',num2str(round(c,2))));

end