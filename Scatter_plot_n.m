
%% Helper function: Plots a scatter and gives a central tendency estimate


% Created by NAVEEN ON 07/05/15 at CUMC




function Central_Tend = Scatter_plot_n(X,Y,MS,Colour)

% clear F;
% F = figure;
plot(X,Y,'o','MarkerSize',MS,'color','k','MarkerFaceColor',Colour);
hold on;

if length(X(~isnan(X))) & length(Y(~isnan(Y))) >=4
    if ( (lillietest(X)==0 & lillietest(Y)==0) )
        Central_Tend = [nanmean(X) nanmean(Y)];
    else
        Central_Tend = [nanmedian(X) nanmedian(Y)];
    end
end

Central_Tend = [nanmean(X) nanmean(Y)];

plot(Central_Tend(1),Central_Tend(2),'+k','MarkerSize',9,'LineWidth',1.5);


set(gca,'FontSize',8,'LineWidth',1);

ref = refline(1,0);           % Line of Unity
set(ref,'linestyle','--','color',[0.4 0.4 0.4],'linewidth',0.7);


end