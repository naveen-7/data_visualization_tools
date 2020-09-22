function varargout = beeswarm_plot_n(varargin)

% beeswarm_plot: Scatter plot of data with random jitter added to the data
% for better visualisation. A notched boxplot is superimposed on the
% scatter plot. A t-test is done between adjoining groups and significant
% p-value is indicated with asterisk. 
% 
% beeswarm_plot({data1}) makes a beeswarm plot of the data in data1.
% 
% beeswarm_plot({data1 data2}) makes a beeswarm plot of data1 & data2.
% 
% beeswarm_plot(...,'color',[a b c]) makes a beeswarm_plot using the color
% mentioned in the color vector
% 
% beeswarm_plot(...,'stats','x') plots the descriptive stats of the data.
% If the string is 'mean' then mean+-sem is drawn. If the string is
% 'boxplot' then a notched box-plot is drawn. the default is box-plot.
% 
% beeswarm_plot(...,'scatter','x') plots the raw data points.
% If the string is 'random' then the raw data points are randomly
% scattered. If the string is 'lined' then the raw data points are plotted
% as a histogram.
% 
% beeswarm_plot(....'distribution','x') scatters the raw data points if
% 'scatter' is set as 'lined'. 
% If the string is 'hist' then the raw data points are randomly scattered
% in each bin of the histogram. If the string is 'leaf' then the data
% points are plotted like a leaf with points having the smallest or largest
% deviation from the median of the mean having the least shift.
% 
% 
% Input variables
% ---------------
% The data to be plotted as an cell array
% 
% 
% Optional input arguments
% ------------------------
% color: RGB values of the colors
% 
% stats: whether 'mean' or 'boxplot'
% 
% scatter: whether 'random' or 'lined'
% 
% distribution: whether 'hist' or 'leaf'
%
% MS: size of the marker
% 
% ****************************************************
% Conjured by Sumitash
% 
% 5 days after a splendid Valentine's Day, 2014
% ****************************************************
% modified by naveen on 04/06/15
% modifications:
% 1) maker size increased
% 2) bug in mean corrected
% 3) output added: mean, std, median



% MEAN   = NaN(1);
% STD    = NaN(1);
% MEDIAN = NaN(1);

% --- Input ---
rtinp = varargin{1,1};

% --- Color ---
iscolor = cellfun(@(x) ischar(x) && strcmp(x, 'color'), varargin);
if any(iscolor)
    varargin = varargin(~iscolor);
    usecolor = cell2mat((varargin(iscolor)));
else
    usecolor = cool(size(rtinp,2));
end

% --- Size ---
isMS = cellfun(@(x) ischar(x) && strcmp(x, 'MS'), varargin);
if any(isMS)
    varargin = varargin(~isMS);
    useMS = cell2mat((varargin(isMS)));
else
    useMS = 5;
end




% --- Stats type ---
isstats = cellfun(@(x) ischar(x) && strcmp(x, 'stats'), varargin);
if any(isstats)
    varargin = varargin(~isstats);
    usestats = cell2mat((varargin(isstats)));
    if strcmp(usestats,'boxplot')==0 & strcmp(usestats,'mean')==0 & strcmp(usestats,'none')==0
        warning('U doofus! Not recognisable stats method! Using boxplot! B-)');
          usestats = 'boxplot';
    end
 else
     usestats = 'none';
end

% --- Scatter type ---
isscatter = cellfun(@(x) ischar(x) && strcmp(x, 'scatter'), varargin);
if any(isscatter)
    varargin = varargin(~isscatter);
    usescatter = varargin{isscatter};
    if strcmp(usescatter,'random')==0 & strcmp(usestats,'lined')==0
        warning('U scatter-brain! Not recognisable scatter type! Using lined points! B-)');
        usecatter = 'lined';
    end
else
    usescatter = 'lined';
end 

% --- Distribution type ---
isdistb = cellfun(@(x) ischar(x) && strcmp(x, 'distribution'), varargin);
if any(isdistb)
    varargin = varargin(~isdistb);
    usedistb = varargin{isdistb};
    if strcmp(usedistb,'hist')==0 & strcmp(usedistb,'leaf')==0
        warning('U block-head! Not recognisable distribution type! Using hist! B-)');
        usedistb = 'hist';
    end
else
    usedistb = 'hist';
end 

% --- Find length of values ----
for i = 1:size(rtinp,2)
    size_rt(i) = length(cell2mat(rtinp(1,i)));
end
rt = NaN(size(rtinp,2),max(size_rt));

% --- Extract values ---
for i = 1:size(rtinp,2)
    rt(i,1:length(rtinp{1,i})) = rtinp{1,i};
end

% --- Calculate histogram ---
no_of_bins = 10;
for i = 1:size(rtinp,2)
    edges = linspace(min(min(rt)),max(max(rt)),no_of_bins);
    [n(i,:),binindx(i,:)] = histc(rt(i,:),edges);
end
while any(any(n>50)) & no_of_bins<40                                       % increase the number of bins if n>40
    no_of_bins = no_of_bins + 1;
    edges = linspace(min(min(rt)),max(max(rt)),no_of_bins);
    n = []; binindx = [];
    for i = 1:size(rtinp,2)
        [n(i,:),binindx(i,:)] = histc(rt(i,:),edges);
    end
end

% --- Find scattering position ---
if strcmp(usescatter,'lined')
    scat_dif = 0.015;
    scat_val = zeros(size(rtinp,2),length(binindx));
    for i = 1:size(rtinp,2) 
        for j = 1:size(n,2)
            if n(i,j)>0
                indx = find(binindx(i,:)==j);
                if strcmp(usedistb,'leaf')==1
                    dummy = sortrows([((rt(i,indx)-nanmedian(rt(i,indx))))' indx']);                % sort rt depending on deviation from median of the bin
                    indx = dummy(:,2)';
                    dum = -scat_dif*floor(length(indx)/2):scat_dif:scat_dif*floor(length(indx)/2);
                    dum = sort(abs(dum));
                    dum(2:2:end) = dum(2:2:end).*-1;
                    if mod(length(indx),2)==0                                                       % remove 0 to fight the demon of dimension mismatch
                        dum = dum(dum~=0);
                    end
                    scat_val(i,indx) = i + dum;
                else
                    scat_val(i,indx) = i + [0.02:0.02:n(i,j)*0.02] - ((n(i,j)+1)*0.02)/2;
                end
            end
        end
    end
end

% --- Find scattering factor ---
if strcmp(usescatter,'random')
    if max(size_rt)<50
        scat_factor = 30;
    elseif max(size_rt)>=50 & max(size_rt)<100
        scat_factor = 20;
    else
        scat_factor = 10;
    end
end

% --- Plot scatter plot ---
for i = 1:size(rtinp,2)
    if strcmp(usescatter,'random')
        scat_val = repmat(i,1,length(rt(i,:)))+randn(1,length(rt(i,:)))./scat_factor;       
        h(i) = scatter(scat_val,rt(i,:),useMS,usecolor(i,1:3));                       % make scatter plot
    else
        h(i) = scatter(scat_val(i,:),rt(i,:),useMS,usecolor(i,1:3));
    end
    hold on
    if strcmp(usestats,'mean')
%         plot([nanmin(scat_val(i)) nanmax(scat_val(i))]+0.09,repmat(nanmean(rt(i,:)),1,2),'k-','LineWidth',2);  % mean
%         plot([nanmin(scat_val(i)) nanmax(scat_val(i))]+0.09,repmat(nanmean(rt(i,:))-nanstd(rt(i,:),1)./sqrt(size_rt(i)),1,2),'k-','LineWidth',1.25); % sem
%         plot([nanmin(scat_val(i)) nanmax(scat_val(i))]+0.09,repmat(nanmean(rt(i,:))+nanstd(rt(i,:),1)./sqrt(size_rt(i)),1,2),'k-','LineWidth',1.25); % sem
        

        MEAN(i) = nanmean(rt(i,:));
        STD(i)  = nanstd(rt(i,:),1);

        plot([i-0.2 i+0.2],repmat(MEAN,1,2),'k-','LineWidth',2);  % mean
        plot([i-0.2 i+0.2],repmat(MEAN-STD./sqrt(size_rt(i)),1,2),'k-','LineWidth',1.25); % sem
        plot([i-0.2 i+0.2],repmat(MEAN+STD./sqrt(size_rt(i)),1,2),'k-','LineWidth',1.25); % sem
    elseif strcmp(usestats,'boxplot')
        % --- stats for box-plot ---
        prctile_rt = prctile(rt(i,:),[25 50 75]);
        ci_median = [prctile_rt(2)-1.7*(1.25*(prctile_rt(3)-prctile_rt(1))/(1.35*sqrt(size_rt(i)))) ...
                     prctile_rt(2)+1.7*(1.25*(prctile_rt(3)-prctile_rt(1))/(1.35*sqrt(size_rt(i))))];
        lowlim = max((prctile_rt(1)-1.5*iqr(rt(i,:))),nanmin(rt(i,:)));             % lower limit for outliers
        upplim = min((prctile_rt(3)+1.5*iqr(rt(i,:))),nanmax(rt(i,:)));             % upper limit for outliers
        
        plot([i-0.15 i+0.15],repmat(prctile_rt(2),1,2),'k-','LineWidth',3);         % median
        
        MEDIAN(i) = prctile_rt(2);
        
        
        % --- CI of median ---
        plot([i-0.2 i-0.15],[ci_median(1) prctile_rt(2)],'k-','LineWidth',1.5);     
        plot([i-0.2 i-0.15],[ci_median(2) prctile_rt(2)],'k-','LineWidth',1.5);
        plot([i+0.2 i+0.15],[ci_median(1) prctile_rt(2)],'k-','LineWidth',1.5);
        plot([i+0.2 i+0.15],[ci_median(2) prctile_rt(2)],'k-','LineWidth',1.5);
        % --- box-plot lines ---
        plot([i-0.2 i+0.2],repmat(prctile_rt(1),1,2),'k-','LineWidth',1.5);
        plot([i-0.2 i+0.2],repmat(prctile_rt(3),1,2),'k-','LineWidth',1.5);
        plot([i-0.2 i-0.2],[ci_median(2) prctile_rt(3)],'k-','LineWidth',1.5);
        plot([i+0.2 i+0.2],[ci_median(2) prctile_rt(3)],'k-','LineWidth',1.5);
        plot([i-0.2 i-0.2],[ci_median(1) prctile_rt(1)],'k-','LineWidth',1.5);
        plot([i+0.2 i+0.2],[ci_median(1) prctile_rt(1)],'k-','LineWidth',1.5);
        % --- errorbars of box-plot ---
        errorbar(i,prctile_rt(3),0,upplim-prctile_rt(3),'k-','LineWidth',1.5)
        errorbar(i,prctile_rt(1),lowlim-prctile_rt(1),0,'k-','LineWidth',1.5)
        elseif strcmp(usestats,'none')
       1;  
    end
end

% --- Put p-value significance marks ---
axes_pos = get(gca,'Position');
if size(rt,1)>2
    start_arrow = axes_pos(1) + (axes_pos(3)/(size(rt,1)-1))*(0.3 + (size(rt,1))*0.05);     % Starting point of double arrow
    add_arrow = axes_pos(1) + (axes_pos(3)/(size(rt,1)-1))*(0.7 + (size(rt,1))*0.05);       % Ending point of double arrow
else
    start_arrow = axes_pos(1) + (axes_pos(3)/(size(rt,1)-1))*0.3;     % Starting point of double arrow
    add_arrow = axes_pos(1) + (axes_pos(3)/(size(rt,1)-1))*0.7;       % Ending point of double arrow
end
for i = 1:size(rt,1)-1
    if (sum(~isnan(rt(i,:)))>4 & sum(~isnan(rt(i+1,:)))>4)
        if ((lillietest(rt(i,:))==0 & lillietest(rt(i+1,:))==0) | (sum(~isnan(rt(i,:)))>30 & sum(~isnan(rt(i+1,:)))>30))
            [~,p] = ttest2(rt(i,:),rt(i+1,:));
        else
            [p,~] = ranksum(rt(i,:),rt(i+1,:));
        end
    else
        [p,~] = ranksum(rt(i,:),rt(i+1,:));
    end
    if p < 0.05 & p >= 0.01
        annotation('textbox',[mean([start_arrow add_arrow])-0.02 axes_pos(2)/3+axes_pos(4) 0.1 0.05],'String','*','FontSize',14,'FontWeight','demi','EdgeColor','none');
    elseif p < 0.01 & p >= 0.001
        annotation('textbox',[mean([start_arrow add_arrow])-0.025 axes_pos(2)/3+axes_pos(4) 0.1 0.05],'String','**','FontSize',14,'FontWeight','demi','EdgeColor','none');
    elseif p < 0.001
        annotation('textbox',[mean([start_arrow add_arrow])-0.03 axes_pos(2)/3+axes_pos(4) 0.1 0.05],'String','***','FontSize',14,'FontWeight','demi','EdgeColor','none');
    else
        annotation('textbox',[mean([start_arrow add_arrow])-0.02 axes_pos(2)/3+axes_pos(4) 0.1 0.05],'String','NS','FontSize',8,'FontWeight','demi','EdgeColor','none');
    end
    annotation('doublearrow',[start_arrow add_arrow],[axes_pos(2)/4+axes_pos(4) axes_pos(2)/4+axes_pos(4)]);
    start_arrow = start_arrow + (axes_pos(3)/(size(rt,1)));
    add_arrow = add_arrow + (axes_pos(3)/(size(rt,1)));
end

% --- Resize y axis ---
ylimits = ylim;
ylim([nanmin(nanmin(rt))-(ylimits(2)-ylimits(1))/10 nanmax(nanmax(rt))+(ylimits(2)-ylimits(1))/10])

% --- Resize x axis ---
xlimits = xlim;
xrange = max(1-xlimits(1),xlimits(2)-size(rt,1));
if xrange<0.5
    xrange = 0.5;
end
xlim([1-xrange size(rt,1)+xrange])

% --- Assign optional output ---
if nargout == 1
    varargout{1} = h;
end

if strcmp(usestats,'boxplot')
    varargout{1} = h;
    varargout{2} = MEDIAN;
end

if strcmp(usestats,'mean')
    varargout{1} = h;
    varargout{2} = MEAN;
    varargout{3} = STD;
end



end