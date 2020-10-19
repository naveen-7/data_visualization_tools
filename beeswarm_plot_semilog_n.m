function  beeswarm_plot_semilog_n(varargin)

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
% modified by naveen on 10/25/17


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
         semilogy(scat_val,rt(i,:),'o','markersize',useMS,'color',usecolor(i,1:3));                       % make scatter plot
    else
         semilogy(scat_val(i,:),rt(i,:),'o','markersize',useMS,'color',usecolor(i,1:3));
    end
    hold on
end
% 
% % --- Resize y axis ---
% ylimits = ylim;
% ylim([nanmin(nanmin(rt))-(ylimits(2)-ylimits(1))/10 nanmax(nanmax(rt))+(ylimits(2)-ylimits(1))/10])

% --- Resize x axis ---
xlimits = xlim;
xrange = max(1-xlimits(1),xlimits(2)-size(rt,1));
if xrange<0.5
    xrange = 0.5;
end
xlim([1-xrange size(rt,1)+xrange])


end
