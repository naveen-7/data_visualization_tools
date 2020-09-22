%% Code to plot a graph with shaded error area
% Created by Naveen on 01/19/17 at cumc



function errorline_n(X,Y,E,n,colour,alpha,points,LW,Style)

% X*        : N x 1 : array of x values
% Y*        : N x 1 : array of y values
% E*        : N x 1 : array of error values
% n         : 1 x 1 : std dev limit, [Default: 1]
% colour    : 1 x 3 : RGB value of raster   [Default: black]
% alpha     : 1 x 1 : alpha value of the shade [Default: 0.2]
% points    : 1 x 1 : to show or hide the actual data points [Default: 0]
% LW        : 1 x 1 : linewidth [Default: 1]
% Style     : 1 x 1 : Marker style [Default '-']


if nargin<3
    error('Incomplete input to the function Raster_n');
elseif nargin==3
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    n           = 1;
    colour      = [0 0 0];
    alpha       = 0.2;
    points      = 0;
    LW          = 1;
    Style       = '-';
elseif nargin==4
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    colour      = [0 0 0];
    alpha       = 0.2;
    points      = 0;
    LW          = 1;
    Style       = '-';
elseif nargin==5
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    varargin{5} = colour;
    alpha       = 0.2;
    points      = 0;
    LW          = 1;
    Style       = '-';
elseif nargin==6
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    varargin{5} = colour;
    varargin{6} = alpha;
    points      = 0;
    LW          = 1;
    Style       = '-';
elseif nargin==7
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    varargin{5} = colour;
    varargin{6} = alpha;
    varargin{7} = points;
    LW          = 1;
    Style       = '-';
elseif nargin==8
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    varargin{5} = colour;
    varargin{6} = alpha;
    varargin{7} = points;
    varargin{8} = LW;
    Style       = '-';
elseif nargin==9
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    varargin{5} = colour;
    varargin{6} = alpha;
    varargin{7} = points;
    varargin{8} = LW;
    varargin{9} = Style;
else
    error('Too many inputs to the function errorline_n');
end


hold on;

clear y1 y2

plot(X,Y,Style,'color',colour,'LineWidth',LW);

try
    y1(1,:) =  Y + n*(E);
    y2(1,:) =  Y - n*(E);
    
    % plot(X,y1,'color',[colour 0.5],'LineWidth',0.5);
    % plot(X,y2,'color',[colour 0.5],'LineWidth',0.5);
    % plot([X(1) X(1)],[y1(1) y2(1)],'color',[colour 0.5],'LineWidth',0.5)
    % plot([X(length(X)) X(length(X))],[y1(length(X)) y2(length(X))],'color',[colour 0.5],'LineWidth',0.5)
    
    X_MERGE = [X X(length(X)) X(length(X)) flip(X) X(1) X(1)];
    Y_MERGE = [y2 y2(length(X)) y1(length(X)) flip(y1) y2(1) y1(1)];
    
    plot(X_MERGE,Y_MERGE,'color',[colour 0.3],'LineWidth',LW/3);
end



if points==1
    plot(X,Y,'o','MarkerSize',4,'MarkerFaceColor',colour,'color','k');
end
set(gca,'TickDir','out');

end