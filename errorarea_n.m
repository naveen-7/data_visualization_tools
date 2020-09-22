%% Code to plot a graph with shaded error area
% Created by Naveen on 01/19/17 at cumc



function errorarea_n(X,Y,E,n,colour,alpha,points)

% X*        : N x 1 : array of x values
% Y*        : N x 1 : array of y values
% E*        : N x 1 : array of error values
% n         : 1 x 1 : std dev limit, [Default: 1]
% colour    : 1 x 3 : RGB value of raster   [Default: black]
% alpha     : 1 x 1 : alpha value of the shade [Default: 0.2]
% points    : 1 x 1 : to show or hide the actual data points [Default: 0]


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
elseif nargin==4
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    colour      = [0 0 0];
    alpha       = 0.2;
    points      = 0;
elseif nargin==5
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    varargin{5} = colour;
    alpha       = 0.2;
    points      = 0;
elseif nargin==6
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    varargin{5} = colour;
    varargin{6} = alpha;
    points      = 0;
elseif nargin==7
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = E;
    varargin{4} = n;
    varargin{5} = colour;
    varargin{6} = alpha;
    varargin{7} = points;
else
    error('Too many inputs to the function Raster_n');
end


hold on;

clear y1 y2
X1=[X,fliplr(X)];


y1(1,:) =  Y + n*(E);
y2(1,:) =  Y - n*(E);
Y1=[y1,fliplr(y2)];

if size(X1)~=size(Y1)
    clear X1
    X1=[X;flip(X)];
end


fi = fill(X1,Y1,colour,'linestyle','none');
set(fi,'FaceAlpha',alpha);

plot(X,Y,'color',colour,'LineWidth',1);

if points==1  
 plot(X,Y,'o','MarkerSize',4,'MarkerFaceColor',colour,'color','k');       
end


end