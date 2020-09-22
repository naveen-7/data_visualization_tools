%% This function plots a customizable PSTH using the input signal
% Created by Naveen on 06/02/15 at CUMC

function [X_bin,Y_bin] = HISTOGRAM_n(Signal,Align_time,Start_time,End_time,Bin,Colour,LW,PLOT)

% let number of trials be N

% Signal*      : N x 1 : a structure (array of trials; each cell being an array of spike timings for that trial)
% Align_time*  : N x 1 : a array of timings, each for a trial
% Start_time*  : 1 x 1 : start time of the PSTH
% End_time*    : 1 x 1 : end time of the PSTH
% Bin*         : 1 x 1 : size of bin to calculate histogram; in ms
% LW           : 1 x 1 : linewidth of the plot [Default: 1]
% Colour       : 1 x 3 : RGB value of raster   [Default: black]
% PLOT         : 1 x 1 : if 1, it plots the histogram. [Default: 1]



if nargin<5
    error('Incomplete input to the function HISTOGRAM_n');
elseif nargin==5
    varargin{1} = Signal;
    varargin{2} = Align_time;
    varargin{3} = Start_time;
    varargin{4} = End_time;
    varargin{5} = Bin;
    LW          = 1;
    Colour      = [0 0 0];
    PLOT        = 1;
elseif nargin==6
    varargin{1} = Signal;
    varargin{2} = Align_time;
    varargin{3} = Start_time;
    varargin{4} = End_time;
    varargin{5} = Bin;
    varargin{6} = LW;
    Colour      = [0 0 0];
    PLOT        = 1;
elseif nargin==7
    varargin{1} = Signal;
    varargin{2} = Align_time;
    varargin{3} = Start_time;
    varargin{4} = End_time;
    varargin{5} = Bin;
    varargin{6} = LW;
    varargin{7} = Colour;
    PLOT        = 1;
    elseif nargin==8
    varargin{1} = Signal;
    varargin{2} = Align_time;
    varargin{3} = Start_time;
    varargin{4} = End_time;
    varargin{5} = Bin;
    varargin{6} = LW;
    varargin{7} = Colour;
    varargin{8} = PLOT;
else
    error('Too many inputs to the function HISTOGRAM_n');
end


time = Start_time:End_time;
Matrix = zeros(length(Signal),End_time-Start_time+1);
NaN_count=0;

for i=1:length(Signal)
    clear Aligned_signal;
    Aligned_spikes = round( Signal{i,1}-Align_time(i,1))+abs(Start_time)+1;
    temp = Aligned_spikes(find(0<Aligned_spikes & Aligned_spikes<End_time-Start_time+1));
    Matrix(i,temp)=1;
    
    if isempty(Signal{i,1})
        Matrix(i,:) = NaN;
        NaN_count=NaN_count+1;
    end
end




for j=1:size(Matrix,2)
    Vertical_mean(j) = nansum(Matrix(:,j))/(length(Signal)-NaN_count);
end


%% Gaussian Kernel :
% Bin = 0.1*LENGTH
LENGTH = End_time-Start_time;
BIN = 1;
HIS = zeros(1,LENGTH+1);

for i=1:BIN:LENGTH
    if i+BIN<=LENGTH
        HIS(i+round(Bin/2)) = (nansum(Vertical_mean(i:i+BIN))/BIN);
    end
end

HIS = HIS*1000;





% % % 
% % % Sigma=Bin;
% % % 
% % % PI=22/7;
% % % Kernel=[-5*Sigma:5*Sigma];
% % % BinSize=length(Kernel);
% % % Half_BW=(BinSize-1)/2;
% % % Kernel=[-BinSize/2:BinSize/2];
% % % Factor=1000/(Sigma*sqrt(2*PI));
% % % Kernel=Factor*(exp(-(0.5*((Kernel./Sigma).^2))));
% % % 
% % % 
% % % Kernel=Kernel; %make Kernel to a column vector
% % % PSTH=convn(Vertical_mean,Kernel,'same');% convolving raw columns of histogram




[X_bin,Y_bin] = stairs_n(time,HIS,Bin,Colour,LW);

% PLOTTING THE HISTOGRAM --------------------------------------------

if PLOT==1
    % F = figure();
    hold on;
    clear ylim;    
%   bar(time,HIS,Bin,'EdgeColor',Colour,'FaceColor','none');
%     [X_bin,Y_bin] = stairs_n(time,HIS,Bin,Colour,LW);
    plot([0 0],ylim,'-','color',[0 0 0],'LineWidth',0.7);
    xlim([time(1) time(length(time))]);
    xlabel('Time in ms','FontSize',10);
    ylabel('Sp/s','FontSize',10);
    hold off;
    box off;
    set(gca,'fontsize',10,'LineWidth',0.2)
    set(gca,'TickDir','out');
end


if PLOT==0
    delete(gcf);
end












end