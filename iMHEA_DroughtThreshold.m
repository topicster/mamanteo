function [Thr_DMA,Thr_MMA,Thr_D30,Thr_FFT,P_DMA,P_MMA] = iMHEA_DroughtThreshold(Date,P,flag1,flag2,varargin)
%iMHEA Calculates different drought thresholds at daily scale using
%      variable threshold approaches (Beyene et al, 2014; Van Loon, 2013).
% [Thr_DMA,Thr_MMA,Thr_D30,Thr_FFT,P_DMA,P_MMA] =
%                               iMHEA_DroughtThreshold(Date,P,flag1,flag2)
%
% Input:
% Date    = dd/mm/yyyy hh:mm:ss [date format]
% P       = evaluated variable: precipitation, streamflow [mm, l/s, m3/s]
% flag1   = 1: smooth series with 30-day moving average; 0: use original
% flag2   = 1: graph plots; 0 or empty: NOT to graph plots
%
% Output:
% Thr_XXX = 1-year threshold (366 values) [mm]
%           Thr_DMA = Moving average of daily quantile (D_MA)
%           Thr_MMA = Moving average of monthly quantile (M_MA)
%           Thr_D30 = 30-day moving window quantile (30D)
%           Thr_FFT = Fast Fourier transform approach(D_FF)
% P_DMA   = 80th quantile daily series [mm]
% P_MMA   = 80th quantile monthly series [mm]
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in September, 2017
% Last edited in November, 2018

%% DEFINITION OF A THRESHOLD

% Check input
if nargin < 3
    flag1 = 1;
    flag2 = 0;
elseif nargin < 4
    flag2 = 0;
end

% Initiliase variables
n = 366;
int = 30;
int2 = ceil(int/2);

%% Smooth interdaily data and consider irregular timing

if flag1 == 1
    % 30-day moving average of precipitation
    nP = length(P);
    P30D = nan(size(P));
    % Extend the tails of the time series for moving average
    P = [nan(int2-1,1);P;nan(int2,1)];
    for i = 1:nP
        P30D(i) = nanmean(P(i:i+int-1));
    end
else
    P30D = P;
end

%% Define a monthly drought threshold (M_MA method)

% Aggregate data at monthly scale
P1M = iMHEA_MonthlyRain(Date,P30D);
Pdm = month(Date);
P1M(1:Pdm(1)-1) = NaN;
P1M(end-(12-Pdm(end))+1:end) = NaN;
nM = length(P1M)/12;
Pmonths = repmat((1:12)',nM,1);
P_MMA = nan(n,1);
monthday = eomday(2000,1:12)';
yearday = cumsum([1;monthday]);

for i = 1:12
    P_MMA(yearday(i):yearday(i+1)) = prctile(P1M(Pmonths==i),20)/monthday(i);
end

%% Define a daily drought threshold (D_MA and D_30 methods)
P_DMA = nan(n,1);
Thr_D30 = nan(n,1);
Pdays(:,1) = day(Date,'dayofyear');
P_ly = leapyear(year(Date));

% Consider continuity of days at the initial extreme of the year
Pdays(:,2) = Pdays(:,1)-n;
Pdays(~P_ly,2) = Pdays(~P_ly,2)+1;
Pdays(Pdays(:,1)<=60+int2,2) = Pdays(Pdays(:,1)<=60+int2,1);

for i = 1:60
    P_DMA(i) = prctile(P30D(Pdays(:,1)==i),20);
    Thr_D30(i) = prctile(P30D(and(Pdays(:,2)>=i-int2-1,Pdays(:,2)<=i+int2)),20);
end

% Consider continuity of days at the final extreme of the year
Pdays(~P_ly,1) = Pdays(~P_ly,1)+1;
Pdays(:,3) = Pdays(:,1)+n;
Pdays(~P_ly,3) = Pdays(~P_ly,3)-1;
Pdays(Pdays(:,1)>=61-int2,3) = Pdays(Pdays(:,1)>=61-int2,1);

for i = 61:n
    P_DMA(i) = prctile(P30D(Pdays(:,1)==i),20);
    Thr_D30(i) = prctile(P30D(and(Pdays(:,3)>=i-int2-1,Pdays(:,3)<=i+int2)),20);
end

%% Smooth the threshold using a 30-day moving average

% Extend the tails of the time series for moving average
P_DMA = [P_DMA(end-int2+2:end);P_DMA;P_DMA(1:int2)];
P_MMA = [P_MMA(end-int2+2:end);P_MMA;P_MMA(1:int2)];

% Use a 30-day moving average
Thr_DMA = nan(n,1);
Thr_MMA = nan(n,1);
for i = 1:n
    Thr_DMA(i) = nanmean(P_DMA(i:i+int-1));
    Thr_MMA(i) = nanmean(P_MMA(i:i+int-1));
end

% Remove the tails
P_DMA(1:int2-1) = [];
P_DMA(end-int2+1:end) = [];
P_MMA(1:int2-1) = [];
P_MMA(end-int2+1:end) = [];

%% Smooth the threshold using a fast Fourier transform

% Calculate the fast Fourier transform coefficients
f = fft(P_DMA);

% Modify the coefficients and invert de modified Fourier transform
Thr_FFTset = nan(n,floor(n/2));
for i = 1:floor(n/2)
    fmod = f;
    fmod(i+1:n-i) = 0;
    Thr_FFTset(:,i) = ifft(fmod);
end

% Remove zeros
Thr_FFTset = real(Thr_FFTset);
Thr_FFTset(Thr_FFTset<0) = 0;

% Keep only the closest threshold to Thr_D30 (use squared error)
Thr_error = sum((Thr_FFTset - Thr_D30).^2);
[~,locat] = min(Thr_error);
Thr_FFT = Thr_FFTset(:,locat);

%% Plot the resulting thresholds

if flag2 == 1
    figure
    hold on
    plot(P_DMA,'DisplayName','80^{th} quantile daily series')
    plot(P_MMA,'DisplayName','80^{th} quantile monthly series')
    plot(Thr_DMA,'DisplayName','Moving average of daily quantile (D_{MA})')
    plot(Thr_MMA,'DisplayName','Moving average of monthly quantile (M_{MA})')
    plot(Thr_D30,'DisplayName','30-day moving window quantile (D_{30})')
    plot(Thr_FFT,'DisplayName','Fast Fourier transform approach(D_{FF})')
    legend show
end