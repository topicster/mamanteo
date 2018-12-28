function [DIndices,Thr] = iMHEA_Drought(Date,P,method,flag1,flag2,varargin)
%iMHEA Calculates drought threshold at daily scale and performes drought
%      analysis and calculates drought indices using a variable threshold
%      approach (Beyene et al, 2014; Van Loon, 2013).
% [DIndices,Thr] = iMHEA_Drought(Date,P,flag1,method,flag1,flag2)
%
% Input:
% Date    = dd/mm/yyyy hh:mm:ss [date format]
% P       = evaluated variable: precipitation, streamflow [mm, l/s, m3/s]
% method  = number from 1 to 4 specifying the method as follows:
% [default] 1: Thr_DMA = Moving average of daily quantile (D_MA)
%           2: Thr_MMA = Moving average of monthly quantile (M_MA)
%           3: Thr_D30 = 30-day moving window quantile (30D)
%           4: Thr_FFT = Fast Fourier transform approach(D_FF)
% flag1   = 1: smooth series with 30-day moving average; 0: use original
% flag2   = 1: graph plots; 0 or empty: NOT to graph plots
%
% Output:
% DIndices = Matriz of drought indices [10 elements].
%            1: Number of years
%            2: Mean annual precipitation
%            3: Standard deviation for annual precipitation
%            4: Number of droughts per year
%            5: Mean drought duration
%            6: Standard deviation of drought duration
%            7: Maximum drought duration
%            8: Mean drought deficit
%            9: Standard deviation of drought deficit
%           10: Maximum drought deficit
% Thr      = 1-year threshold (366 values) for the selected method [mm]
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in September, 2017
% Last edited in November, 2018

%% DEFINITION OF A THRESHOLD

% Check input
if nargin < 3
    method = 1;
    flag1 = 1;
    flag2 = 0;
elseif nargin < 4
    flag1 = 1;
    flag2 = 0;
elseif nargin < 5
    flag2 = 0;
end

[Thr_DMA,Thr_MMA,Thr_D30,Thr_FFT,P_DMA,P_MMA] = iMHEA_DroughtThreshold(Date,P,flag1);

%% CALCULATION OF DROUGHT INDICES

if method == 4
    Thr = Thr_FFT;
elseif method == 3
    Thr = Thr_D30;
elseif method == 2
    Thr = Thr_MMA;
else
    Thr = Thr_DMA;
end

[DIndices] = iMHEA_DroughtIndices(Date,P,Thr,flag1,flag2);

%% Plot results

if flag2 == 1
    
    % Plot thresholds
    subplot(3,3,8)
    a = plot(P_DMA,'DisplayName','80^{th} quantile daily series','LineWidth',2);
    b =plot(P_MMA,'DisplayName','80^{th} quantile monthly series','LineWidth',2);
    c = plot(Thr_DMA,'DisplayName','Moving average of daily quantile (D_{MA})','LineWidth',3);
    d = plot(Thr_MMA,'DisplayName','Moving average of monthly quantile (M_{MA})','LineWidth',3);
    e = plot(Thr_D30,'DisplayName','30-day moving window quantile (D_{30})','LineWidth',3);
    f = plot(Thr_FFT,'DisplayName','Fast Fourier transform approach(D_{FF})','LineWidth',3);
    % Legend
    legend([a,b,c,d,e,f])
end
