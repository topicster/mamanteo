function [DDate,BQ,SQ,BFI,k] = iMHEA_BaseFlowUK(Date,Q,varargin)
%iMHEA Baseflow separation following Gustard et al (1992).
% [DDate,BQ,SQ,BFI] = iMHEA_BaseFlowUK(Date,Q,flag).
%
% Input:
% Date = dd/mm/yyyy hh:mm:ss [date format].
% Q    = Daily Discharge [l/s].
%        Time series will be added at daily timescale.
% flag1 = leave empty NOT to calculate the recession constant.
% flag2 = leave empty NOT to graph plots.
%
% Output:
% DDate = dd/mm/yyyy hh:mm:ss [date format, daily time scale].
% BQ    = Baseflow [l/s].
% SQ    = Stormflow [l/s].
% BFI   = Baseflow index [-].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in February, 2018
% Modified in February, 2018

%% INITIALISE VARIABLES

Daycheck = 5; % Continuous days for recessions.

%% PROCESS
% Average data at daily basis.
[DDate,DQ1] = iMHEA_Average(Date,Q,1440);
DDate = datenum(DDate);
% Consider only periods when data exist
DQ = DQ1;
DQ(isnan(DQ1)) = inf;

% Divide the time series into non-overlapping blocks of n days.
DI = min(DDate); % Initial date
DF = max(DDate); % Final date
nDate = (DI:Daycheck:DF);
k = length(nDate);
nQmin = NaN(size(nDate));

% Calculate the minima for each block.
for i = 2:k
    nQmin(i) = nanmin(DQ(DDate>=nDate(i-1)&DDate<nDate(i)));
end
nQmin(1) = [];
nDate(1) = [];

% Determine the turning points of the baseflow line based on the minima.
nBQ = nQmin;
for i = 2:(k-2)
    if or(0.9*nQmin(i)>nQmin(i-1),0.9*nQmin(i)>nQmin(i+1))
        nBQ(i) = NaN;
    end
end
nDate(isnan(nBQ)) = [];
nBQ(isnan(nBQ)) = [];
nBQ(end) = [];
nDate(end) = [];

% By linear interpolation of nBQ, estimate the daily baseflow BQ.
BQ = interp1(nDate,nBQ,DDate,'linear','extrap');
BQ(BQ>DQ) = DQ(BQ>DQ);
BQ(isnan(DQ1)) = NaN;
SQ = DQ1 - BQ;

%% CALCULATE THE RECESSION CONSTANT
if nargin>=3
    % Initialise variables
    lim = 0.8; % Minimum R2 for linear fit.
    n = length(DDate);
    R = zeros(n,1);
    M = zeros(n,1);
    LogBQ = log(BQ);
    
    h = waitbar(0,'Calculating recession constant...');
    % Identifying residence time and k.
    for i = 1:n
        Today = DDate(i);
        X = datenum(DDate(and(DDate>=Today,DDate<Today+Daycheck)));
        Y = LogBQ(and(DDate>=Today,DDate<Today+Daycheck));
        [R(i),M(i)] = regression(X',Y');
        waitbar(i/n)
    end
    close(h);
    R = R.^2;
    M = datenum(DDate(2)-DDate(1))*M;
    
    DateTau = DDate(and(R>=lim,M<0));
    RTau = R(and(R>=lim,M<0));
    MTau = M(and(R>=lim,M<0));
    K = exp(MTau);
    k = max(K);
end

%% CALCULATE THE BASEFLOW INDEX BFI
Vb = nansum(BQ(DDate>=nDate(1) & DDate<=nDate(end)));
Va = nansum(DQ1(DDate>=nDate(1) & DDate<=nDate(end)));
BFI = Vb/Va;
% Rescale the data
Date = datetime(Date,'ConvertFrom','datenum');
DDate = datetime(DDate,'ConvertFrom','datenum');

%% PLOT THE RESULTS
if nargin >= 4
    LogQ = log(Q);
    DateTau = datetime(DateTau,'ConvertFrom','datenum');
    
    figure
    subplot(3,1,1)
    hold on
    plot(Date,Q,DDate,DQ,DDate,BQ,DDate,SQ)
    xlabel('Date')
    ylabel('Discharge [l/s]')
    legend('Discharge','Daily Discharge','Baseflow','Stormflow','Location','NorthWest')
    box on
    
    subplot(3,1,2)
    hold on
    plot(Date,LogQ,DDate,log(DQ),DDate,LogBQ)
    xlabel('Date')
    ylabel('Log(Discharge) log[l/s]')
    legend('Discharge','Daily Discharge','Baseflow','Location','NorthWest')
    box on
    
    subplot(3,1,3)
    plot(DDate(1:n),R,DDate(1:n),M,Date,LogQ,...
        DateTau,LogBQ(and(R>=lim,M<0)),'o',...
        DateTau,RTau,DateTau,MTau);
    legend('Coeff. R^2','Regression slope','Log Discharge',...
        'Identified linear','Behavioural R^2','Behavioural Slope',...
        'Location','NorthWest')
    box on
end