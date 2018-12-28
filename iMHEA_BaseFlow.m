function [BQ2,SQ2,BFI2,k] = iMHEA_BaseFlow(Date,Q,varargin)
%iMHEA Baseflow separation following Chapman (1999).
% [BQ,SQ,BFI,k] = iMHEA_BaseFlow(Date,Q,flag).
%
% Input:
% Date = dd/mm/yyyy hh:mm:ss [date format].
% Q    = Daily Discharge [l/s].
%        A regular interval in the Q time series is needed.
% flag = leave empty NOT to graph plots.
%
% Output:
% BQ  = Baseflow [l/s].
% SQ  = Stormflow [l/s].
% BFI = Baseflow Index [-].
% k   = Recession constant [-].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in June, 2014
% Modified in November, 2017

%% INITIALISE VARIABLES
n = length(Date);
Daycheck = 7; % Continuous days for recessions.
lim = 0.8; % Minimum R2 for linear fit.

% Initialise variables
R = zeros(n,1);
M = zeros(n,1);
LogQ = log(Q);

%% PROCESS
h = waitbar(0,'Calculating baseflow...');
% Identifying residence time and k.
for i = 1:n
    Today = Date(i);
    X = datenum(Date(and(Date>=Today,Date<Today+Daycheck)));
    Y = LogQ(and(Date>=Today,Date<Today+Daycheck));
    [R(i),M(i)] = regression(X',Y');
    waitbar(i/n)
end
close(h);
R = R.^2;
M = datenum(Date(2)-Date(1))*M;

DateTau = Date(and(R>=lim,M<0));
RTau = R(and(R>=lim,M<0));
MTau = M(and(R>=lim,M<0));
K = exp(MTau);
k = max(K);
% T = -1/log(K);
% mT = min(T);
% MT = max(T);
% k = 0.949 to 0.993;
% C = 1-k, or 0.018 to 0.085;
% alpha = -0.01 to -0.81;
% C = 1-k, or 0.011 to 0.197; with alpha
C = datenum(Date(2)-Date(1))*.085;
alpha = -0.1;

% Separate baseflow [BQ,SQ] = par3(Q,k,C,alpha).
if isempty(k)
    BQ2 = [];
    SQ2 = [];
    BFI2 = [];
    disp('These time series do not allow the determination of base flow.')
    return
else
    [BQ1] = par3(Q,k,1-k,0);
    [BQ2,SQ2] = par3(Q,k,C,0);
    [BQ3] = par3(Q,k,C*2,alpha);
    LogBQ1 = log(BQ1);
    LogBQ2 = log(BQ2);
    LogBQ3 = log(BQ3);
end
 
% BFI1 = sum(BQ1)/sum(Q);
BFI2 = sum(BQ2)/sum(Q);
% BFI3 = sum(BQ3)/sum(Q);

if nargin >= 3
    figure
    subplot(3,1,1)
    hold on
    plot(Date,Q,Date,BQ1,Date,BQ2,Date,BQ3,Date,SQ2)
    xlabel('Date')
    ylabel('Discharge (l/s)')
    legend('Discharge','Baseflow Linear','Baseflow 2par','Baseflow 3 par',...
        'Stormflow 2 par','Location','NorthWest')
    box on

    subplot(3,1,2)
    hold on
    plot(Date,LogQ,Date,LogBQ1,Date,LogBQ2,Date,LogBQ3)
    xlabel('Date')
    ylabel('Log(Discharge) log(l/s)')
    legend('Discharge','Baseflow Linear','Baseflow 2par','Baseflow 3 par',...
        'Location','NorthWest')
    box on

    subplot(3,1,3)
    plot(Date(1:n),R,Date(1:n),M,Date,LogQ,...
        DateTau,LogQ(and(R>=lim,M<0)),'o',...
        DateTau,RTau,DateTau,MTau);
    legend('Coeff. R^2','Regression slope','Log Discharge',...
        'Identified linear','Behavioural R^2','Behavioural Slope',...
        'Location','NorthWest')
    box on
end

%% AUXILIARY FUNCTIONS

% function [BQ,SQ] = par1(Q,k)
% % Initialise variables.
% BQ = zeros(size(Q));
% SQ = zeros(size(Q));
% BQ(1) = Q(1);
% for i = 2:length(Q)
%     BQ(i) = min(k/(2-k)*BQ(i-1)+(1-k)/(2-k)*Q(i),Q(i));
%     SQ(i) = Q(i) - BQ(i);
% end
% 
% function [BQ,SQ] = par2(Q,k,C)
% % Initialise variables.
% BQ = zeros(size(Q));
% SQ = zeros(size(Q));
% BQ(1) = Q(1);
% for i = 2:length(Q)
%     BQ(i) = min(k/(1+C)*BQ(i-1)+(C)/(1+C)*Q(i),Q(i));
%     SQ(i) = Q(i) - BQ(i);
% end

function [BQ,SQ] = par3(Q,k,C,alpha)
% Initialise variables.
BQ = zeros(size(Q));
SQ = zeros(size(Q));
BQ(1) = Q(1);
% h = waitbar(0,'Calculating baseflow...');
for i = 2:length(Q)
    BQ(i) = min(k/(1+C)*BQ(i-1)+(C)/(1+C)*(Q(i)+alpha*Q(i-1)),Q(i));
    SQ(i) = Q(i) - BQ(i);
    % waitbar(i/length(Q))
end
% close(h);