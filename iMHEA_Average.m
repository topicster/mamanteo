function [NewDate,NewQ,CumQ,VoidQ,MeanQ,MaxQ,MinQ] = iMHEA_Average(Date,Q,scale,varargin)
%iMHEA Agregation of hydrological data (average within an interval).
% [NewDate,NewQ,CumQ,VoidQ,MeanQ,MaxQ,MinQ] =
% iMHEA_Average(Date,Q,scale,flag) averages discharge data.
%
% Input:
% Date  = dd/mm/yyyy hh:mm:ss [date format].
% Q     = Stage or Discharge [l/s, m3/s, mm].
% scale = Agregation interval [min].
% flag  = leave empty NOT to run the data voids assessment and plots.
%
% Output:
% NewDate   = dd/mm/yyyy hh:mm:ss [date format] at specified interval.
% NewQ      = Average stage or discharge [l/s, m3/s, mm].
% CumQ      = Cumulative discharge [l/s, m3/s, mm].
% VoidP     = Void intervals [l/s, m3/s, mm].
% MeanQ     = Mean value for specified interval [l/s, m3/s, mm].
% MaxQ      = Maximum value for specified interval [l/s, m3/s, mm].
% MinQ      = Minimum value for specified interval [l/s, m3/s, mm].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in May, 2014
% Last edited in November, 2017

%% INITIALISE VARIABLES AND IDENTIFY DATA GAPS
% Move date by 0.25 seconds to avoid numerical errors.
Date = Date - 0.25/86400;
if nargin > 3
    % Run data gap assessment and print inventory.
    [Voids] = iMHEA_Voids(Date,Q,1);
else
    % Run data gap assessment.
    [Voids] = iMHEA_Voids(Date,Q);
end
% Transform date variables for easier processing.
Date = datenum(Date);

%% AVERAGE
% Numeric value of 1 minute: 1/1440
nd = 1440/scale; % Number of intervals per day
DI = ceil(min(Date)*nd); % Initial date
DF = ceil(max(Date)*nd); % Final date
NewDate = (DI:DF)'; % Equally spaced time interval
n = length(NewDate); % Number of intervals
NewQ = zeros(size(NewDate)); % Initialise agregation
% Delete NaN events to help process relevant data only.
Date(isnan(Q)) = [];
Q(isnan(Q)) = [];
k = length(Q); % Length of input data
% Set initial counter
if nd*(Date(1)) == NewDate(1)
    j = 2; % Data counter
    NewQ(1) = Q(1);
else
    j = 1; % Data counter
end
for i = j:n
    l = 0; % Interval data counter
    % Aggregate values
    while j<=k && nd*Date(j)<=NewDate(i) % && nd*Date(j)>NewDate(i-1) 
        NewQ(i) = NewQ(i) + Q(j);
        j = j+1;
        l = l+1;
    end
    NewQ(i) = NewQ(i)/l;
end
% Fill gaps between data when there is only one value missing.
for i = 2:(n-1)
    if isnan(NewQ(i))
        NewQ(i) = mean([NewQ(i-1),NewQ(i+1)]);
    end
end
% Fill remaining gaps with zeros to calculate cumQ
NewQ(isnan(NewQ))= 0;

%% PREPARE THE DATA VECTORS AT THE SPECIFIED SCALE
CumQ = cumsum(NewQ); % Initialise accumulation
NewDate = datetime(NewDate/nd,'ConvertFrom','datenum'); % Rescale the date
Date = datetime(Date,'ConvertFrom','datenum'); % Restore the original date
% Placing data gaps again in the aggregated vectors.
VoidQ = NewQ;
for i = 1:size(Voids,1)
    CumQ(NewDate>Voids(i,1) & NewDate<Voids(i,2)) = NaN;
    NewQ(NewDate>Voids(i,1) & NewDate<Voids(i,2)) = NaN;
end
VoidQ(~isnan(NewQ)) = NaN;
% Check initial and final values of Q for data existence
if NewQ(1) == 0 && NewQ(2)~=0
    VoidQ(1) = NewQ(1);
    NewQ(1) = NaN;
    CumQ(1) = NaN;
end
if NewQ(end) == 0 && NewQ(end-1)~=0
    VoidQ(end) = NewQ(end);
    NewQ(end) = NaN;
    CumQ(end) = NaN;
end
MeanQ = nanmean(NewQ); % Mean averaged value
MaxQ = nanmax(NewQ); % Maximum averaged value
MinQ = nanmin(NewQ); % Minimum averaged value

%% PLOT THE RESULTS
if nargin > 3
    fprintf('Plot: Data1 (Original), Data2 (Averaged), Data3 (Voids).\n')
    fprintf('\n')
    iMHEA_Plot2(Date,Q,NewDate,NewQ,NewDate,VoidQ);  % Plot data
    iMHEA_Plot3(NewDate,NewQ,VoidQ,CumQ);  % Plot data
end

%% ASSIGN SINGLE OUTPUT IF NEEDED
if nargout == 1
    NewDate = [datenum(NewDate),NewQ,CumQ,VoidQ];
end