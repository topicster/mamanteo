function [NewDate,NewP,CumP,VoidP,MaxP] = iMHEA_Aggregation(Date,P,scale,varargin)
%iMHEA Agregation of rainfall (agregation within an interval).
% [NewDate,NewP,CumP,VoidP,MaxP] = iMHEA_Aggregation(Date,P,scale,flag)
% aggregates precipitation data using tipping bucket counting.
%
% Input:
% Date  = dd/mm/yyyy hh:mm:ss [date format].
% P     = Precipitation [mm].
% scale = Agregation interval [min].
% flag  = leave empty NOT to run the data voids assessment and plots.
%
% Output:
% NewDate   = dd/mm/yyyy hh:mm:ss [date format] at specified interval.
% NewP      = Agregated Precipitation [mm].
% CumP      = Cumulative rainfall [mm].
% VoidP     = Void intervals [mm].
% MaxP      = Maximum intensity for specified interval [mm].
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
    [Voids] = iMHEA_Voids(Date,P,1);
else
    % Run data gap assessment.
    [Voids] = iMHEA_Voids(Date,P);
end
% Allocate 0 mm to data gaps temporarily.
P(isnan(P)) = 0;
% Transform date variables for easier processing.
Date = datenum(Date);

%% AGGREGATION
% Numeric value of 1 minute: 1/1440
nd = 1440/scale; % Number of intervals per day
DI = ceil(min(Date)*nd); % Initial date
DF = ceil(max(Date)*nd); % Final date
NewDate = (DI:DF)'; % Equally spaced time interval
n = length(NewDate); % Number of intervals
NewP = zeros(size(NewDate)); % Initialise agregation
% Delete zero events to help process relevant data only.
Date(P==0) = [];
P(P==0) = [];
k = length(P); % Length of input data

if nd*(Date(1)) == NewDate(1)
    j = 2; % Data counter
    NewP(1) = P(1);
else
    j = 1; % Data counter
end
for i = j:n
    % Aggregate values
    while j<=k && nd*Date(j)<=NewDate(i) % && nd*Date(j)>NewDate(i-1) 
        NewP(i) = NewP(i) + P(j);
        j = j+1;
    end
end
% Fill gaps between data when there is only one value missing.
for i = 2:(n-1)
    if isnan(NewP(i))
        NewP(i) = 0;
    end
end

%% PREPARE THE DATA VECTORS AT THE SPECIFIED SCALE
CumP = cumsum(NewP); % Initialise accumulation
NewDate = datetime(NewDate/nd,'ConvertFrom','datenum'); % Rescale the date
Date = datetime(Date,'ConvertFrom','datenum'); % Restore the original date
% Placing data gaps again in the aggregated vectors.
VoidP = NewP;
for i = 1:size(Voids,1)
    CumP(NewDate>Voids(i,1) & NewDate<Voids(i,2)) = NaN;
    NewP(NewDate>Voids(i,1) & NewDate<Voids(i,2)) = NaN;
end
VoidP(~isnan(NewP)) = NaN;
% Correct the last row.
if NewP(end) == 0 && isnan(NewP(end-1))
    VoidP(end) = NewP(end);
    NewP(end) = NaN;
    CumP(end) = NaN;
end
MaxP = nanmax(NewP); % Maximum intensity

%% PLOT THE RESULTS
if nargin > 3
    fprintf('Plot: Data1 (Original), Data2 (Agregated), Data3 (Voids).\n')
    fprintf('\n')
    iMHEA_Plot2(Date,P,NewDate,NewP,NewDate,VoidP);  % Plot data
    iMHEA_Plot3(NewDate,NewP,VoidP,CumP);  % Plot data
end

%% ASSIGN SINGLE OUTPUT IF NEEDED
if nargout == 1
    NewDate = [datenum(NewDate),NewP,CumP,VoidP];
end