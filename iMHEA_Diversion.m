function [Qdiv,Qdiv_1yr] = iMHEA_Diversion(Date,Q,QMin,QMax,WSDates,varargin)
%iMHEA Calculates the diversion flows from Q using limits QMin-QMax.
% [Qdiv,Qdiv_1yr] = iMHEA_Diversion(Date,Q,QMin,QMax,WSDates)
%
% Input:
% Date     = dd/mm/yyyy hh:mm:ss [date format]
% Q        = Time series of original flows [l/s or l/s/km2] for Date
% QMin     = Minimum flows left in the stream [vector] [l/s/km2]
% QMax     = Maximum diversion intake capacity [vector] [l/s/km2]
% WSDates  = Wet seasons for diversion [matrix] [datetime]
%            each row must be [initial date, final date]
%
% Output:
% Qdiv     = Time series of diverted flows [l/s or l/s/km2] for Date
%            for the last combination of QMin & Qmax only
% Qdiv_1yr = Diverted flows for each combination [matrix] for each ws
%            column 1-2:   wet season years
%            column 3-end: combination QMax & QMin
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in November, 2018
% Last edited in November, 2018

%% Initialise variables

% QMin limits
if nargin < 3
    QMin = 0;
end

% QMax limits
if nargin < 4
    QMax = inf;
end

% Wet season limits
if nargin < 5
    WSDlims = datenum([Date(1) Date(end)]);
    WSDates = datetime(WSDlims,'ConvertFrom',datenum);
end

% Initiliase the matrix
WSDlims = datenum(WSDates);
nW = size(WSDates,1);
nM = length(QMax);
nm = length(QMin);
Qdiv_1yr = zeros(nW,2+nM*nm);

% Assign years in the wet seasons as indices
for i = 1:nW
    Qdiv_1yr(i,1:2) = [year(WSDates(i,1)) year(WSDates(i,2))];
end

% Convert Date to numbers if needed
Date = datenum(Date);
int = median(diff(Date));

%% Diverted flows
disp('Estimating diverted flows')

for i = 1:length(QMax)
    for j = 1:length(QMin)
        
        % Flow limits
        Qmax = QMax(i); % in l/s
        Qmin = QMin(j); % in l/s/km2
        
        % Leave a minimum flow in the stream
        Qdiv = Q-Qmin;
        % Correct data if the minimum flow is greater than the available
        Qdiv(Qdiv<0) = 0;
        % Limit the capacity of the system to the maximum
        Qdiv(Qdiv>Qmax) = Qmax;
        
        % Only implement mamanteo during the defined wet seasons
        % Create axuliar column
        Qdiv(:,2) = zeros(size(Qdiv(:,1)));
        for k = 1:nW
            % Fill with data only between wet season limits
            Qdiv(and(Date>WSDlims(k,1),Date<WSDlims(k,2)),2) = Qdiv(and(Date>WSDlims(k,1),Date<WSDlims(k,2)),1);
        end
        % Discard original values
        Qdiv(:,1) = [];
        
        % Calculate potential diverted flows per wet season
        for k = 1:nW
            Qdiv_1yr(k,2+i+(j-1)*nM) = nansum(Qdiv(Date>=WSDlims(k,1) & Date<=WSDlims(k,2)))*int*86400/1000000;
        end
        
    end
end

% Flows remaining in the stream for the last combination QMin & QMax
disp('Estimating flows remaining in the stream after diversion')
Qdiv(:,2) = Q - Qdiv(:,1);
