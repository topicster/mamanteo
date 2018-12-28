function [Qmod] = iMHEA_ModifiedFlows(Qrem,Qdiv,lags,frr)
%iMHEA Calculates the modified hydrograph after flow diversion.
% [Qmod] = iMHEA_ModifiedFlows(Qrem,Qdiv,factors,frr)
%
% Input:
% Qrem  = Time series of remaining flows in the stream [l/s or l/s/km2]
% Qdiv  = Time series of diverted flows [l/s or l/s/km2]
% lags  = Residence times used to modify the hydrograph at Q time scale
% frr   = Recovery rate (default 50%) [scalar]
%
% Output:
% Qmod  = Time series of diverted flows [l/s or l/s/km2]
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in November, 2018
% Last edited in December, 2018

%% Calculate delayed flows

% Counters
n1 = length(Qrem);
n2 = length(lags);

% Start from the remaining flows in the stream
Qmod = Qrem;

% Modified flows
h = waitbar(0,'Modifying hydrograph...');
for i = 1:n1
    % Use only those flows diverted to the mamanteo
    if Qdiv(i)>0
        % Calculate delayed flows using the residence times and recovery rate
        delayed = Qdiv(i) * frr * lags;
        % Add the delayed flows to the remaining flow time series
        Qmod(i:min(i+n2-1,n1)) = Qmod(i:min(i+n2-1,n1))+delayed(1:min(n1-i+1,n2));
        waitbar(i/n1)
    end
end
close(h);
