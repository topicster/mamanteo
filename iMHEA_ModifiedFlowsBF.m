function [Qmod] = iMHEA_ModifiedFlowsBF(Qrem,Qdiv,Qbf,tts,frr)
%iMHEA Calculates the modified hydrograph after flow diversion scaled by
%      the base flow.
% [Qmod] = iMHEA_ModifiedFlowsBF(Qrem,Qdiv,Qbf,factors,frr)
%
% Input:
% Qrem  = Time series of remaining flows in the stream [l/s or l/s/km2]
% Qdiv  = Time series of diverted flows [l/s or l/s/km2]
% Qbf   = Time series of base flows to scale diversion [l/s or l/s/km2]
% tts   = Resicence times used to modify the hydrograph at Q time scale
% frr   = Recovery rate (default 50%) [scalar]
%
% Output:
% Qmod  = Time series of diverted flows [l/s or l/s/km2]
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in November, 2018
% Last edited in December, 2018

%% Calculate modified flows

% Counters
n1 = length(Qrem);
n2 = length(tts);
tts = tts/nansum(tts);

% Start from the remaining flows in the stream
Qmod = Qrem;

% Modified flows
h = waitbar(0,'Modifying hydrograph...');
for i = 1:n1
    % Use only those flows diverted to the mamanteo
    if Qdiv(i)>0
        fbf = nan(size(tts));
        fbf(1:min(n1-i+1,n2)) = Qbf(i:min(i+n2-1,n1));
        fbf = fbf/nansum(fbf);
        factors = (tts .* fbf)/nansum(tts .* fbf);
        % Calculate delayed flows using the residence times and recovery rate
        delayed = Qdiv(i) * frr .* factors;
        % Add the delayed flows to the remaining flow time series
        Qmod(i:min(i+n2-1,n1)) = Qmod(i:min(i+n2-1,n1))+delayed(1:min(n1-i+1,n2));
        waitbar(i/n1)
    end
end
close(h);
