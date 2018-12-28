function [FDC,R2FDC,IRH,Ptile,FDC100] = iMHEA_FDC(Q,varargin)
%iMHEA Flow Duration Curve.
% [FDC,R2FDC,IRH,Ptile,FDC100] = iMHEA_FDC(Q) calculates FDC, slope of the
% FDC, IRH, percentiles, and summarised FDC.
%
% Input:
% Q      = Discharge [l/s, l/s/km2, m3/s, mm, etc.].
%
% Output:
% FDC    = Flow duration curve complete values [Q vs %].
% R2FDC  = Slope of the FDC between 33%-66% / Mean flow [-].
% IRH    = Hydrological Regulation Index [-].
% Ptile  = Percentiles from the FDC, including:
%          Q95   = 05th Percentile [l/s].
%          Q75   = 25th Percentile [l/s].
%          Q66   = 33rd Percentile [l/s].
%          Q50   = 50th percentile [l/s].
%          Q33   = 66th Percentile [l/s].
%          Q25   = 75th Percentile [l/s].
%          Q10   = 90th Percentile [l/s].
% FDC100 = Flow duration curve summarised in 100 values [Q vs %].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in May, 2014
% Last edited in November, 2018

%% PROCESS
Q = Q(~isnan(Q)); % Only when data exists
k = length(Q); % Number of elements in Q

% Flow Duration Curve.
FDC(:,1) = 100*(1-((1:k)-.44)/(k+.12)); % Gringorten (1963)
FDC(:,2) = sort(Q); % Vector of sorted discharge

% Percentiles from the FDC.
% Interpolation needed because frequency does not match exactly with the
% required percentile.
Ptile(1) = prctile(FDC(:,2),5);
Ptile(2) = prctile(FDC(:,2),25);
Ptile(3) = prctile(FDC(:,2),33);
Ptile(4) = prctile(FDC(:,2),50);
Ptile(5) = prctile(FDC(:,2),66);
Ptile(6) = prctile(FDC(:,2),75);
Ptile(7) = prctile(FDC(:,2),90);

% Slope of the FDC between 33 and 66 %.
R2FDC = (log10(Ptile(3)) - log10(Ptile(5)))/(.66-.33);
if ~isreal(R2FDC)
    R2FDC = -inf;
end

% Hydrological Regulation Index.
auxFDC = FDC(:,2);
auxFDC(FDC(:,1)<50) = Ptile(4);
IRH = sum(auxFDC) / sum(FDC(:,2));

% FDC divided in 100 values
FDC100 = nan(100,2);
FDC100(:,1) = 100*(1-((1:100)-.44)/(100+.12));
for i = 1:100
    FDC100(i,2) =  prctile(FDC(:,2),i);
end

%% PLOT RESULTS
if nargin >= 2
    % Plot in logarithmic scale.
    figure
    subplot(1,2,1)
    semilogy(FDC(:,1),FDC(:,2))
    hold on
    semilogy([66 33],[Ptile(3) Ptile(5)],'r','LineWidth',1.5)
    xlabel('Exceedance probability')
    ylabel('Discharge')
    title('Flow Duration Curve')
    legend('FDC','Slope 33-66%')
    
    subplot(1,2,2)
    plot(FDC(:,1),FDC(:,2))
    hold on
    area(FDC(:,1),auxFDC)
    xlabel('Exceedance probability')
    ylabel('Discharge')
    title('Flow Duration Curve')
    legend('FDC','Area under 50%')
end