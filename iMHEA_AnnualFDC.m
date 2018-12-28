function [FDCs] = iMHEA_AnnualFDC(Date,Q,HYDates)
%iMHEA Annual Flow Duration Curves.
% [FDCs,FDCunc] = iMHEA_AnnualFDC(Q,HYDates) calculates annual FDCs
% aggregated for the dates indicated in HYDates.
%
% Input:
% Date    = dd/mm/yyyy hh:mm:ss [date format]
% Q       = Discharge [l/s, l/s/km2, m3/s, mm, etc.]
% HYDates = Hydrological years for annual FDC [matrix] [datetime]
%            each row must be [initial date, final date]
%
% Output:
% FDCs    = Annual flow duration curves summarised in 100 values [Q vs %]
%           between dates defined in HYDates.
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in November, 2018
% Last edited in November, 2018

%% PROCESS

nY = size(HYDates,1);
FDCs = nan(100,nY+1);

for i = 1:nY
    [~,~,~,~,FDC100] = iMHEA_FDC(Q(and(Date>=HYDates(i,1),Date<=HYDates(i,2))));
    FDCs(:,i+1) = FDC100(:,2);
end
FDCs(:,1) = FDC100(:,1);
