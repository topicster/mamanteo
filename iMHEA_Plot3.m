function iMHEA_Plot3(Date,varargin)
%iMHEA Plot of rainfall-runoff data.
% iMHEA_Plot3(Date,Variable,...)
%
% Input:
% Date = dd/mm/yyyy hh:mm:ss (number format).
% Variable: Precipitation [mm] or Discharge [l/s, m3/s, mm].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in May, 2014
% Last edited in November, 2017

%% PLOT

if nargin < 2
    disp('Not enough arguments.')
    return
end

n = nargin-1;
figure
XLIM = [datetime([inf inf inf]) datetime([0 0 0])];

for i = 1:n
    subplot(n,1,i)
    plot(Date,varargin{i})
    Xlimaux = get(gca,'XLim');
    XLIM = [min(XLIM(1),Xlimaux(1)),max(XLIM(2),Xlimaux(2))];
end

for i = 1:n
    subplot(n,1,i)
    set(gca,'XLim',XLIM)
    grid on
    box on
end