function [P_Month,P_Year,P_Avg_Month,P_Avg_Year,P_Matrix,P_YMin,P_YMax] = iMHEA_MonthlyRain(Date,P,varargin)
%iMHEA Calculation of monthly and annual Precipitation averages.
% [P_Month,P_Year,P_Avg_Month,P_Avg_Year,P_Matrix] =
% iMHEA_MonthlyRain(Date,P,flag).
%
% Input:
% Date = dd/mm/yyyy hh:mm:ss [date format].
% P    = Precipitation [mm].
% flag = leave empty NOT to graph plots.
%
% Output:
% P_Month     = Time series of monthly precipitation [mm].
% P_Year      = Time series of annual precipitation [year and mm].
% P_Avg_Month = 12 average monthly precipitation values [mm].
% P_Avg_Year  = Annual precipitation value [mm].
% P_Matrix    = Matrix of precipitation data (Year vs Months) [mm].
% P_Min_Year  = Time series of minimum annual discharge [year and l/s or l/s/km2].
% P_Max_Year  = Time series of maximum annual discharge [year and l/s or l/s/km2].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in September, 2017
% Last edited in November, 2018

%% PROCESS
Years = year(Date);
n = max(Years)-min(Years)+1; % Number of years
Months = month(Date);

P_Year = zeros(n,1);
P_YMax = zeros(n,2);
P_YMin = zeros(n,2);
matrixPM1 = zeros(12,n);
sizePM1 = zeros(12,n);

for i = 1:n
    % Annual accumulation
    P_Year(i) = nansum(P(Years==min(Years)+i-1));
    % Position of the annual minimum
    [P_YMin(i,1),MinPos] = nanmin(P(Years==min(Years)+i-1));
    P_YMin(i,2) = day(Date(MinPos),'dayofyear');
    % Position of the annual maximum
    [P_YMax(i,1),MaxPos] = nanmax(P(Years==min(Years)+i-1));
    P_YMax(i,2) = day(Date(MaxPos),'dayofyear');
    for j = 1:12
        matrixPM1(j,i) = nansum(P(and(Years==min(Years)+i-1,Months==j)));
        sizePM1(j,i) = length(P(and(Years==min(Years)+i-1,Months==j)));
    end
end

%% GENERATE OUTPUT VARIABLES
P_Avg_Month = nansum(matrixPM1.*sizePM1,2)./nansum(sizePM1,2);
P_Avg_Year = nanmean(P_Year);

% Discard extremes without data
Pdm = month(Date);
matrixPM1(1:Pdm(1)-1) = NaN;
matrixPM1(end-(12-Pdm(end))+1:end) = NaN;
% Reorganise outputs
P_Month = matrixPM1(:);
P_Matrix = matrixPM1';

P_Year = [(min(Years):max(Years))' , P_Year];
P_YMin = [(min(Years):max(Years))' , P_YMin];
P_YMax = [(min(Years):max(Years))' , P_YMax];

%% PLOT RESULTS
if nargin >= 3
    figure
    subplot(3,1,1)
    bar((1:length(P_Month))',P_Month,'DisplayName',inputname(2));
    hold on
    grid on
    box on
    title('Monthly data')
    legend('show')
    ylabel('Precipitation [mm]')
    set(gca,'XTick',(1:length(P_Month)),...
        'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});

    subplot(3,1,2)
    bar((1:12)',P_Avg_Month,'DisplayName',inputname(2));
    hold on
    P_MatrixPlot = cat(1,nan(1,12),P_Matrix,nan(1,12));
    boxplot(P_MatrixPlot,'PlotStyle','compact')
    grid on
    box on
    title('Average monthly data')
    legend('show')
    ylabel('Precipitation [mm]')
    set(gca,'Xlim',[0 13],...
        'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},...
        'XTick',(1:12));

    subplot(3,1,3)
    bar(P_Year(:,1),P_Year(:,2),'DisplayName',inputname(2));
    hold on
    box on
    title('Annual data')
    legend('show')
    ylabel('Precipitation [mm]')

    drawnow
end