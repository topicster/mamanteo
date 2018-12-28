function [Q_Month,Q_Year,Q_Avg_Month,Q_Avg_Year,Q_Matrix,Q_YMin,Q_YMax] = iMHEA_MonthlyFlow(Date,Q,varargin)
%iMHEA Calculation of monthly and annual Discharge averages.
% [Q_Month,Q_Year,Q_Avg_Month,Q_Avg_Year,Q_Matrix] =
% iMHEA_MonthlyRain(Date,Q,flag).
%
% Input:
% Date = dd/mm/yyyy hh:mm:ss [date format].
% Q    = Discharge [l/s or l/s/km2].
% flag = leave empty NOT to graph plots.
%
% Output:
% Q_Month     = Time series of monthly discharge [l/s or l/s/km2].
% Q_Year      = Time series of annual discharge [year and l/s or l/s/km2].
% Q_Avg_Month = 12 average monthly discharge values [l/s or l/s/km2].
% Q_Avg_Year  = Annual discharge value [l/s or l/s/km2].
% Q_Matrix    = Matrix of discharge data (Year vs Months) [l/s or l/s/km2].
% Q_Min_Year  = Time series of minimum annual precipitation [year and mm].
% Q_Max_Year  = Time series of maximum annual precipitation [year and mm].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in September, 2017
% Last edited in November, 2018

%% PROCESS
Years = year(Date);
n = max(Years)-min(Years)+1; % Number of years
Months = month(Date);

Q_Year = zeros(n,1);
Q_YMin = zeros(n,2);
Q_YMax = zeros(n,2);
matrixQM1 = zeros(12,n);
sizeQM1 = zeros(12,n);

for i = 1:n
    % Annual mean
    Q_Year(i) = nanmean(Q(Years==min(Years)+i-1));
    % Position of the annual minimum
    [Q_YMin(i,1),MinPos] = nanmin(Q(Years==min(Years)+i-1));
    Q_YMin(i,2) = day(Date(MinPos),'dayofyear');
    % Position of the annual maximum
    [Q_YMax(i,1),MaxPos] = nanmax(Q(Years==min(Years)+i-1));
    Q_YMax(i,2) = day(Date(MaxPos),'dayofyear');
    for j = 1:12
        matrixQM1(j,i) = nanmean(Q(and(Years==min(Years)+i-1,Months==j)));
        sizeQM1(j,i) = length(Q(and(Years==min(Years)+i-1,Months==j)));
    end
end

%% GENERATE OUTPUT VARIABLES
Q_Avg_Month = nansum(matrixQM1.*sizeQM1,2)./nansum(sizeQM1,2);
Q_Avg_Year = nanmean(Q_Year);

% Discard extremes without data
Qdm = month(Date);
matrixQM1(1:Qdm(1)-1) = NaN;
matrixQM1(end-(12-Qdm(end))+1:end) = NaN;
% Reorganise outputs
Q_Month = matrixQM1(:);
Q_Matrix = matrixQM1';

Q_Year = [(min(Years):max(Years))' , Q_Year];
Q_YMin = [(min(Years):max(Years))' , Q_YMin];
Q_YMax = [(min(Years):max(Years))' , Q_YMax];

%% PLOT RESULTS
if nargin >= 3
    figure
    subplot(3,1,1)
    bar((1:length(Q_Month))',Q_Month,'DisplayName',inputname(2));
    hold on
    grid on
    box on
    title('Monthly data')
    legend('show')
    ylabel('Discharge [l/s]')
    set(gca,'XTick',(1:length(Q_Month)),...
        'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});

    subplot(3,1,2)
    bar((1:12)',Q_Avg_Month,'DisplayName',inputname(2));
    hold on
    Q_MatrixPlot = cat(1,nan(1,12),Q_Matrix,nan(1,12));
    boxplot(Q_MatrixPlot,'PlotStyle','compact')
    grid on
    box on
    title('Average monthly data')
    legend('show')
    ylabel('Discharge [l/s]')
    set(gca,'Xlim',[0 13],...
        'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},...
        'XTick',(1:12));

    subplot(3,1,3)
    bar(Q_Year(:,1),Q_Year(:,2),'DisplayName',inputname(2));
    hold on
    box on
    title('Annual data')
    legend('show')
    ylabel('Discharge [l/s]')

    drawnow
end