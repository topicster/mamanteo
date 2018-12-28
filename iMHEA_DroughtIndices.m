function [DIndices] = iMHEA_DroughtIndices(Date,P,Thr,flag1,flag2,varargin)
%iMHEA Calculates drought indices using the threshold specified as input
%      (Beyene et al, 2014; Van Loon, 2013).
% [DIndices] = iMHEA_DroughtIndices(Date,P,Thr,flag)
%
% Input:
% Date    = dd/mm/yyyy hh:mm:ss [date format]
% P       = evaluated variable: precipitation, streamflow [mm, l/s, m3/s]
% Thr      = 1-year threshold (366 values) [mm]
% flag1   = 1: smooth series with 30-day moving average; 0: use original
% flag2   = 1: graph plots; 0 or empty: NOT to graph plots
%
% Output:
% DIndices = Matriz of drought indices [10 elements].
%            1: Number of years
%            2: Mean annual precipitation
%            3: Standard deviation for annual precipitation
%            4: Number of droughts per year
%            5: Mean drought duration
%            6: Standard deviation of drought duration
%            7: Maximum drought duration
%            8: Mean drought deficit
%            9: Standard deviation of drought deficit
%           10: Maximum drought deficit
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in September, 2017
% Last edited in November, 2018

%% CALCULATION OF DROUGHT INDICES

% Check input
if nargin < 4
    flag1 = 1;
    flag2 = 0;
elseif nargin < 5
    flag2 = 0;
end

% Initiliase variables
n = 366;
int = 30;
int2 = ceil(int/2);

%% Smooth interdaily data and consider irregular timing

if flag1 == 1
    % 30-day moving average of precipitation
    nP = length(P);
    P30D = nan(size(P));
    % Extend the tails of the time series for moving average
    P = [nan(int2-1,1);P;nan(int2,1)];
    for i = 1:nP
        P30D(i) = nanmean(P(i:i+int-1));
    end
    % Remove the tails
    P(1:int2-1) = [];
    P(end-int2+1:end) = [];
else
    P30D = P;
end

%% Replicating the threshold for the assessed data
PThr = nan(size(P30D));
Pdays = day(Date,'dayofyear');
Pyears = year(Date);
P_ly = leapyear(Pyears);

for i = 1:60
    PThr(Pdays==i) = Thr(i);
end
Pdays(~P_ly) = Pdays(~P_ly)+1;
for i = 61:n
    PThr(Pdays==i) = Thr(i);
end

%% Drought occurrence when variable is below the threshold
Drought = P30D<PThr;

% Assigning individual numbers to each non-drought event
nod = zeros(size(Drought));
counter = 0;
if ~Drought(1); nod(1)=1; counter=1; else; nod(1)=0; end
for i = 2:length(Drought)
    if Drought(i)==Drought(i-1) && ~Drought(i)
        nod(i)=nod(i-1);
    elseif ~Drought(i)
        nod(i)=counter+1;
        counter=counter+1;
    end
end

%% Pooling droughts if inter-event period is 10 days or less
mind = 10;
% Calculate inter-event durations
nod2 = nod;
nodur = zeros(max(nod),1);
for i = 1:max(nod)
    nodur(i) = sum(nod(nod==i))/i;
    % Pooling droughts when inter-event period is less than 10 days
    if nodur(i)<mind
        % Removing inter-event period
        nod2(nod2==i)=0;
        nod2(nod2>i)=nod2(nod2>i)-1;
        % Assigning drought events
        Drought(nod==i)=1;
    end
end
% Alternative pooling for 2-day inter-event period
% for i = 3:length(Drought)
%     % Pooling droughts if inter-event period is 2 days or less
%     % if Drought(i)==0 && (Drought(i-1)==1 && (Drought(i+1)==1 || Drought(i+2)==1))
%     if 
%         Drought(i)=true;
%     end
% end

% Assigning individual numbers to each drought event
d = zeros(size(Drought));
counter = 0;
if Drought(1); d(1)=1; counter=1; else; d(1)=0; end
for i = 2:length(Drought)
    if Drought(i)==Drought(i-1) && Drought(i)
        d(i)=d(i-1);
    elseif Drought(i)
        d(i)=counter+1;
        counter=counter+1;
    end
end

%% Discard droughts of less than 10 days of duration
% Calculate duration of drought events
d2 = d;
dur = zeros(max(d),1);
for i = 1:max(d)
    dur(i) = sum(d(d==i))/i;
    % Discarding droughts of less than 10 days of duration
    if dur(i)<mind
        % Removing drought events
        d2(d2==i)=0;
        d2(d2>i)=d2(d2>i)-1;
        % Discarding droughts events
        Drought(d==i)=0;
    end
end
dur2 = dur;
dur2(dur<10)=[];

%% Calculate deficit volume of drought events
def = zeros(max(d),1);
DVolume = P30D-PThr;
for i = 1:max(d)
    def(i) = sum(DVolume(d==i&DVolume<0));
end
def2 = def(dur>=10);

%% Drought indices

% Annual rainfall
[~,PYear] = iMHEA_MonthlyRain(Date,P);

% Number of years
DIndices(1) = length(Date)/365.25;
% DIndices(1) = PYear(end,1)-PYear(1,1)+1;
% n = max(year(Date))-min(year(Date))+1;

% Average annual precipitation
DIndices(2) = nanmean(PYear(:,2));

% Standard deviation of annual precipitation
DIndices(3) = nanstd(PYear(:,2));

% Number of droughts per year
DIndices(4) = length(dur2)/DIndices(1);

% Average duration of droughts
DIndices(5) = nanmean(dur2);

% Standard deviation of duration of droughts
DIndices(6) = nanstd(dur2);

% Maximum duration of droughts
DIndices(7) = nanmax(dur2);

% Average deficit volume during droughts
DIndices(8) = nanmean(def2);

% Standard deviation of deficit volume during droughts
DIndices(9) = nanstd(def2);

% Maximum deficit volume during droughts
DIndices(10) = nanmin(def2);

%% Print results

fprintf('\n')
fprintf('DROUGHT ANALYSIS PERIOD: %4i to %4i\n',PYear(1,1),PYear(end,1))
fprintf('\n')
fprintf('Drought indices:\n')
fprintf('Number of years: %5.2f\n',DIndices(1))
fprintf('Annual mean: %5.2f (+/-) %5.2f [mm]\n',DIndices(2),DIndices(3))
fprintf('Number of droughts: %5.2f [per year]\n',DIndices(4))
fprintf('Mean duration: %5.2f (+/-) %5.2f [day]\n',DIndices(5),DIndices(6))
fprintf('Max duration: %5.2f [day]\n',DIndices(7))
fprintf('Mean deficit: %5.2f (+/-) %5.2f [mm]\n',DIndices(8),DIndices(9))
fprintf('Max deficit: %5.2f [mm]\n\n',DIndices(10))

%% Plot results

if flag2 == 1
    
    figure
    
    % Plot time series
    subplot(3,3,1:3)
    hold on
    plot(Date,P30D,'LineWidth',1,'DisplayName','Assessed variable time series')
    plot(Date,PThr,'LineWidth',1.5,'DisplayName','Drought thrershold')
    % Plot settings
    box on
    grid on
    % Labels
    title('DROUGHT ANALYSIS')
    ylabel('variable [mm]')
    % Legend
    legend('location','Northwest')
    
    % Plot drought events
    subplot(3,3,4:6)
    hold on
    plot(Date,nod,'LineWidth',2,'DisplayName','Original inter-event periods')
    plot(Date,nod2,'LineWidth',2,'DisplayName','New inter-event periods')
    plot(Date,d,'LineWidth',2.5,'DisplayName','Original drought events')
    plot(Date,d2,'LineWidth',2.5,'DisplayName','New drought events')
    % Plot settings
    set(gca,'yscale','log')
    box on
    grid on
    % Labels
    xlabel('date')
    ylabel('event id')
    % Legend
    legend('location','Northwest')
    
    % Plot drought duration and deficit
    subplot(3,3,7)
    hold on
    plot(dur,'LineWidth',1,'DisplayName','Drought duration total')
    plot(def,'LineWidth',1,'DisplayName','Drought deficit total')
    % Plot settings
    box on
    grid on
    % Labels
    xlabel('drought number')
    ylabel('(-)deficit [mm], duration [day]')
    title('Original drought events')
    % Legend
    legend('location','Northeast')
    
    subplot(3,3,9)
    hold on
    plot(dur2,'LineWidth',1,'DisplayName','Drought duration')
    plot(def2,'LineWidth',1,'DisplayName','Drought deficit')
    % Plot settings
    box on
    grid on
    % Labels
    xlabel('drought number')
    ylabel('(-)deficit [mm], duration [day]')
    title('Pooled and depured drought events')
    % Legend
    legend('location','Northeast')
    
    % Plot threshold
    subplot(3,3,8)
    hold on
    for i = Pyears(1):Pyears(end)
        plot(P30D(Pyears==i),'LineWidth',1,'color',0.75-((Pyears(end)-i+1)/(Pyears(end)-Pyears(1)+1))*0.5*[1 1 1])
    end
    a = plot(Thr,'r','LineWidth',5,'DisplayName','Threshold');
    % Plot settings
    box on
    grid on
    set(gca,'XLim',[0 367])
    % Labels
    xlabel('[day]')
    ylabel('variable [mm]')
    title('Annual data and threshold')
    % Legend
    legend(a)
end
