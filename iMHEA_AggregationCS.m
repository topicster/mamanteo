function [NewDate,NewP,CumP,Single] = iMHEA_AggregationCS(Event_Date,Event_mm,scale,bucket,mintip,halves,varargin)
%iMHEA Cubic spline interpolation rainfall aggregation.
% [NewDate,NewP,CumP,Single] =
% iMHEA_AggregationCS(Event_Date,Event_mm,scale,Vb,mintip,halves,flag)
% aggregates precipitation data using cublic spline interpolation.
%
% Input:
% Event_Date= dd/mm/yyyy hh:mm:ss [date format].
% Event_mm  = Precipitation [mm].
% scale     = Agregation interval [min].
% bucket    = Rain gauge bucket volume [mm].
% mintip    = Aggregate at 1-min before interpolation [default: true].
% halves    = Add zero rates at estimated event endpoints [default: true].
% flag1     = leave empty NOT to graph data inventory and event plots.
%
% Output:
% NewDate   = dd/mm/yyyy hh:mm:ss [date format] at specified interval.
% NewP      = Agregated precipitation [mm].
% CumP      = Cumulative rainfall [mm].
% Single    = Single tip rainfall rates [mm].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in February, 2016
% Last edited in November, 2017

%% USER-DEFINED VARIABLES
% Minimum intensity to separate events: 0.2 mm h^-1 [Padron et al, 2015].
Minint = 0.2/1; % 0.2 mm over 1 hour
% Maximum intensity to merge events: 12.7 cm h^-1 [Onset, 2013],
% or 300 mm h^-1 [Manz, personal communication].
Maxint = 127;
% Intensity to distribute single tips: 3 mm h^{-1} [Wang et al, 2008].
Meanint = 3;
% Threshold intensity above which data is kept:
% 0.10 mm h^{-1} [Wang et al, 2008] or 1/2 Minint.
Lowint = min(0.1/60,Minint/120); % [mm min^{-1}]
% Move date by 0.25 seconds to avoid numerical or exportation errors.
Event_Date = Event_Date - 0.25/86400;
% Event_Date = Event_Date + 0.25/86400;

% Add zero rates at estimated initial and final event times, using (Vb/2):
% 1 half rate and 1 half partial tip [Sadler and Busscher, 1989].
% halves = true; % [default: true]

%% INITIALISE VARIABLES
fprintf('\n')
fprintf('RAINFALL AGGREGATION USING CUBIC SPLINE INTERPOLATION.\n')
if nargin < 2 || isempty(Event_mm)
    % Bucket volume in the rain gauge assumed to be 0.2 if not given.
    Event_mm = 0.2*ones(size(Event_Date));
end
if nargin < 3 || isempty(scale)
    % Aggregation time interval assumed to be 1 min
    scale = 1;
end
if nargin < 4 || isempty(bucket)
    % Bucket volume in the rain gauge assumed to be 0.2 if not given.
    bucket = 0.2;
end
if nargin < 5 || isempty(mintip)
    % Do not aggregate at 1-min before interpolation.
    mintip = true;
end
if nargin < 6 || isempty(halves)
    % Add zero rates at estimated event endpoints.
    halves = true;
end
% IT IS RECOMMENDED TO DEPURE UNREALISTIC EVENTS BEFORE RUNNING THE CODE
% Unrealistic events are those which happen consecutively after 1 second.
% [Event_mm] = iMHEA_Depure(Event_Date,Event_mm);

%% IDENTIFY DATA GAPS
if nargin > 6
    % Run data gap assessment and plot inventory.
    [Voids] = iMHEA_Voids(Event_Date,Event_mm,1,1);
    % Plot tips
    plot(Event_Date(Event_mm~=0),Event_mm(Event_mm~=0),'o','DisplayName','Events')
    set(gca,'YTickLabel',[],'Ylim',[-2 max(Event_mm)+1])
    legend('off'); legend('show');
    legend('location','SouthWest'); legend('boxoff')
else
    % Run data gap assessment.
    [Voids] = iMHEA_Voids(Event_Date,Event_mm,1);
end

%% TRANFORM DATES TO NUMBERS FOR EASIER PROCESSING
% Modified variables to process.
Event_Date = datenum(Event_Date);
NewEvent_Date = datenum(Event_Date);
NewEvent_mm = Event_mm;
% Allocate 0 mm to data gaps temporarily.
NewEvent_mm(isnan(Event_mm)) = 0;
% Delete zero events to help process relevant data only.
NewEvent_Date(NewEvent_mm==0) = [];
NewEvent_mm(NewEvent_mm==0) = [];

%% PREPROCESS RAINFALL EVENTS
nd = 1440; % Number of minutes per day or numeric value of 1 minute: 1/1440
% Maximum tip interval to separate events.
MaxT = 60*(1/nd)*bucket/Minint;
% Minimum tip interval to merge events.
MinT = 60*(1/nd)*bucket/Maxint;
% Aggregate events to avoid large intensities
if mintip == true
    % Aggregate data at 1-min scale.
    [NewEvent_Date,NewEvent_mm] = AggregateEvents(NewEvent_Date,NewEvent_mm);
else
    % Merge rainfall tips occurring at extremely short periods.
    [NewEvent_Date,NewEvent_mm] = MergeEvents(NewEvent_Date,NewEvent_mm,MinT);
end
% Adding a supporting initial extreme to avoid crashing the code later.
NewEvent_Date = cat(1,Event_Date(1)-MaxT,NewEvent_Date);
NewEvent_mm = cat(1,0,NewEvent_mm);
% Redistribute rainfall tips occurring at relatively long periods.
[NewEvent_Date,NewEvent_mm] = DivideEvents(NewEvent_Date,NewEvent_mm,MaxT);
% Redistribute rainfall over relatively long periods slightly shorter.
% [NewEvent_Date,NewEvent_mm] = DivideEvents(NewEvent_Date,NewEvent_mm,MaxT/2);
% Remove initial extreme to avoid crashing the code later.
NewEvent_Date(1) = [];
NewEvent_mm(1) = [];

if nargin > 6
    % Plot the new half tips
    plot(datetime(NewEvent_Date(NewEvent_mm~=bucket),'ConvertFrom','datenum'),...
        NewEvent_mm(NewEvent_mm~=bucket),'o','DisplayName','Auxiliary events')
    plot(datetime(Event_Date(Event_mm==0),'ConvertFrom','datenum'),...
        Event_mm(Event_mm==0),'xk','DisplayName','Zero intensity')
    set(gca,'Ylim',[-2 max(NewEvent_mm)+1])
    legend('off'); legend('show')
    legend('location','SouthWest'); legend('boxoff')
end

% Identify events from tips separated by more than the maximum time MaxT.
NewEventDiff = diff(NewEvent_Date) > MaxT;
NewEventDiff = cat(1,true,NewEventDiff);
% Array containing the event pointers.
indx = find(NewEventDiff);
% Number of elements per event.
n = diff(indx)-1;
n(end+1) = length(NewEvent_Date) - indx(end);
% Duration of the events in minutes.
D = (NewEvent_Date(indx(2:end)-1) - NewEvent_Date(indx(1:end-1)))*1440;
D(end+1) = NewEvent_Date(end) - NewEvent_Date(indx(end));
n1 = indx(n<1); % Index of events with only 1 point
fprintf('Number of rainfall events identified: %6i.\n',length(indx))
fprintf('Average duration of the events: %8.2f min.\n',mean(D(D>0)))
fprintf('Rainfall events consisting of 1 tip only: %6i.\n',length(n1))
% fprintf('\n')

%% FIT EVENTS AND AGGREGATING AT 1-min INTERVAL
% Build a 1 minute cumulative rainfall curve
DI = floor(min(Event_Date))*nd; % Initial date [day]
DF = ceil(max(Event_Date))*nd; % Final date [day]
NewDate_1min = (DI:DF)'; % Equally spaced time interval
CumP_1min = zeros(size(NewDate_1min)); % Initialise accumulation
Single_1min = zeros(size(NewDate_1min)); % Initialise single tip counting
biased = zeros(size(n)); % Initialise bias vector
bEvent = zeros(size(n)); % Initialise biased events counter.

h = waitbar(0,'Interpolating data...');
for i = 1:length(n)
    % Events with more than 2 points (fit a CS) [Wang et al, 2008].
    % Events with only 2 points (fit a line) [Ciach, 2003].
    % Events with only 1 points (distribute at a rate of 3 mm h^{-1}) [Wang et al, 2008].
    if n(i) >= 1
        % Relative time in seconds from the beggining of the event.
        x = (NewEvent_Date(indx(i)+(0:n(i)))-NewEvent_Date(indx(i)))*86400;
        % Cumulative rainfall during the event.
        y = cumsum(NewEvent_mm(indx(i):indx(i)+n(i)));
        if halves ~= false || halves ~= 0
            % Estimate initial point of the rainfall event.
            % Reduce half a sec only to ensure correct initial date calculation.
            x0 = bucket*(x(2)-x(1))/(y(2)-y(1))-0.5;
            xf = bucket*(x(end)-x(end-1))/(y(end)-y(end-1));
            % Allocate only 1-half tip at the start and end of event.
            x = x + x0;
            y = y - bucket/2;
            y = cat(1,0,y,y(end)+bucket/2); x = cat(1,0,x,x(end)+xf);
            x = round(x);
            % Aggregating data at 1 min interval starting at :00.
            DI = max(DI,floor((NewEvent_Date(indx(i))-x0/86400)*nd)); % Initial date in [min]
            DF = ceil((NewEvent_Date(indx(i)+n(i))+xf/86400)*nd); % Final date  in [min]
            x1m = (DI:DF)' - NewEvent_Date(indx(i))*nd+x0/60; % Equally spaced time interval.
            x1m = round(60*x1m); % Convert to seconds
        else
            % x0 = -0.5; % Only to ensure correct initial date calculation.
            % Aggregating data at 1 min interval starting at :00.
            DI = max(DI,floor((NewEvent_Date(indx(i))+0.5/86400)*nd)); % Initial date in [min]
            DF = ceil(NewEvent_Date(indx(i)+n(i))*nd); % Final date  in [min]
            x1m = (DI:DF)' - NewEvent_Date(indx(i))*nd; % Equally spaced time interval.
            x1m = round(60*x1m); % Convert to seconds
        end
        % CS fitted to the current event and interpolated at 1-sec.
        % pp = spline(x,y); % yy = spline(x,y,xx);
        if halves ~= false || halves ~= 0
            % Set the estimated zero rate endpoints first derivatives to 0
            % [Sadler and Busscher, 1989]. 
            pp = spline(x,[0;y;0]);
        else
            % Set the endpoints second derivatives to 0 [Wang et al, 2008].
        	pp = csape(x,y,'second');
        end
        y1m = fnval(pp,x1m); % Cumulative rainfall at each x1m
        if halves ~= false || halves ~= 0
            % Zero rainfall rates at borders.
            y1m(1) = 0; y1m(end) = y1m(end-1);
        end
        r1m = [y1m(1);diff(y1m)]; % Rainfall rate at each x1m [mm min^{-1}]
        % Correction for negative intensities and biased volumes.
        [r2m,y2m,biased(i),bEvent(i)] = intCorrection(r1m,y,Lowint,halves,x,x1m);
        
        % Assemblement of the cumulative rainfall curve.
        CumP_1min(NewDate_1min>=DI & NewDate_1min<=DF) = ...
            CumP_1min(NewDate_1min>=DI & NewDate_1min<=DF) + y2m;
        CumP_1min(NewDate_1min>DF) = CumP_1min(NewDate_1min==DF);
        % Plot events if selected.
        if nargin > 6 || any(isnan(r2m)) %|| bEvent(i)>0
            % Calculations for event plots.
            % xx = (0:x(end))';
            % yy = fnval(pp,xx);
            % Aggregating data at 1 min interval
            % Alternatively starting at x(1).
            % x1m = 60*(1:floor(length(xx)/60))'-1; % Equally spaced time interval.
            % y1m = yy(60*(1:floor(length(xx)/60)));
            % r1m = 60*(yy(60*(1:floor(length(xx)/60)))-yy(60*(1:floor(length(xx)/60))-59));
            
            % Linear sections of the current event and interpolated at 1-sec.
            y3m = interp1(x,y,x1m,'linear',0); % Cumulative rainfall at each x1m
            if halves ~= false || halves ~= 0
                % Zero rainfall rates at borders.
                y3m(1)= 0; y3m(end) = y3m(end-1);
            end
            r3m = [y3m(1);diff(y3m)]; % Rainfall rate at each x1m [mm min^{-1}]
            % Correction for biased volumes.
            [r4m,y4m] = intCorrection(r3m,y,Lowint,halves,x,x1m);
            
            % Tip counting
            r5m = zeros(size(x1m)); % Initialise aggregation
            p = [y(1);diff(y)];
            if x(1)==x1m(1)
                j = 2; % Data counter
                r5m(1) = p(1);
            else
                j = 1; % Data counter
            end
            for itb = 2:length(x1m)
                % Agregate values.
                while j<=length(p) && x(j)>x1m(itb-1) && x(j)<=x1m(itb)
                    r5m(itb) = r5m(itb) + p(j);
                    j = j+1;
                end
            end
            y5m = cumsum(r5m); % Accumulation
            
            % Plot events.
            figure(199)
            subplot(2,1,1)
            plot(x,y,'o',x1m,y5m,':',x1m,y4m,'-.',x1m,y1m,'-.',x1m,y2m,'-.')
            set(gca,'Xlim',[-60 x(end)+60])
            title(['Rainfall event number ',num2str(i)])
            xlabel('Time in seconds from the beggining of the event [s]')
            ylabel('Cumulative rainfall [mm]')
            legend('Rain gauge tip','1-min tip counting','1-min linear corrected',...
                '1-min CS interpolated','1-min CS corrected',...
                'location','NorthWest')
            legend('boxoff')
            subplot(2,1,2)
            plot(x,60*[y(1);diff(y)],'o',x1m,60*r5m,':',x1m,60*r4m,'-.',x1m,60*r1m,'-.',x1m,60*r2m,'-.')
            set(gca,'Xlim',[-60 x(end)+60])
            xlabel('Time in seconds from the beggining of the event [s]')
            ylabel('Rainfall rate [mm h^{-1}]')
            legend('Rain gauge tip','1-min tip counting','1-min linear corrected',...
                '1-min CS interpolated','1-min CS corrected',...
                'location','NorthWest')
            legend('boxoff')
            % Add breakpoint to check plots.
        end
    else
        % Aggregating data at 1 min interval starting at :00.
        x0 = NewEvent_mm(indx(i))/Meanint*60-1; % Time interval in [min]
        xf = NewEvent_Date(indx(i))*nd; % Final date in [min]
        x = (xf-x0*nd/1440:xf)'; % Equally spaced divided tip
        DI = floor((xf-x0*nd/1440)); % Initial date in [min]
        DF = ceil(xf); % Final date in [min]
        x1m = (DI:DF)'; % Equally spaced time interval
        y = NewEvent_mm(indx(i))*ones(size(x))/(x0+1);
        % Tip counting
        r1m = zeros(size(x1m)); % Initialise aggregation
        if x(1)==x1m(1)
            j = 2; % Data counter
            r1m(1) = y(1);
        else
            j = 1; % Data counter
        end
        for itb = 2:length(x1m)
            % Aggregate values.
            while j<=length(y) && x(j)>x1m(itb-1) && x(j)<=x1m(itb)
                r1m(itb) = r1m(itb) + y(j);
                j = j+1;
            end
        end
        y1m = cumsum(r1m); % Cumulative rainfall at each x1mim
        % Assemblement of the cumulative rainfall curve.
        CumP_1min(NewDate_1min>=DI & NewDate_1min<=DF) = ...
            CumP_1min(NewDate_1min>=DI & NewDate_1min<=DF) + y1m;
        CumP_1min(NewDate_1min>DF) = CumP_1min(NewDate_1min==DF);
        % Assemblement of the single tip rainfall vector.
        Single_1min(NewDate_1min>=DI & NewDate_1min<=DF) = ...
            Single_1min(NewDate_1min>=DI & NewDate_1min<=DF) + y1m;
        Single_1min(NewDate_1min>DF) = Single_1min(NewDate_1min==DF);
    end
    waitbar(i/length(n))
end
close(h);
fprintf('Maximum bias corrected in event interpolation: %5.2f%%.\n',100*nanmax(biased))
fprintf('%2i event(s) with bias >25%% interpolated linearly.\n',nansum(bEvent))
fprintf('\n')
% drawnow
% delete(findall(0,'Type','figure'))

%% PREPARE THE DATA VECTORS AT THE SPECIFIED SCALE

% Equally spaced time interval
NewDate = (NewDate_1min(1):scale:NewDate_1min(end))';
if halves ~= false || halves ~= 0
    % Rainfall rate at scale interval obtained from fitted cumulative rainfall.
    NewP = [bucket/2;CumP_1min(scale+1:scale:end) - CumP_1min(1:scale:end-scale)];
    % Aggregate single tip counting at scale interval.
    Single = [bucket/2;Single_1min(scale+1:scale:end) - Single_1min(1:scale:end-scale)];
else
    % Rainfall rate at scale interval obtained from fitted cumulative rainfall.
    NewP = [0;CumP_1min(scale+1:scale:end) - CumP_1min(1:scale:end-scale)];
    % Aggregate single tip counting at scale interval.
    Single = [0;Single_1min(scale+1:scale:end) - Single_1min(1:scale:end-scale)];
end
% Cumulative rainfall at scale interval.
CumP = cumsum(NewP);
% Correct numerical errors in the calculations of zero rain rates.
NewP(round(NewP,8)==0) = 0; Single(round(Single,8)==0) = 0;
% Cut the vectors to the actual initial and final date.
nd = 1440/scale; % Number of intervals per day
DI = ceil(min(Event_Date)*nd)*scale; % Initial date in [min]
DF = ceil(max(Event_Date)*nd)*scale; % Final date in [min]
CumP(NewDate<DI | NewDate>DF) = [];
NewP(NewDate<DI | NewDate>DF) = [];
Single(NewDate<DI | NewDate>DF) = [];
NewDate(NewDate<DI | NewDate>DF) = [];
NewDate = NewDate/1440; % Rescale the date
% Example (or validation) for 1 minute aggregation.
% NewDate = NewDate_1min/1440;
% CumP = CumP_1min;
% NewP = [CumP(1);diff(CumP)];

% Restoring dates to date format
NewDate = datetime(NewDate,'ConvertFrom','datenum');
% Placing data gaps again in the aggregated vectors.
for i = 1:size(Voids,1)
    CumP(NewDate>Voids(i,1) & NewDate<Voids(i,2)) = NaN;
    NewP(NewDate>Voids(i,1) & NewDate<Voids(i,2)) = NaN;
    Single(NewDate>Voids(i,1) & NewDate<Voids(i,2)) = NaN;
end
if nargout == 1
    NewDate = [datenum(NewDate),NewP,CumP,Single];
end
fprintf('Rainfall volume before aggregation: %8.2f mm.\n',nansum(Event_mm))
fprintf('Rainfall volume after aggregation: %8.2f mm.\n',nansum(NewP))
fprintf('\n')


%% AUXILIARY FUNCTIONS
function [NewDate_1min,NewP_1min] = AggregateEvents(Event_Date,Event_mm)
% Agregate rainfall at 1-min intervals.
k = length(Event_mm); % Length of input data
nd = 1440; % Number of minutes per day or numeric value of 1 minute: 1/1440
% Build a 1 minute cumulative rainfall curve
DI = floor(min(Event_Date))*nd; % Initial date [day]
DF = ceil(max(Event_Date))*nd; % Final date [day]
NewDate_1min = (DI:DF)'; % Equally spaced time interval
n = length(NewDate_1min); % Number of 1-min intervals
NewP_1min = zeros(size(NewDate_1min)); % Initialise aggregation
if nd*Event_Date(1)==NewDate_1min(1)
    j = 2; % Data counter
    NewP_1min(1) = Event_mm(1);
else
    j = 1; % Data counter
end
for i = 2:n
    % Aggregate values.
    while j<=k && nd*Event_Date(j)<=NewDate_1min(i) % && nd*Event_Date(j)>NewDate_1min(i-1)
        NewP_1min(i) = NewP_1min(i) + Event_mm(j);
        j = j+1;
    end
end
NewDate_1min = NewDate_1min/nd; % Rescale the date
% Delete zero events to help process relevant data only.
NewDate_1min(NewP_1min==0) = [];
NewP_1min(NewP_1min==0) = [];
fprintf('Routine for aggregating tips at 1-min time interval.\n')
fprintf('New number of data points: %4i.\n',length(NewP_1min))
fprintf('Rainfall volume before aggregation: %8.2f mm.\n',nansum(Event_mm))
fprintf('Rainfall volume after aggregation: %8.2f mm.\n',nansum(NewP_1min))
fprintf('\n')

function [NewEvent_Date,NewEvent_mm] = MergeEvents(Event_Date,Event_mm,MinT)
% Delete tips for small periods.
% Calculate the time between tips.
Diff_Event_Date = diff(Event_Date);
% Identify tips separated by less than the minimum time MinT.
EventDiff = Diff_Event_Date <= MinT;
% Modified variables to process.
NewEvent_Date = Event_Date;
NewEvent_mm = Event_mm;
i = 0; j = 0; % Initialise counters
while any(EventDiff)
    % Reassign rainfall to the following interval.
    j = 1; n =length(NewEvent_Date);
    for i = 2:n
        j = j+1;  % Index for NewEvent_mm.
        % When the time between tips is less than MaxT.
        if EventDiff(i-1)
            % Aggregate the tip volume to the following time stamp.
            NewEvent_mm(j) = NewEvent_mm(j)+NewEvent_mm(j-1);
            % Include these data in the rainfall tip time series.
            NewEvent_mm(j-1) = [];
            NewEvent_Date(j-1) = [];
            j = j-1;
        end
    end
    % Calculate the time between tips.
    Diff_Event_Date = diff(NewEvent_Date);
    % Identify tips separated by less than the minimum time MinT.
    EventDiff = Diff_Event_Date <= MinT;
end
fprintf('Routine for merging tips occurring faster than MinT = %6.2f seconds.\n',MinT*86400)
fprintf('Number of tips removed: %4i.\n',i-j)
fprintf('Rainfall volume before merging: %8.2f mm.\n',nansum(Event_mm))
fprintf('Rainfall volume after merging: %8.2f mm.\n',nansum(NewEvent_mm))
fprintf('\n')

function [NewEvent_Date,NewEvent_mm] = DivideEvents(Event_Date,Event_mm,MaxT)
% Add additional tips for long periods [Wang et al, 2008].
% Calculate the time between tips.
Diff_Event_Date = diff(Event_Date);
% Identify tips separated by more than the maximum time MaxT.
EventDiff = Diff_Event_Date > MaxT;
% Redistribute rainfall over relatively long periods but lower than MaxT.
% but greater than half MaxT.
HalfEventDiff = Diff_Event_Date > MaxT/2;
% Modified variables to process.
NewEvent_Date = Event_Date;
NewEvent_mm = Event_mm;
j = 1;
for i = 2:length(EventDiff)
    j = j+1;  % Index for NewEvent_mm.
    % When the time between tips is between half and one complete MaxT.
    % Be aware that ~EventDiff(i-2) may conflict with i = 2.
    if HalfEventDiff(i-1) && ~EventDiff(i-1) && (~EventDiff(i) || ~EventDiff(i-2))
        % Divide the following tip volume in two and assign a time stamp.
        Halftip_mm = Event_mm(i)/2;
        t0 = Event_Date(i) - Diff_Event_Date(i-1)/2;
        % Include these data in the rainfall tip time series.
        NewEvent_Date = cat(1,NewEvent_Date(1:j-1),t0,NewEvent_Date(j:end));
        NewEvent_mm = cat(1,NewEvent_mm(1:j-1),Halftip_mm,Halftip_mm,NewEvent_mm(j+1:end));
        j = j+1;
    end
end
fprintf('Routine for spreading tips occurring between %6.2f and %6.2f minutes.\n',MaxT*1440/2,MaxT*1440)
fprintf('Number of tips added: %4i.\n',j-i)
fprintf('Rainfall volume before spreading: %8.2f mm.\n',nansum(Event_mm))
fprintf('Rainfall volume after spreading: %8.2f mm.\n',nansum(NewEvent_mm))
fprintf('\n')

function [r2m,y2m,biased,bEvent] = intCorrection(r1m,y,Lowint,halves,x,x1m)
% Only keep rain rates above the threshold.
r2m = r1m;
bEvent = 0;
biased = abs(y(end)-sum(r1m(r1m>0)))/y(end);
if biased > 0.25
    bEvent = 1;
    % Use linear interpolation instead
    y2m = interp1(x,y,x1m,'linear',0); % Cumulative rainfall at each x1m
    if halves ~= false || halves ~= 0
        % Zero rainfall rates at borders.
        y2m(1) = 0; y2m(end) = y2m(end-1);
    end
    r2m = [y2m(1);diff(y2m)]; % Rainfall rate at each x1m [mm min^{-1}]
    biased = abs(y(end)-sum(r2m(r2m>0)))/y(end);
    % Correction for negative intensities and biased volumes.
end
iter = 0;
while iter<= 10 && (abs(y(end)-sum(r2m))>Lowint || any(round(r2m(r2m~=0),8) < Lowint))
    iter = iter + 1;
    r2m(r2m < 0) = 0; % Set the negative rainfall rates to zero
    r2m(r2m > 0 & r2m < Lowint) = Lowint; % Replace the lowest rates
    r2m(r2m >= Lowint) = r2m(r2m >= Lowint)*...
        (y(end)-sum(r2m(r2m < Lowint)))/...
        (sum(r2m)-sum(r2m(r2m < Lowint))); % Correction of biased rainfall
end
% Calculate the cumulative rainfall curve.
y2m = cumsum(r2m);
% Correct slight differences in totals.
y2m(end) = y(end);
if halves ~= false || halves ~= 0
    % Zero rainfall rates at borders.
    y2m(end-1) = y(end);
end
% Recalculate the intensities.
r2m = [y2m(1);diff(y2m)];