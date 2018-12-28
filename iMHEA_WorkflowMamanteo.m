function [DataHRes] = iMHEA_WorkflowMamanteo(AREA,DateQ,Q,bucket,varargin)
%iMHEA entire workflow for a catchment, mamanteo exercise.
% [DataHRes] =
% iMHEA_Workflow(AREA,DateQ,Q,bucket,DateP1,P1,DateP2,P2,...) 
% processes the raw precipitation and discharge data to get processed data.
%Rainfal
% Input:
% AREA   = Catchment area [km2].
% DateQ  = dd/mm/yyyy hh:mm:ss [date format] for discharge.
% Q      = Discharge [l/s].
% bucket = Rain gauge bucket volume [mm].
% DatePi = dd/mm/yyyy hh:mm:ss [date format] for rain gauge i.
% Pi     = Precipitation [mm] from rain gauge i.
%
% Output:
% DataHres = [Date, P, Q] at max resolution, e.g. 5min.
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in November, 2017
% Last edited in February, 2018

%% START PROCESS
fprintf('\n')
fprintf('iMHEA WORKFLOW FOR CATCHMENT %s',inputname(1))
fprintf('\n')

%% DETERMINE NUMBER OF RAIN GAUGES
nrg = floor((nargin-4)/2);
if nrg>=1 && rem(nrg,1)==0
    fprintf('The input data includes: %4i rain gauges.\n',nrg)
else
    fprintf('There are incomplete rain gauge data to process.\n')
    DataHRes = [];
    fprintf('\n')
    return
end

%% DEPURATE PRECIPITATION EVENTS
NewEvent_mm = cell(nrg,1);
for i = 1:nrg
    [NewEvent_mm{i}] = iMHEA_Depure(varargin{2*i-1},varargin{2*i});
end

%% AGGREGATE PRECIPITATION DATA TO MATCH DISCHARGE INTERVAL
% Determine discharge interval
int_HRes = diff(datenum(DateQ))*1440;
int_HRes = round(nanmedian(int_HRes)); % Worst discharge interval defines the max resolution
nd = 1440/int_HRes; % Number of intervals per day
% Interpolate each rain gauge data at max resolution using the CS algorithm
PrecHRes = cell(nrg,1);
for i = 1:nrg
    [PrecHRes{i}] = iMHEA_AggregationCS(varargin{2*i-1},NewEvent_mm{i},int_HRes,bucket);
end

%% FILL PRECIPITATION GAPS AND OBTAIN SINGLE MATRIX
if nrg > 1
    % Fill Precipitation gaps between all combinations of rain gauges
    combinations = combnk(1:nrg,2);
    combin_index = size(combinations,1);
    PrecHResFill = cell(combin_index,1);
    c = nan(combin_index,1);
    d = nan(combin_index,1);
    for i = 1:combin_index
        a = PrecHRes{combinations(i,1)};
        b = PrecHRes{combinations(i,2)};
        [PrecHResFill{i}] = iMHEA_FillGaps(a(:,1),a(:,2),b(:,1),b(:,2));
        c(i) = PrecHResFill{i}(1,1);
        d(i) = PrecHResFill{i}(end,1);
    end
    % Extend start and end of vectors to cover the same date period.
    date_start = round(min(c)*nd);
    date_end   = round(max(d)*nd);
    % Create high resolution matrix
    DateP_HRes = (date_start:date_end)';
    Precp_Fill_Compiled = nan(length(DateP_HRes),2*combin_index);
    for i = 1:combin_index
        % Compile precipitation data in a single matrix.
        DateAux = round(PrecHResFill{i}(:,1)*nd);
        Precp_Fill_Compiled(ismember(DateP_HRes,DateAux),2*i-1:2*i) = PrecHResFill{i}(:,2:3);
    end
    % Obtain average Precipitation
    P_HRes = nanmean(Precp_Fill_Compiled,2);
    % Rescale the date
    DateP_HRes = DateP_HRes/nd;
else
    % Obtain average Precipitation and create high resolution matrix
    DateP_HRes = PrecHRes{1}(:,1);
    P_HRes = PrecHRes{1}(:,2);
end

%% AVERAGE DISCHARGE DATA TO THE MAXIMUM RESOLUTION
[DateQ_HRes,Q_HRes] = iMHEA_Average(DateQ,Q,int_HRes);
DateQ_HRes = datenum(DateQ_HRes);

%% COMPILE AVERAGE PRECIPITATION AND DISCHARGE DATA IN A SINGLE MATRIX
% Extend start and end of vectors to cover the same date period
date_start = round(min([DateP_HRes(1),DateQ_HRes(1)])*nd);
date_end   = round(max([DateP_HRes(end),DateQ_HRes(end)])*nd);
% Create high resolution matrix
DateAux = date_start:date_end;
DataHRes = nan(length(DateAux),3+nrg);
DataHRes(:,1) = DateAux;
DataHRes(ismember(DataHRes(:,1),round(DateP_HRes*nd)),2) = P_HRes;
DataHRes(ismember(DataHRes(:,1),round(DateQ_HRes*nd)),3) = Q_HRes;
for i=1:nrg
    DataHRes(ismember(DataHRes(:,1),round(DateP_HRes*nd)),3+i) = Precp_Fill_Compiled(:,i);
end
% Rescale the date
DataHRes(:,1) = DataHRes(:,1)/nd;
% Normalise the flow
DataHRes(:,3) = DataHRes(:,3)/AREA;

% %% AGGREGATE DATA AT 1 HOUR AND 1 DAY TEMPORAL RESOLUTIONS
% % Aggregate precipitation
% [Data1day] = iMHEA_Aggregation(DataHRes(:,1),DataHRes(:,2),1440);
% [Data1hr] = iMHEA_Aggregation(DataHRes(:,1),DataHRes(:,2),60);
% Data1day(:,3:end) = [];
% Data1hr(:,3:end) = [];
% % Average stream flows
% [~,Data1day(:,3)] = iMHEA_Average(DataHRes(:,1),DataHRes(:,3)/AREA,1440);
% [~,Data1hr(:,3)] = iMHEA_Average(DataHRes(:,1),DataHRes(:,3)/AREA,60);
% % Calculate daily baseflow
% [~,Data1day(:,4)] = iMHEA_BaseFlowUK(DataHRes(:,1),DataHRes(:,3)/AREA);

%% PLOT RESULTING TIME SERIES
iMHEA_Plot3(datetime(DataHRes(:,1),'ConvertFrom','datenum'),DataHRes(:,2),DataHRes(:,3),DataHRes(:,4),DataHRes(:,5))

% %% CALCULATE HYDROLOGICAL AND CLIMATE INDICES
% % Tranform to date format.
% DateFormatHRes = datetime(DataHRes(:,1),'ConvertFrom','datenum');
% % Calculate indices.
% [Climate,Indices] = iMHEA_IndicesTotal(DateFormatHRes,DataHRes(:,2),DataHRes(:,3),AREA,1);
% % Normalise discharge data [l/s/km2].
% DataHRes(:,3) = DataHRes(:,3)/AREA;

%% END PROCESS
fprintf('\n')
fprintf('Calculations finished')
fprintf('\n')