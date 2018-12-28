function [DataHRes] = iMHEA_WorkflowPairMamanteo(DataHRes1,DataHRes2)
%iMHEA workflow for paired catchments, mamanteo exercise.
% [DataHRes] = iMHEA_WorkflowPair(DataHRes1,DataHRes2) assimilates the
% processed precipitation and discharge data to provide a pairwise
% assessment.
%
% Input:
% DataHres1 = Catchment 1 [Date, P, Q] at max resolution, e.g. 5min.
% DataHres1 = Catchment 2 [Date, P, Q] at max resolution, e.g. 5min.
%             P = Average Precipitation [mm].
%             Q = Normalised Discharge [l/s/km2].
%
% Output:
% DataHres = [Date, P1, Q1, P2, Q2] at max resolution, e.g. 5min.
% Data1day = [Date, P1, Q1, BQ1 P2, Q2, BQ2] at 1 day resolution.
% Data1hr  = [Date, P1, Q1, P2, Q2] at 1 hour resolution.
%             P1, Q1, BQ1 = Catchment 1.
%             P1, Q1, BQ2 = Catchment 2.
%             BQ = Normalised Baseflow [l/s/km2].
% Indices  = Matrix of hydrological indices from streamflow.
% Climate  = Matrix of climate indices from precipitation.
% PM = Monthly precipitation (mm) per month number [Jan=1, Dec=12].
% QM = Monthly Mean Daily flow (l/s) per month number [Jan=1, Dec=12].
% IDCi  = Maximum Intensity - Duration Curve [mm/h v time] for catchment i.
% FDCi  = Flow Duration Curve [l/s v %] for catchment i.
% Comparative plots.
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in November, 2017
% Last edited in February, 2018

%% START PROCESS
fprintf('\n')
fprintf('iMHEA WORKFLOW FOR PAIRED CATCHMENTS')
fprintf('\n')

%% CHECK IF PRECIPITATION AND DISCHARGE DATA IN BOTH CATCHMENTS MATCH.
% Determine aggregation interval
int_HRes1 = diff(DataHRes1(:,1))*1440;
int_HRes2 = diff(DataHRes2(:,1))*1440;
% Worst discharge interval defines the max resolution
int_HRes = round(max([nanmedian(int_HRes1),nanmedian(int_HRes2)]));
nd = 1440/int_HRes; % Number of intervals per day
% Interpolate data at max resolution
if round(nanmedian(int_HRes1)) > round(nanmedian(int_HRes2))
    % If interval in catchment 1 is larger than in catchment 2,
    % aggregate catchment 2
    [PrecHRes2] = iMHEA_Aggregation(DataHRes2(:,1),DataHRes2(:,2),int_HRes);
    [DiscHRes2] = iMHEA_Average(DataHRes2(:,1),DataHRes2(:,3),int_HRes);
    PrecHRes2(:,3:end) = [];
    PrecHRes2(:,3) = DiscHRes2(:,2);
    DataHRes2 = PrecHRes2;
elseif round(nanmedian(int_HRes1)) < round(nanmedian(int_HRes2))
    % If interval in catchment 2 is larger than in catchment 1,
    % aggregate catchment 1
    [PrecHRes1] = iMHEA_Aggregation(DataHRes1(:,1),DataHRes1(:,2),int_HRes);
    [DiscHRes1] = iMHEA_Average(DataHRes1(:,1),DataHRes1(:,3),int_HRes);
    PrecHRes1(:,3:end) = [];
    PrecHRes1(:,3) = DiscHRes1(:,2);
    DataHRes1 = PrecHRes1;
end

%% FILL PRECIPITATION GAPS BETWEEN CATCHMENTS
% Fill Precipitation gaps between catchments 1 and 2.
% (not recommended for discharge).
[DataHRes] = iMHEA_FillGaps(DataHRes1(:,1),DataHRes1(:,2),DataHRes2(:,1),DataHRes2(:,2));

%% COMPILE AVERAGE PRECIPITATION AND DISCHARGE DATA IN A SINGLE MATRIX
% Extend start and end of vectors to cover the same date period.
date_start = round(DataHRes(1,1)*nd);
date_end   = round(DataHRes(end,1)*nd);
% Create high resolution matrix
DataHRes(:,4:5) = nan;
DataHRes(:,1) = date_start:date_end;
DataHRes(:,4) = DataHRes(:,3); DataHRes(:,3) = nan;
% Add NaN values at the start and end of Q and P vector to even sizes.
DataHRes(ismember(DataHRes(:,1),round(DataHRes1*nd)),3) = DataHRes1(:,3);
DataHRes(ismember(DataHRes(:,1),round(DataHRes2*nd)),5) = DataHRes2(:,3);
% Rescale the date
DataHRes(:,1) = DataHRes(:,1)/nd;

% %% AGGREGATE DATA AT 1 HOUR AND 1 DAY TEMPORAL RESOLUTIONS
% % Aggregate precipitation
% [Data1day] = iMHEA_Aggregation(DataHRes(:,1),DataHRes(:,2),1440);
% Data1day(:,3:end) = nan;
% [~,Data1day(:,5)] = iMHEA_Aggregation(DataHRes(:,1),DataHRes(:,4),1440);
% [Data1hr] = iMHEA_Aggregation(DataHRes(:,1),DataHRes(:,2),60);
% Data1hr(:,3:end) = nan;
% [~,Data1hr(:,4)] = iMHEA_Aggregation(DataHRes(:,1),DataHRes(:,4),60);
% 
% % Average stream flows
% [~,Data1day(:,3)] = iMHEA_Average(DataHRes(:,1),DataHRes(:,3),1440);
% [~,Data1day(:,6)] = iMHEA_Average(DataHRes(:,1),DataHRes(:,5),1440);
% [~,Data1hr(:,3)] = iMHEA_Average(DataHRes(:,1),DataHRes(:,3),60);
% [~,Data1hr(:,5)] = iMHEA_Average(DataHRes(:,1),DataHRes(:,5),60);
% 
% % Calculate daily baseflow
% [~,Data1day(:,4)] = iMHEA_BaseFlowUK(DataHRes(:,1),DataHRes(:,3));
% [~,Data1day(:,7)] = iMHEA_BaseFlowUK(DataHRes(:,1),DataHRes(:,5));
% 
% %% CALCULATE HYDROLOGICAL AND CLIMATE INDICES
% Date1 = datetime(DataHRes(:,1),'ConvertFrom','datenum');
% [Indices,Climate,PM,QM,IDC1,FDC1,IDC2,FDC2] = iMHEA_Pair(Date1,DataHRes(:,2),DataHRes(:,3),1,Date1,DataHRes(:,4),DataHRes(:,5),1);

%% END PROCESS
fprintf('\n')
fprintf('Calculations finished')
fprintf('\n')