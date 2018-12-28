function [NewDate1,NewP1,NewDate2,NewP2] = iMHEA_FillGaps(Date1,P1,Date2,P2,cutend,varargin)
%iMHEA Data gaps fill.
% [NewDate1,NewP1,NewDate2,NewP2] = iMHEA_FillGaps(Date1,P1,Date1,P2,cutend,flag).
%
% Input:
% Date1, Date2 = dd/mm/yyyy hh:mm:ss [date format].
%                Both should be at the same temporal scale.
% P1, P2       = Precipitation [mm].
% cutend = Cut vectors not to fill gaps after the end [default: false].
% flag   = leave empty NOT to run the data voids assessment and plots.
%
% Output:
% NewDate1, 2   = dd/mm/yyyy hh:mm:ss [date format].
% NewP1, NewP2  = Filled precipitation data [mm].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in March, 2016
% Last edited in November, 2017

%% INITIALISE VARIABLES
fprintf('\n')
fprintf('FILL DATA GAPS BY LINEAR REGRESSION.\n')
fprintf('\n')
if nargin < 5 || isempty(cutend)
    % Do not cut end of the vectors by default.
    cutend = false;
end

%% CHECK IF DATA ARE AT THE SAME TEMPORAL RESOLUTION
scale1 = diff(datenum(Date1))*1440;
scale2 = diff(datenum(Date2))*1440;
if nanmedian(scale1) > nanmedian(scale2)
    fprintf('Input data are not at the same temporal scale.\n')
    scale = round(nanmedian(scale1)); % Same temporal resolution
    [Date1,P1] = iMHEA_Aggregation(Date1,P1,scale);
    [Date2,P2] = iMHEA_Aggregation(Date2,P2,scale);
elseif nanmedian(scale2) > nanmedian(scale1)
    fprintf('Input data are not at the same temporal scale.\n')
    scale = round(nanmedian(scale2)); % Same temporal resolution
    [Date1,P1] = iMHEA_Aggregation(Date1,P1,scale);
    [Date2,P2] = iMHEA_Aggregation(Date2,P2,scale);
else
    scale = round(nanmedian(scale1)); % Same temporal resolution
end

%% VOID ASSESSMENT
% Run data gap assessment and print inventory.
fprintf('Data gap assessment of P1.\n')
[~] = iMHEA_Voids(Date1,P1,1);
fprintf('Data gap assessment of P2.\n')
[~] = iMHEA_Voids(Date2,P2,1);

%% CREATE UNIFIED DATE VECTOR AND ASSIGN CORRESPONDENT INPUT DATA
% Numeric value of 1 minute: 1/1440
nd = 1440/scale; % Number of intervals per day
% Convert Dates to integers to avoid precision errors
Date1 = round(nd*datenum(Date1));
Date2 = round(nd*datenum(Date2));
% Define initial and end dates and create single vector
DI = min(Date1(1),Date2(1));
DF = max(Date1(end),Date2(end));
NewDate = (DI:DF)';
% Assign data when they correspond
NewP1 = nan(size(NewDate));
NewP1(ismember(NewDate,Date1)) = P1;
NewP2 = nan(size(NewDate));
NewP2(ismember(NewDate,Date2)) = P2;
% Optionally, cut vectors not to fill gaps after the end.
if cutend
    % Identify the last non-NaN data in both vectors.
    indexnP1 = find(~isnan(P1),1,'last');
    indexnP2 = find(~isnan(P2),1,'last');
    indexndate = min(Date1(indexnP1),Date2(indexnP2));
    % Cut vectors after the minimum of the identified dates.
    NewP1(NewDate>indexndate) = [];
    NewP2(NewDate>indexndate) = [];
    NewDate(NewDate>indexndate) = [];
end

%% TEST IF OVERLAPPING DATA EXIST
% Extract all the sections where NaN data exist in any of the vectors.
auxP1 = NewP1; auxP2 = NewP2;
auxP1(isnan(NewP1)|isnan(NewP2)) = [];
auxP2(isnan(NewP1)|isnan(NewP2)) = [];
% Check if any of the vectors is empty
if isempty(auxP1) || length(auxP1) == 1
    fprintf('There is not date coincidence between the input data.\n')
    fprintf('\n')
    % Restore data if cut before
    if cutend
        NewDate = (DI:DF)';
        % Assign data when they correspond
        NewP1 = nan(size(NewDate));
        NewP1(ismember(NewDate,Date1)) = P1;
        NewP2 = nan(size(NewDate));
        NewP2(ismember(NewDate,Date2)) = P2;
    end
    % Produce outputs
    if nargout == 1
        NewDate1 = [NewDate/nd,NewP1,NewP2];
    else
        % Otherwise transform dates to date format
        NewDate1 = datetime(NewDate/nd,'ConvertFrom','datenum');
        NewDate2 = datetime(NewDate/nd,'ConvertFrom','datenum');
    end
    iMHEA_Plot3(datetime(NewDate/nd,'ConvertFrom','datenum'),NewP1,NewP2)
    return
end

%% FILL DATA GAPS
% Calculate cumulative rainfall curves.
auxCumP1 = cumsum(auxP1);  
auxCumP2 = cumsum(auxP2);
[R,M,~] = regression(auxCumP1',auxCumP2');
% Fill gaps only if the correlation is almost perfect.
if R < 0.99
    fprintf('The correlation is not significant as to fill the data, with R2 = %6.4f.',R)
    fprintf('\n')
    figure
    plotregression(auxCumP1,auxCumP2,'Regression')
    % Restore data if cut before
    if cutend
        NewDate = (DI:DF)';
        % Assign data when they correspond
        NewP1 = nan(size(NewDate));
        NewP1(ismember(NewDate,Date1)) = P1;
        NewP2 = nan(size(NewDate));
        NewP2(ismember(NewDate,Date2)) = P2;
    end
    % Produce outputs
    if nargout == 1
        NewDate1 = [NewDate/nd,NewP1,NewP2];
        % NewDate1 = [NewDate1,NewP1,NewP2,cumsum(NewP1),cumsum(NewP2)];
    else
        % Otherwise transform dates to date format
        NewDate1 = datetime(NewDate/nd,'ConvertFrom','datenum');
        NewDate2 = datetime(NewDate/nd,'ConvertFrom','datenum');
    end
    return
end
if nargin >= 6
    % Plot the regression.
    figure
    plotregression(auxCumP1,auxCumP2,'Regression')
end
NewP1(isnan(NewP1)) = NewP2(isnan(NewP1))/M;
NewP2(isnan(NewP2)) = NewP1(isnan(NewP2))*M;

%% RESTORE THE DATA AND THE END OF THE VECTORS
if cutend
    % Assign data when they correspond
    NewP1(ismember(NewDate,Date1)) = P1;
    NewP2(ismember(NewDate,Date2)) = P2;
end

%% GENERATE OUTPUTS
% Restore Dates from integers made to avoid precision errors
Date1 = datetime(Date1/nd,'ConvertFrom','datenum');
Date2 = datetime(Date2/nd,'ConvertFrom','datenum');
NewDate1 = datetime(NewDate/nd,'ConvertFrom','datenum');
NewDate2 = datetime(NewDate/nd,'ConvertFrom','datenum');

%% PRINT RESULTS
fprintf('Rainfall volumes before filling gaps: %8.2f and %8.2f mm.\n',nansum(P1),nansum(P2))
fprintf('Rainfall volumes after filling gaps: %8.2f and %8.2f mm.\n',nansum(NewP1),nansum(NewP2))
fprintf('\n')

% Plot regressions
if nargin >= 6
    % Extract all the sections where NaN data exists in any of the vectors.
    auxP1 = NewP1; auxP2 = NewP2;
    auxP1(isnan(NewP1)|isnan(NewP2)) = [];
    auxP2(isnan(NewP1)|isnan(NewP2)) = [];
    % Calculate cumulative rainfall curves.
    auxCumP1 = cumsum(auxP1);  
    auxCumP2 = cumsum(auxP2);
    % Plot the results
    figure
    plotregression(auxCumP1,auxCumP2,'Regression')
    
    % Plot time series
    figure
    subplot(2,1,1)
    plot(Date1,P1,Date2,P2)
    xlabel('Date')
    ylabel('Precipitation [mm]')
    title('Before gap fill')
    box on
    subplot(2,1,2)
    plot(NewDate1,NewP1,NewDate2,NewP2)
    xlabel('Date')
    ylabel('Precipitation [mm]')
    title('After gap fill')
    box on
end

if nargout == 1
    NewDate1 = [NewDate/nd,NewP1,NewP2];
end