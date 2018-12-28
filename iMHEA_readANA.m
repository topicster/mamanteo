function [HQ_1day,HQ_DDATE] = iMHEA_readANA(prefix,yrindx)
%iMHEA Reads flow data downloaded from ANA's SNIRH website.
% [ANA_CHO_HQ_1day,ANA_CHO_HQ_DDATE] = iMHEA_readANA(prefix,yrindx)
%
% Input:
% prefix   = filename prefix used in all the files to be read [string]
% yrindx   = year index to read the data [vector double]
%            if range is a row vector of size [1x2] use them as ini-end
%            if range is a column vector or single value use it as specific
%
% Output:
% HQ_1day  = time series of daily flows [usually m3/s]
% HQ_DDATE = corresponding date dd/mm/yyyy [date format]
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in November, 2018
% Last edited in November, 2018

%% Procedure

% Check specific years or create a vector given the range
[m,n] = size(yrindx);
if n==2 && m==1
    years = yrindx(1):yrindx(2);
else
    years = yrindx(:);
end

% Initialise the variables
N = length(years);
HQ_DDATE = datetime([1 1 1]);
HQ_1day = nan(1,1);

% Read the files for each specified year and merge them
for i = 1:N
    % Define filename and read table
    filename = [prefix,'_',num2str(years(i)),'.csv'];
    datatable = readtable(filename);
    
    % Determine amount of data (files come in months of 31 days all)
    M = length(datatable{:,1});
    appenddate = datetime([repmat(years(i),M,1),repmat((1:M/31)',31,1),datatable{:,1}]);
    appendflow = datatable{:,3};
    % Discard inexistent values
    appenddate(isnan(appendflow)) = [];
    appendflow(isnan(appendflow)) = [];
    % Sort data by date
    [sortdate,sortindx] = sort(appenddate);
    sortflow = appendflow(sortindx);
    % Append the new read data
    HQ_DDATE = [HQ_DDATE;sortdate];
    HQ_1day = [HQ_1day;sortflow];
end

% Remove auxiliar points
HQ_DDATE(1) =[];
HQ_1day(1) =[];