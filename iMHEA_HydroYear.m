function [mondata_hy,yeardata_hy] = iMHEA_HydroYear(mondata,years,ax)
%iMHEA Reorganises data into hydrological years.
% [data_hy] = iMHEA_HydroYear(mondata,ax)
%
% Input:
% mondata = Matrix of monthly data [12 columns, one for each month],
%           months are defined as 1:Jan to 12:Dec
% years   = Vector of 1 column with the correspondent years for mondata
% ax      = Vector of 2 elements dd/mm/yyyy [date format]
%
% Output:
% mondata_hy  = Matrix of monthly data [12 columns, one for each month],
%               using months defined in input ax.
% yeardata_hy = Matrix of yearly data [1 column:hydro_year, 2 column: data]
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in November, 2018
% Last edited in November, 2018

%% REORGANISE DATA

% Months and years that define the hydrological year
refmon_hy = month(ax(1));
refyrs_hy = year(ax);
% Regorganise the matrix
mondata_hy = [mondata(and(years(:,1)>=refyrs_hy(1),years(:,1)<=refyrs_hy(2)-1),refmon_hy:12),mondata(and(years(:,1)>=refyrs_hy(1)+1,years(:,1)<=refyrs_hy(2)),1:refmon_hy-1)];
yeardata_hy(:,1) = (refyrs_hy(1)+1:refyrs_hy(2))';
yeardata_hy(:,2) = nansum(mondata_hy,2);