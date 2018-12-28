function [Voids,NoVoids] = iMHEA_Voids(Date,Data,varargin)
%iMHEA Determining void intervals.
% [Voids,NoVoids] = iMHEA_Voids(Date,Data).
%
% Input:
% Date  = dd/mm/yyyy hh:mm:ss [date format].
% Data  = Precipitation [mm] or Discharge [l/s, m3/s, mm].
% flag1 = Leave empty not to print interval results.
% flag2 = Leave empty not to plot data inventory.
%
% Output:
% Voids   = A two column matrix containing: [Initial date, Final date].
% NoVoids = A two column matrix containing: [Initial date, Final date].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in May, 2014
% Last edited in November, 2017

%% INITIALISE VARIABLES
k = length(Data); % Length of input data
Date = datenum(Date);

v = 1; nv = 1;  % Voids counter
Voids = zeros(1,1); % Initialise matrix
NoVoids = zeros(1,1); % Initialise matrix
Date(end+1) = Date (end); % Extend the vector in one element

%% DETERMINE DATA GAPS

% h = waitbar(0,'Identifying voids in data...');
for jv = 1:k
    if isnan(Data(jv))
        v = v + 1;
        Voids(v,1) = Date(jv); % Initial date of void interval
        Voids(v,2) = Date(jv+1); % Final date of void interval
        if Voids(v,1) == Voids(v-1,2)
            v = v - 1; % Agregate continuous voids
            Voids(v,2) = Voids(v+1,2);
            Voids(v+1,:) = [];
        end
    else
        nv = nv + 1;
        NoVoids(nv,1) = Date(jv); % Initial date of no-void interval
        NoVoids(nv,2) = Date(jv+1); % Final date of no-void interval
        if NoVoids(nv,1) == NoVoids(nv-1,2)
            nv = nv - 1; % Agregate continuous interval
            NoVoids(nv,2) = NoVoids(nv+1,2);
            NoVoids(nv+1,:) = [];
        end
    end
%     waitbar(jv/k)
end
% close(h);

%% GENERATE OUTPUT VARIABLES
% Cut the auxuliar element in each vector
Voids(1,:) = []; v = v-1;
NoVoids(1,:) = []; nv = nv-1;
if v==0
    Voids = [];
end
Voids = datetime(Voids,'ConvertFrom','datenum');
NoVoids = datetime(NoVoids,'ConvertFrom','datenum');

%% PRINT RESULTS
if nargin >= 3
    fprintf('\n')
    fprintf('IDENTIFICATION OF DATA GAPS.\n')
    fprintf('\n')
    fprintf('Data:\tInitial Date\t\tFinal Date\n')
    for i = 1:nv
    fprintf('\t')
    fprintf(datestr(NoVoids(i,1)))
    fprintf('\t')
    fprintf(datestr(NoVoids(i,2)))
    fprintf('\n')
    end
    fprintf('\n')
    fprintf('Voids:\tInitial Date\t\tFinal Date\n')
    for i = 1:v
    fprintf('\t')
    fprintf(datestr(Voids(i,1)))
    fprintf('\t')
    fprintf(datestr(Voids(i,2)))
    fprintf('\n')
    end
    fprintf('\n')
end

%% PLOT RESULTS
if nargin >= 4
    figure
    hold on
    Y = -0.5*ones(2,nv);
    Z = -ones(2,v);
    plot(NoVoids',Y,'b','LineWidth',1.5,'DisplayName','Data');
    plot(Voids',Z,'r','LineWidth',1.5,'DisplayName','Voids');
    datetick('x','dd/mm/yyyy','keepticks')
    title('Data Assessment')
    legend('show')
    legend('location','SouthWest'); legend('boxoff')
    set(gca,'YTickLabel',[],'Ylim',[-2 0.5])
    box on
end