function [NewEvent_mm] = iMHEA_Depure(Event_Date,Event_mm)
%iMHEA Depuration of repetitive tips above maximum intensity.
% [NewEvent_Date,NewEvent_mm] =
% iMHEA_AggregationDepure(Event_Date,Event_mm).
%
% Input:
% Event_Date    = dd/mm/yyyy hh:mm:ss [date format].
% Event_mm      = Precipitation tips [mm].
%
% Output:
% NewEvent_mm   = Precipitation tips, repetitions replaced by zero [mm].
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in August, 2017
% Last edited in November, 2017


%% INITIALISE VARIABLES
fprintf('\n')
fprintf('DEPURATION OF REPETITIVE RAINFALL TIPS ABOVE MAXIMUM INTENSITY.\n')
nd = 86400; % Number of seconds per day
% Minimum tip interval to merge events (slightly greater than 1 second).
MinT = 1.1/nd;

%% IDENTIFY AND REMOVE REPETITIVE EVENTS
% Calculate the time between tips.
Diff_Event_Date = diff(Event_Date);
Diff_Event_Date = [MinT*100;Diff_Event_Date];
% Identify tips separated by less than the minimum time MinT.
EventDiff = Diff_Event_Date <= MinT;
NewEvent_mm = Event_mm;
NewEvent_mm(EventDiff) = 0;

%% PRINT RESULTS
fprintf('Removing tips occurring faster than MinT = %6.2f seconds.\n',MinT*86400)
fprintf('Number of tips identified: %4i.\n',length(find(EventDiff)))
fprintf('Rainfall volume before depuration: %8.2f mm.\n',nansum(Event_mm))
fprintf('Rainfall volume after depuration: %8.2f mm.\n',nansum(NewEvent_mm))
fprintf('\n')