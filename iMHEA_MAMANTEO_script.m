%% SCRIPT TO PROCESS THE MAMANTEO DATA
%
% This script follows a workflow that process the hydrometeorological and
% tracer data in Huamantanga and Lima. It carries on a drought analysis on
% the long-term data from Peru's SENAMHI's and ANA's stations, followed by
% the processing of high-resolution time series from the iMHEA network. It
% estimates residence times in the mamanteo system and potential diversion
% volumes, both locally and regionally, and resulting modified hydrographs.
%
% Boris Ochoa Tocachi
% Imperial College London
% Created in November, 2017
% Last edited in March, 2019

close all
clc
disp('CUSTOM SCRIPT TO PROCESS MAMANTEO DATA')

disp(' ')
%% 0. Input data and user-defined parameters
disp('0. Loading iMHEA HMT (Huamantanga) and SENAMHI data')

load iMHEA_MAMANTEO_data

disp('Assimilating user-defined input information')
disp(' ')

% Define hydrological years
iMHEA_HMT_hy = [datetime('Sep 01, 2014') datetime('Aug 31, 2015');...
                datetime('Sep 01, 2015') datetime('Aug 31, 2016')];
disp('The considered hydrological years for iMHEA_HMT are:')
disp(iMHEA_HMT_hy)

% Define wet seasons
iMHEA_HMT_ws = [datetime('Dec 01, 2014') datetime('Apr 30, 2015');...
                datetime('Dec 01, 2015') datetime('Apr 30, 2016')];
disp('The considered wet seasons for iMHEA_HMT are:')
disp(iMHEA_HMT_ws)

% Define potential mamanteo diversion flow capacities (be aware of units)
% Minimum flows left in the stream (0: take all, 4: long-term daily mean)
userinput_QMin = [0 4]; % Units: l/s/km2
% Estimated intake maxima (inf: unlimited, 75: best estimate)
userinput_QMax = [inf 300 200 100 75]; % Units: l/s

% Define an estimated recovery rate (between 0 and 1)
userinput_frr = 0.50;

% Define a river flow diversion factor for the regional analysis as a
% precentage of wet season discharge in the highlands (between 0 and 1)
userinput_fdv = 0.50;

% Define a river flow reliability factor (between 0 and 1)
userinput_frl = 1;

% Display input data
disp('The considered system functioning parameters are:')
disp('Minimum flows left in the stream [l/s/km2]')
disp(userinput_QMin)
disp('Maximum system intake capacity [l/s]')
disp(userinput_QMax)
disp('Estimated recovery rate [fraction]')
disp(userinput_frr)
disp('Regional flow diversion factor [fraction]')
disp(userinput_fdv)
disp('Flow reliability factor [fraction]')
disp(userinput_frl)

% Define dates for drought and regional analysis
userinput_ax1 = [datetime([1989, 9, 1]) datetime([2018, 8, 31])];
userinput_ax2 = [datetime([1964, 9 ,1]) datetime([2018, 8, 31])];
userinput_ax3 = [datetime([1920, 9, 1]) datetime([1960, 8, 31])];
userinput_ax4 = [datetime([1960, 9, 1]) datetime([2018, 8, 31])];
disp('The considered periods for the regional, long-term analysis are:')
disp('Period when all statios have data (after SENAMHI_CHO)')
disp(userinput_ax1)
disp('Long-term Huamantanga rainfall data period (SENAMHI_HMT)')
disp(userinput_ax2)
disp('Natural flow regime period of the Rimac river (ANA_CHO)')
disp(userinput_ax3)
disp('Altered flow regime period of the Rimac river (ANA_CHO)')
disp(userinput_ax4)

disp('Done')
disp(' ')
%% 1. Process long-term hydrometeorological data from SENAMHI and ANA
disp('1. Processing long-term hydrometeorological data')

% Matrix SENAMHI_XXX_PP_1day
% column 1: datetime [number format]
% column 2: daily rainfall volumes [mm]

% Matrix ANA_XXX_HQ_1day
% column 1: datetime [number format]
% column 2: daily streamflow [m3/s]
% column 3: daily discharge volumes [million m3 /day]
% column 4: daily discharge volumes [mm]
% column 5: daily baseflow volumes [mm]

% Extract daily data
% Convert date vector to date time data
SENAMHI_CHO_PP_DDATE = SENAMHI_CHO_PP_table_1day{:,1};
SENAMHI_HMT_PP_DDATE = SENAMHI_HMT_PP_table_1day{:,1};
ANA_CHO_HQ_DDATE     = ANA_CHO_HQ_table_1day{:,1};
% Convert dates to numbers
SENAMHI_CHO_PP_1day(:,1) = datenum(SENAMHI_CHO_PP_DDATE);
SENAMHI_HMT_PP_1day(:,1) = datenum(SENAMHI_HMT_PP_DDATE);
ANA_CHO_HQ_1day(:,1)     = datenum(ANA_CHO_HQ_DDATE);
% Extract rainfall [mm] and streamflow [m3/s] data
SENAMHI_CHO_PP_1day(:,2) = SENAMHI_CHO_PP_table_1day{:,2};
SENAMHI_HMT_PP_1day(:,2) = SENAMHI_HMT_PP_table_1day{:,2};
ANA_CHO_HQ_1day(:,2)     = ANA_CHO_HQ_table_1day{:,2};
% Fill streamflow gaps using linear interpolation if needed
ANA_CHO_HQ_1day(:,2) = interp1(ANA_CHO_HQ_1day(~isnan(ANA_CHO_HQ_1day(:,2)),1),ANA_CHO_HQ_1day(~isnan(ANA_CHO_HQ_1day(:,2)),2),ANA_CHO_HQ_1day(:,1));
% Calculate discharge volumes [million m3 /day]
ANA_CHO_HQ_1day(:,3)     = ANA_CHO_HQ_1day(:,2)*86400/1e6;
% Volumes in mm
ANA_CHO_HQ_1day(:,4)     = ANA_CHO_HQ_1day(:,3)*1e9/(ANA_CHO_AREA_station*1e6);
% Calculate baseflow
[~,ANA_CHO_HQ_1day(:,5)] = iMHEA_BaseFlowUK(ANA_CHO_HQ_DDATE,ANA_CHO_HQ_1day(:,4));

% Matrix SENAMHI_XXX_PP_1yr or ANA_XXX_HQ_1day
% column 1: year
% column 2 P: annual rainfall volumes [mm]
% column 2 Q: annual discharge volumes [million m3 /day]

% Matrix SENAMHI_XXX_PP_Mmon or ANA_XXX_HQ_Mmon
% columns 1-12 P: monthly rainfall volumes, 1:Jan to 12:Dec [mm]
% columns 1-12 Q: monthly discharge volumes, 1:Jan to 12:Dec [million m3 /day]

% Create monthly and yearly data
[~,SENAMHI_CHO_PP_1yr,~,~,SENAMHI_CHO_PP_Mmon] = iMHEA_MonthlyRain(SENAMHI_CHO_PP_DDATE,SENAMHI_CHO_PP_1day(:,2));
[~,SENAMHI_HMT_PP_1yr,~,~,SENAMHI_HMT_PP_Mmon] = iMHEA_MonthlyRain(SENAMHI_HMT_PP_DDATE,SENAMHI_HMT_PP_1day(:,2));
[~,ANA_CHO_HQ_1yr,~,~,ANA_CHO_HQ_Mmon] = iMHEA_MonthlyRain(ANA_CHO_HQ_DDATE,ANA_CHO_HQ_1day(:,3));

% Matrix SENAMHI_XXX_PP_XXXXYYYY_Hmon or ANA_XXX_HQ_XXXXYYYY_Hmon
% columns 1-12 P: monthly rainfall using XXXX-YYYY hydrological years [mm]
% columns 1-12 Q: monthly discharge using XXXX-YYYY hydrological years [mm]

% Matrix SENAMHI_XXX_PP_XXXXYYYY_Hmon or ANA_XXX_HQ_XXXXYYYY_Hmon
% column 1: hydrological year
% column 2 P: annual rainfall using hydrological years [mm]
% column 2 Q: annual discharge using hydrological years [mm]

% Extract data from hydrological years
[ANA_CHO_HQ_19902018_Hmon,ANA_CHO_HQ_19902018_Hyr] = iMHEA_HydroYear(ANA_CHO_HQ_Mmon,ANA_CHO_HQ_1yr(:,1),userinput_ax1);
[SENAMHI_HMT_PP_19902018_Hmon,SENAMHI_HMT_PP_19902018_Hyr] = iMHEA_HydroYear(SENAMHI_HMT_PP_Mmon,SENAMHI_HMT_PP_1yr(:,1),userinput_ax1);
[SENAMHI_CHO_PP_19902018_Hmon,SENAMHI_CHO_PP_19902018_Hyr] = iMHEA_HydroYear(SENAMHI_CHO_PP_Mmon,SENAMHI_CHO_PP_1yr(:,1),userinput_ax1);

% Water demand 
FT_RIMAC_HQ_dem = FT_RIMAC_HQ_table_1mon{9,3:14};

disp('Done')
disp(' ')
%% 2. Drought Analysis
disp('2. Drought analysis of SENAMHIs HMT and CHO data')

drought_method = 3;
% method =  number from 1 to 4 specifying the method as follows:
%           1: Thr_DMA = Moving average of daily quantile (D_MA)
%           2: Thr_MMA = Moving average of monthly quantile (M_MA)
%           3: Thr_D30 = 30-day moving window quantile (30D)
%           4: Thr_FFT = Fast Fourier transform approach(D_FF)

drought_aggr_P = 1;
drought_aggr_Q = 1;
% aggreg = number 0 or 1 specyfing further method details as follows:
%           1: smooth series with 30-day moving average 
%           0: use original series (to evaluate Q only)

% Precipitation [mm] drought indices and threshold
[SENAMHI_CHO_PP_DIndices_19902018,SENAMHI_CHO_PP_DThr_19902018] = iMHEA_Drought(SENAMHI_CHO_PP_DDATE(and(SENAMHI_CHO_PP_DDATE>=userinput_ax1(1),SENAMHI_CHO_PP_DDATE<=userinput_ax1(2))),SENAMHI_CHO_PP_1day(and(SENAMHI_CHO_PP_DDATE>=userinput_ax1(1),SENAMHI_CHO_PP_DDATE<=userinput_ax1(2)),2),drought_method,drought_aggr_P,1);
[SENAMHI_HMT_PP_DIndices_19902018,SENAMHI_HMT_PP_DThr_19902018] = iMHEA_Drought(SENAMHI_HMT_PP_DDATE(and(SENAMHI_HMT_PP_DDATE>=userinput_ax1(1),SENAMHI_HMT_PP_DDATE<=userinput_ax1(2))),SENAMHI_HMT_PP_1day(and(SENAMHI_HMT_PP_DDATE>=userinput_ax1(1),SENAMHI_HMT_PP_DDATE<=userinput_ax1(2)),2),drought_method,drought_aggr_P,1);
[SENAMHI_HMT_PP_DIndices_19652018,SENAMHI_HMT_PP_DThr_19652018] = iMHEA_Drought(SENAMHI_HMT_PP_DDATE(and(SENAMHI_HMT_PP_DDATE>=userinput_ax2(1),SENAMHI_HMT_PP_DDATE<=userinput_ax2(2))),SENAMHI_HMT_PP_1day(and(SENAMHI_HMT_PP_DDATE>=userinput_ax2(1),SENAMHI_HMT_PP_DDATE<=userinput_ax2(2)),2),drought_method,drought_aggr_P,1);
% Streamflow [mm] drought indices and threshold
[ANA_CHO_HQ_DIndices_19902018,ANA_CHO_HQ_DThr_19902018] = iMHEA_Drought(ANA_CHO_HQ_DDATE(and(ANA_CHO_HQ_DDATE>=userinput_ax1(1),ANA_CHO_HQ_DDATE<=userinput_ax1(2))),ANA_CHO_HQ_1day(and(ANA_CHO_HQ_DDATE>=userinput_ax1(1),ANA_CHO_HQ_DDATE<=userinput_ax1(2)),4),drought_method,drought_aggr_Q,1);
[ANA_CHO_HQ_DIndices_19211960,ANA_CHO_HQ_DThr_19211960] = iMHEA_Drought(ANA_CHO_HQ_DDATE(and(ANA_CHO_HQ_DDATE>=userinput_ax3(1),ANA_CHO_HQ_DDATE<=userinput_ax3(2))),ANA_CHO_HQ_1day(and(ANA_CHO_HQ_DDATE>=userinput_ax3(1),ANA_CHO_HQ_DDATE<=userinput_ax3(2)),4),drought_method,drought_aggr_Q,1);

% Estimate monthly drought threshold (omit leap day)
ax18 = (datetime([2018, 1, 1]):datetime([2018, 12, 31]));
[SENAMHI_CHO_PP_DThrmon_19902018,SENAMHI_CHO_PP_DThryr_19902018] = iMHEA_MonthlyRain(ax18,SENAMHI_CHO_PP_DThr_19902018([1:59,61:366]));
[SENAMHI_HMT_PP_DThrmon_19902018,SENAMHI_HMT_PP_DThryr_19902018] = iMHEA_MonthlyRain(ax18,SENAMHI_HMT_PP_DThr_19902018([1:59,61:366]));
[ANA_CHO_HQ_DThrmon_19902018,ANA_CHO_HQ_DThryr_19902018] = iMHEA_MonthlyRain(ax18,ANA_CHO_HQ_DThr_19902018([1:59,61:366]));
[ANA_CHO_HQ_DThrmon_19211960,ANA_CHO_HQ_DThryr_19211960] = iMHEA_MonthlyRain(ax18,ANA_CHO_HQ_DThr_19211960([1:59,61:366]));

% Clear auxiliar variables
clear ax18

disp('Done')
disp(' ')
%% 3. Figure: Long-term historic records
disp('3. Creating Figure of long-term historic records in HMT and CHO')

% Initialise variables
offset_loc = [-0.2 0.1 0.02]; % To locate plot texts and scale axis
refmon_hy = month(userinput_ax1(1));
refyrs_hy = year(userinput_ax1);

fig_S1 = figure;
set(fig_S1,'Renderer','painters')

% Monthly precipitation in the highlands
disp('Ploting monthly precipitation in the highlands')
sp1 = subplot(2,2,1);
hold on
bar(mean(SENAMHI_HMT_PP_19902018_Hmon),'facecolor',0.8*[1 1 1],'DisplayName','Highlands'' precipitation (1990-2018)')
boxplot(SENAMHI_HMT_PP_19902018_Hmon,'PlotStyle','compact','colors',0*[1 1 1],'Labels',MonthLabels)
plot(SENAMHI_HMT_PP_DThrmon_19902018([refmon_hy:12,1:refmon_hy-1]),'-r','Linewidth',1,'DisplayName','Highlands'' precipitation drought threshold')
% Plot settings
set(gca,'Xlim',[0.5 12.5],...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*300,...
    'XTick',1:2:12,...
    'XTickLabel',...
    MonthLabels(1:2:end),...
    'YGrid','on',...
    'TickDir','out') % ,'FontSize',16)
box off
% Labels
text(offset_loc(1),1+offset_loc(2),'a','Units','normalized','FontWeight','bold','FontSize',16)
ylabel('mm month^{-1}')
% Legend
legend('location','NorthWest')
legend('boxoff')
legend('off')

% Monthly precipitation in Lima
disp('Ploting monthly precipitation in Lima')
sp2 = subplot(2,2,2);
hold on
bL = bar(mean(SENAMHI_CHO_PP_19902018_Hmon),'facecolor',0.3*[1 1 1],'DisplayName','Lowlands'' precipitation (1990-2018)');
boxplot(SENAMHI_CHO_PP_19902018_Hmon,'PlotStyle','compact','colors',0*[1 1 1],'Labels',MonthLabels)
dL   = plot(SENAMHI_CHO_PP_DThrmon_19902018([refmon_hy:12,1:refmon_hy-1]),'--r','Linewidth',1,'DisplayName','Lowlands'' precipitation drought threshold');
%-- For legend
bH = bar(-12:-1,mean(SENAMHI_HMT_PP_19902018_Hmon),'facecolor',0.8*[1 1 1],'DisplayName','Highlands'' precipitation (1990-2018)');
dH = plot(-12:-1,SENAMHI_HMT_PP_DThrmon_19902018([refmon_hy:12,1:refmon_hy-1]),'-r','Linewidth',1,'DisplayName','Highlands'' precipitation drought threshold');
%--
% Plot settings
set(gca,'Xlim',[0.5 12.5],...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*300,...
    'XTick',1:2:12,...
    'XTickLabel',...
    MonthLabels(1:2:end),...
    'YGrid','on',...
    'TickDir','out') % ,'FontSize',16)
box off
% Labels
text(offset_loc(1),1+offset_loc(2),'b','Units','normalized','FontWeight','bold','FontSize',16)
ylabel('mm month^{-1}')
% Legend
legend('location','NorthWest')
legend([bH,dH,bL,dL])
legend('boxoff')
clear bH dH bL dL

% Supply - Demand plot for discharge in Lima
disp('Supply - Demand plot for discharge in Lima')
sp3 = subplot(2,2,3);
hold on
% Variability ranges
patch([1:12,12:-1:1],[prctile(userinput_frl*ANA_CHO_HQ_19902018_Hmon,5),prctile(userinput_frl*ANA_CHO_HQ_19902018_Hmon(:,end:-1:1),95)],0.5*[1 1 1],'Edgecolor','none','FaceAlpha',.1,'DisplayName','95% uncertainty bounds')
% Mean volumes
m = plot(userinput_frl*mean(ANA_CHO_HQ_19902018_Hmon),'Linewidth',1,'color',0.5*[1 1 1],'DisplayName','Mean Rimac river flow (1990-2018)');
d = plot(FT_RIMAC_HQ_dem,'--k','Linewidth',1,'DisplayName','Total Lima''s water demand');
% plot(userinput_frl*ANA_CHO_HQ_DThrmon_19211960([refmon:12,1:refmon-1])*ANA_CHO_AREA_station/1e3,'-.m','Linewidth',1.5,'DisplayName','Natural Rimac river flow drought threshold')
D = plot(userinput_frl*ANA_CHO_HQ_DThrmon_19902018([refmon_hy:12,1:refmon_hy-1])*ANA_CHO_AREA_station/1e3,':r','Linewidth',2,'DisplayName','Altered Rimac river flow drought threshold');
y1 = plot(userinput_frl*ANA_CHO_HQ_19902018_Hmon(1990-refyrs_hy(1),:),'color',[0.85 0.33 0.10],'Linewidth',1,'DisplayName','1989-1990 Rimac river flow');
y2 = plot(userinput_frl*ANA_CHO_HQ_19902018_Hmon(1991-refyrs_hy(1),:),'color',[0.93 0.69 0.13],'Linewidth',1,'DisplayName','1990-1991 Rimac river flow');
y3 = plot(userinput_frl*ANA_CHO_HQ_19902018_Hmon(1992-refyrs_hy(1),:),'color',[0.64 0.08 0.18],'Linewidth',1,'DisplayName','1991-1992 Rimac river flow');
% Plot settings
set(gca,'Xlim',[0.5 12.5],...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*300,...
    'XTick',1:2:12,...
    'XTickLabel',...
    MonthLabels(1:2:end),...
    'YGrid','on',...
    'TickDir','out') % ,'FontSize',16)
box off
% Labels
text(offset_loc(1),1+offset_loc(2),'c','Units','normalized','FontWeight','bold','FontSize',16)
ylabel('million m^{3} month^{-1}')
% Legend
legend('location','NorthWest')
legend([m,d,D,y1,y2,y3])
legend('boxoff')
clear m d D y1 y2 y3

% Annual precipitation in Lima and in the highlands (1990 - 2013)
disp('Annual precipitation in Lima and in the highlands')
sp4 = subplot(2,2,4);
hold on
b = bar(SENAMHI_HMT_PP_19902018_Hyr(:,1),...
    [SENAMHI_HMT_PP_19902018_Hyr(:,2),...
    SENAMHI_CHO_PP_19902018_Hyr(:,2)],...
    'facecolor',0.8*[1 1 1],'barwidth',1);
b(2).FaceColor = 0.3*[1 1 1];
b(1).DisplayName = 'Highlands'' annual precipitation';
b(2).DisplayName = 'Lowlands'' annual precipitation';
% Thresholds
plot([refyrs_hy(1) refyrs_hy(2)+1],[1 1]*SENAMHI_HMT_PP_DThryr_19902018(2),'-r','Linewidth',1,'DisplayName','Highlands'' drought threshold')
plot([refyrs_hy(1) refyrs_hy(2)+1],[1 1]*SENAMHI_CHO_PP_DThryr_19902018(2),'--r','Linewidth',1,'DisplayName','Lowlands'' drought threshold')
% plot([refyrs_hy(1) refyrs_hy(2)+1],[1 1]*userinput_frl*ANA_CHO_HQ_DThryr_19211960),'-.m','Linewidth',1.5,'DisplayName','Natural Rimac river flow drought threshold')
plot([refyrs_hy(1) refyrs_hy(2)+1],[1 1]*userinput_frl*ANA_CHO_HQ_DThryr_19902018(2),':r','Linewidth',2,'DisplayName','Altered Rimac river flow drought threshold')
% Plot settings
set(gca,'Xlim',[refyrs_hy(1) refyrs_hy(2)]+0.5,...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*800,...
    'Xtick',(refyrs_hy(1)+1:4:refyrs_hy(2)),...
    'YGrid','on',...
    'TickDir','out') % ,'FontSize',16)
box off
% Labels
text(offset_loc(1),1+offset_loc(2),'d','Units','normalized','FontWeight','bold','FontSize',16)
ylabel('mm year^{-1}')
% Legend
legend('location','NorthWest')
legend('boxoff')
legend('off')
clear b

% Fix location
sp1.Position(1) = sp3.Position(1);
sp2.Position(1) = sp4.Position(1);
sp2.Position(2) = sp1.Position(2);
sp4.Position(2) = sp3.Position(2);
% Fix width
sp1.Position(3) = sp3.Position(3);
sp2.Position(3) = sp3.Position(3);
sp4.Position(3) = sp3.Position(3);
% Fix height
sp1.Position(4) = sp3.Position(4);
sp2.Position(4) = sp3.Position(4);
sp4.Position(4) = sp3.Position(4);

% Clear
clear sp1 sp2 sp3 sp4

% Export historic data Figure
disp('Exporting the Figure to pdf')
% Set the paper size [width height]
set(fig_S1,'PaperSize',[30 30]);
% set(fig_S1,'PaperOrientation','landscape')
print(fig_S1,'FigS1_export','-dpdf','-fillpage')

% Plot rainfall timeseries for Figure 1
fig_1 = figure;
set(fig_1,'Renderer','painters')

% Monthly precipitation in the highlands
disp('Ploting monthly precipitation in the highlands')
subplot(2,1,1);
hold on
bar(mean(SENAMHI_HMT_PP_19902018_Hmon),'facecolor',0.8*[1 1 1],'DisplayName','Highlands'' precipitation (1990-2018)')
boxplot(SENAMHI_HMT_PP_19902018_Hmon,'PlotStyle','compact','colors',0*[1 1 1],'Labels',MonthLabels)
plot(SENAMHI_HMT_PP_DThrmon_19902018([refmon_hy:12,1:refmon_hy-1]),'-r','Linewidth',1,'DisplayName','Highlands'' precipitation drought threshold')
% Plot settings
set(gca,'Xlim',[0.5 12.5],...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*300,...
    'XTick',1:2:12,...
    'XTickLabel',...
    MonthLabels(1:2:end),...
    'YGrid','on',...
    'TickDir','out') % ,'FontSize',16)
box off
% Labels
text(offset_loc(1)/1.5,1+offset_loc(2),'a','Units','normalized','FontWeight','bold','FontSize',16)
ylabel('mm month^{-1}')
% Legend
legend('location','NorthWest')
legend('boxoff')
legend('off')

% Monthly precipitation in Lima
disp('Ploting monthly precipitation in Lima')
subplot(2,1,2);
hold on
bL = bar(mean(SENAMHI_CHO_PP_19902018_Hmon),'facecolor',0.3*[1 1 1],'DisplayName','Lowlands'' precipitation (1990-2018)');
boxplot(SENAMHI_CHO_PP_19902018_Hmon,'PlotStyle','compact','colors',0*[1 1 1],'Labels',MonthLabels)
dL   = plot(SENAMHI_CHO_PP_DThrmon_19902018([refmon_hy:12,1:refmon_hy-1]),'--r','Linewidth',1,'DisplayName','Lowlands'' precipitation drought threshold');
%-- For legend
bH = bar(-12:-1,mean(SENAMHI_HMT_PP_19902018_Hmon),'facecolor',0.8*[1 1 1],'DisplayName','Highlands'' precipitation (1990-2018)');
dH = plot(-12:-1,SENAMHI_HMT_PP_DThrmon_19902018([refmon_hy:12,1:refmon_hy-1]),'-r','Linewidth',1,'DisplayName','Highlands'' precipitation drought threshold');
%--
% Plot settings
set(gca,'Xlim',[0.5 12.5],...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*300,...
    'XTick',1:2:12,...
    'XTickLabel',...
    MonthLabels(1:2:end),...
    'YGrid','on',...
    'TickDir','out') % ,'FontSize',16)
box off
% Labels
text(offset_loc(1)/1.5,1+offset_loc(2),'b','Units','normalized','FontWeight','bold','FontSize',16)
ylabel('mm month^{-1}')
% Legend
legend('location','NorthWest')
legend([bH,dH,bL,dL])
legend('boxoff')
clear bH dH bL dL

% Export Figure
disp('Exporting the Figure to pdf')
% Set the paper size [width height]
set(fig_1,'PaperSize',[15 30]);
% set(fig_S1,'PaperOrientation','landscape')
print(fig_1,'Fig1_export','-dpdf','-fillpage')

% Clear auxiliar variables
clear fig_S1 fig_1 offset_loc refmon_hy refyrs_hy

disp('Done')
disp(' ')
%% 4. Process high-resolution micro-catchment data from iMHEA
disp('4. Processing high-resolution micro-catchment data')

% Workflow to preprocess raw data and generate high resolution time series
disp('Processing catchment C1')
[iMHEA_HMT_01_DataHRes] = iMHEA_WorkflowMamanteo(iMHEA_HMT_01_AREA,iMHEA_HMT_01_HI_01_raw{:,1},iMHEA_HMT_01_HI_01_raw{:,3},0.2,iMHEA_HMT_01_PO_01_raw{:,1},iMHEA_HMT_01_PO_01_raw{:,2},iMHEA_HMT_01_PO_02_raw{:,1},iMHEA_HMT_01_PO_02_raw{:,2});
disp('Processing catchment C2')
[iMHEA_HMT_02_DataHRes] = iMHEA_WorkflowMamanteo(iMHEA_HMT_02_AREA,iMHEA_HMT_02_HI_01_raw{:,1},iMHEA_HMT_02_HI_01_raw{:,3},0.2,iMHEA_HMT_02_PO_01_raw{:,1},iMHEA_HMT_02_PO_01_raw{:,2},iMHEA_HMT_02_PO_02_raw{:,1},iMHEA_HMT_02_PO_02_raw{:,2});

disp('Done')
disp(' ')
%% 5. Consolidate paired catchments and fill gaps
disp('5. Consolidating paired catchment data')

% Merge both catchments in a single paired catchment dataset
[iMHEA_HMT_Pair_DataHRes] = iMHEA_WorkflowPairMamanteo(iMHEA_HMT_01_DataHRes(:,1:3),iMHEA_HMT_02_DataHRes(:,1:3));

% Fill precipitation gaps using paired catchment correlation if needed
disp('Filling precipitation gaps')
iMHEA_HMT_Pair_DataHRes(:,8:9) = iMHEA_HMT_Pair_DataHRes(:,4:5);
iMHEA_HMT_Pair_DataHRes(:,4:5) = iMHEA_HMT_Pair_DataHRes(:,2:3);
[~,iMHEA_HMT_Pair_DataHRes(:,2),~,~] = iMHEA_FillGaps(iMHEA_HMT_01_DataHRes(:,1),iMHEA_HMT_01_DataHRes(:,4),iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_Pair_DataHRes(:,4));
[~,iMHEA_HMT_Pair_DataHRes(:,3),~,~] = iMHEA_FillGaps(iMHEA_HMT_01_DataHRes(:,1),iMHEA_HMT_01_DataHRes(:,5),iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_Pair_DataHRes(:,4));
[~,iMHEA_HMT_Pair_DataHRes(:,6),~,~] = iMHEA_FillGaps(iMHEA_HMT_02_DataHRes(:,1),iMHEA_HMT_02_DataHRes(:,4),iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_Pair_DataHRes(:,8));
[~,iMHEA_HMT_Pair_DataHRes(:,7),~,~] = iMHEA_FillGaps(iMHEA_HMT_02_DataHRes(:,1),iMHEA_HMT_02_DataHRes(:,5),iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_Pair_DataHRes(:,8));

% Fill streamflow gaps using linear interpolation if needed
disp('Filling streamflow gaps')
iMHEA_HMT_Pair_DataHRes(:,5) = interp1(iMHEA_HMT_Pair_DataHRes(~isnan(iMHEA_HMT_Pair_DataHRes(:,5)),1),iMHEA_HMT_Pair_DataHRes(~isnan(iMHEA_HMT_Pair_DataHRes(:,5)),5),iMHEA_HMT_Pair_DataHRes(:,1));
iMHEA_HMT_Pair_DataHRes(:,9) = interp1(iMHEA_HMT_Pair_DataHRes(~isnan(iMHEA_HMT_Pair_DataHRes(:,9)),1),iMHEA_HMT_Pair_DataHRes(~isnan(iMHEA_HMT_Pair_DataHRes(:,9)),9),iMHEA_HMT_Pair_DataHRes(:,1));

% Extract vector with datetime values in numbers
iMHEA_HMT_Pair_datevec = datetime(iMHEA_HMT_Pair_DataHRes(:,1),'ConvertFrom','datenum');

% Plot resulting datasets
iMHEA_Plot3(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,2),iMHEA_HMT_Pair_DataHRes(:,3),iMHEA_HMT_Pair_DataHRes(:,4),iMHEA_HMT_Pair_DataHRes(:,5));
iMHEA_Plot3(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,6),iMHEA_HMT_Pair_DataHRes(:,7),iMHEA_HMT_Pair_DataHRes(:,8),iMHEA_HMT_Pair_DataHRes(:,9));

disp('Done')
disp(' ')
%% 6. Process data and calculate hydrological indices
disp('6. Processing data and calculating hydrological indices')

% Use rainfall aggregation and streamflow average at daily scale (1440 min)
disp('Aggregate and average data')
% Processing C1
iMHEA_HMT_01_PO_1day_mm = iMHEA_Aggregation(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_Pair_DataHRes(:,4),1440);
iMHEA_HMT_01_PO_1day_mm(:,3:4) = [];
iMHEA_HMT_01_HQ_1day_mm = iMHEA_Average(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_Pair_DataHRes(:,5),1440);
iMHEA_HMT_01_HQ_1day_mm(:,3:4) = [];
[iMHEA_HMT_01_HQ_1day_mm(:,2)] = iMHEA_HMT_01_HQ_1day_mm(:,2)*86400/1000000;
% Processing C2
iMHEA_HMT_02_PO_1day_mm = iMHEA_Aggregation(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_Pair_DataHRes(:,8),1440);
iMHEA_HMT_02_PO_1day_mm(:,3:4) = [];
iMHEA_HMT_02_HQ_1day_mm = iMHEA_Average(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_Pair_DataHRes(:,9),1440);
iMHEA_HMT_02_HQ_1day_mm(:,3:4) = [];
[iMHEA_HMT_02_HQ_1day_mm(:,2)] = iMHEA_HMT_02_HQ_1day_mm(:,2)*86400/1000000;

disp('Calculating baseflows in the streams')
[~,iMHEA_HMT_01_HQ_1day_mm(:,3)] = iMHEA_BaseFlowUK(iMHEA_HMT_01_HQ_1day_mm(:,1),iMHEA_HMT_01_HQ_1day_mm(:,2));
[~,iMHEA_HMT_02_HQ_1day_mm(:,3)] = iMHEA_BaseFlowUK(iMHEA_HMT_02_HQ_1day_mm(:,1),iMHEA_HMT_02_HQ_1day_mm(:,2));
% Interpolate baseflow at 5-min resolution
iMHEA_HMT_Pair_DataHRes(:,10) = interp1(iMHEA_HMT_01_HQ_1day_mm(:,1),iMHEA_HMT_01_HQ_1day_mm(:,3)/86400*1000000,iMHEA_HMT_Pair_DataHRes(:,1),'linear','extrap');
iMHEA_HMT_Pair_DataHRes(:,11) = interp1(iMHEA_HMT_02_HQ_1day_mm(:,1),iMHEA_HMT_02_HQ_1day_mm(:,3)/86400*1000000,iMHEA_HMT_Pair_DataHRes(:,1),'linear','extrap');

disp('Estimating runoff ratios from hydrological years')
% Processing C1
iMHEA_HMT_01_PHYear = nansum(iMHEA_HMT_01_PO_1day_mm(and(iMHEA_HMT_01_PO_1day_mm(:,1)>=datenum(iMHEA_HMT_hy(1)),iMHEA_HMT_01_PO_1day_mm(:,1)<datenum(iMHEA_HMT_hy(end))),2))/2;
iMHEA_HMT_01_PMaxi  = nanmax(iMHEA_HMT_01_PO_1day_mm(and(iMHEA_HMT_01_PO_1day_mm(:,1)>=datenum(iMHEA_HMT_hy(1)),iMHEA_HMT_01_PO_1day_mm(:,1)<datenum(iMHEA_HMT_hy(end))),2));
iMHEA_HMT_01_QHYear = nansum(iMHEA_HMT_01_HQ_1day_mm(and(iMHEA_HMT_01_HQ_1day_mm(:,1)>=datenum(iMHEA_HMT_hy(1)),iMHEA_HMT_01_HQ_1day_mm(:,1)<datenum(iMHEA_HMT_hy(end))),2))/2;
iMHEA_HMT_01_QMean = nanmean(iMHEA_HMT_01_HQ_1day_mm(and(iMHEA_HMT_01_HQ_1day_mm(:,1)>=datenum(iMHEA_HMT_hy(1)),iMHEA_HMT_01_HQ_1day_mm(:,1)<datenum(iMHEA_HMT_hy(end))),2))/86400*1000000;

% Processing C2
iMHEA_HMT_02_PHYear = nansum(iMHEA_HMT_02_PO_1day_mm(and(iMHEA_HMT_02_PO_1day_mm(:,1)>=datenum(iMHEA_HMT_hy(1)),iMHEA_HMT_02_PO_1day_mm(:,1)<datenum(iMHEA_HMT_hy(end))),2))/2;
iMHEA_HMT_02_PMaxi  = nanmax(iMHEA_HMT_02_PO_1day_mm(and(iMHEA_HMT_02_PO_1day_mm(:,1)>=datenum(iMHEA_HMT_hy(1)),iMHEA_HMT_02_PO_1day_mm(:,1)<datenum(iMHEA_HMT_hy(end))),2));
iMHEA_HMT_02_QHYear = nansum(iMHEA_HMT_02_HQ_1day_mm(and(iMHEA_HMT_02_HQ_1day_mm(:,1)>=datenum(iMHEA_HMT_hy(1)),iMHEA_HMT_02_HQ_1day_mm(:,1)<datenum(iMHEA_HMT_hy(end))),2))/2;
iMHEA_HMT_02_QMean = nanmean(iMHEA_HMT_02_HQ_1day_mm(and(iMHEA_HMT_02_HQ_1day_mm(:,1)>=datenum(iMHEA_HMT_hy(1)),iMHEA_HMT_02_HQ_1day_mm(:,1)<datenum(iMHEA_HMT_hy(end))),2))/86400*1000000;

% 1st estimation: runoff ratio = annual discharge (HY) / annual precipitation (HY)
iMHEA_HMT_01_RR(:,1) = iMHEA_HMT_01_QHYear/iMHEA_HMT_01_PHYear;
iMHEA_HMT_02_RR(:,1) = iMHEA_HMT_02_QHYear/iMHEA_HMT_02_PHYear;

% Calculate average precipitation and discharge in Huamantanga for reference
iMHEA_HMT_PHYear = (iMHEA_HMT_01_PHYear+iMHEA_HMT_02_PHYear)/2;
iMHEA_HMT_QHYear = (iMHEA_HMT_01_QHYear+iMHEA_HMT_02_QHYear)/2;
iMHEA_HMT_RRHYear = iMHEA_HMT_QHYear/iMHEA_HMT_PHYear;

% Calculate average dicharge from satellite precipitation climatology
ANA_CHO_TRMM_4000{:,9} = ANA_CHO_TRMM_4000{:,4}*(iMHEA_HMT_02_RR(1)+iMHEA_HMT_01_RR(1))/2;
ANA_CHO_TRMM_4000{:,10} = ANA_CHO_TRMM_4000{:,9}/1000*ANA_CHO_TRMM_4000{1,4}*10000;

disp('Estimating runoff ratios from daily data during the monitoring period')
% 2nd estimation: runoff ratio = long-term discharge / long-term precipitation
iMHEA_HMT_01_RR(:,2) = nanmean(iMHEA_HMT_01_HQ_1day_mm(:,2))/nanmean(iMHEA_HMT_01_PO_1day_mm(:,2));
iMHEA_HMT_02_RR(:,2) = nanmean(iMHEA_HMT_02_HQ_1day_mm(:,2))/nanmean(iMHEA_HMT_02_PO_1day_mm(:,2));

disp('Processing data in chronological years')
% Matrix iMHEA_HMT_0x_XCYear
% column 1:     chronological year
% column 2 P:   rainfall at yearly scale [mm]
% column 2 Q:   streamflow at yearly scale[mm]

% Annual rainfall
[~,iMHEA_HMT_01_PCYear] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,4));
[~,iMHEA_HMT_02_PCYear] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,8));

% Annual discharge
[~,iMHEA_HMT_01_QCYear] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,4)*5*60/1000000);
[~,iMHEA_HMT_02_QCYear] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,9)*5*60/1000000);

% 3rd estimation: runoff ratio = annual discharge (CY) / annual precipitation (CY)
iMHEA_HMT_01_RRYear = iMHEA_HMT_01_QCYear./iMHEA_HMT_01_PCYear;
iMHEA_HMT_02_RRYear = iMHEA_HMT_02_QCYear./iMHEA_HMT_02_PCYear;
iMHEA_HMT_01_RRYear(:,1) = iMHEA_HMT_01_QCYear(:,1);
iMHEA_HMT_02_RRYear(:,1) = iMHEA_HMT_02_QCYear(:,1);

disp('Estimating runoff ratios from average site data')
% 4th estimation: yearly runoff ratio (CY)
iMHEA_HMT_PYear = (iMHEA_HMT_01_PCYear+iMHEA_HMT_02_PCYear)/2;
iMHEA_HMT_QYear = (iMHEA_HMT_01_QCYear+iMHEA_HMT_02_QCYear)/2;
% Calculate runoff ratios
iMHEA_HMT_RR = iMHEA_HMT_QYear./iMHEA_HMT_PYear;
iMHEA_HMT_RR(:,1) = iMHEA_HMT_QYear(:,1);

disp('Processing data in hydrological years')
% Initialise the variable
iMHEA_HMT_1mon = zeros(48,6);

% Monthly rainfall
[iMHEA_HMT_1mon(:,1)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,2));
[iMHEA_HMT_1mon(:,2)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,3));
[iMHEA_HMT_1mon(:,3)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,6));
[iMHEA_HMT_1mon(:,4)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,7));

% Monthly streamflow
[iMHEA_HMT_1mon(:,5)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,5)*5*60/1000000);
[iMHEA_HMT_1mon(:,6)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,9)*5*60/1000000);

% Rainfall and streamflow in hydrological years
% Matrix iMHEA_HMT_1yr:
% column 1: hydrological year
% column 2: HMT_01_PO_01
% column 3: HMT_01_PO_02
% column 4: HMT_02_PO_01
% column 5: HMT_02_PO_02
% column 6: HMT_01_HQ_01
% column 7: HMT_02_HQ_01

% Initiliase the matrix
iMHEA_HMT_1yr = [20142105 ; 20152016];
% Allocate data
refmon_hy = month(iMHEA_HMT_hy(1));
for i =1:6
    iMHEA_HMT_1yr(1,i+1) = sum(iMHEA_HMT_1mon(refmon_hy:refmon_hy+12,i));
    iMHEA_HMT_1yr(2,i+1) = sum(iMHEA_HMT_1mon(refmon_hy+13:refmon_hy+24,i));
end

% Clear auxiliar variables
clear i refmon_hy

disp('Done')
disp(' ')
%% 7. Diverted water volumes in Huamantanga
disp('7. Estimating potential diverted water volumes in Huamantanga')

% Matrix iMHEA_HMT_0x_MAMFLOW_1yr:
% column 1-2:  wet season year
% column 3-end: combination QMax & QMin
[iMHEA_HMT_01_MAMFLOW_HRes,iMHEA_HMT_01_MAMFLOW_1yr] = iMHEA_Diversion(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,5),userinput_QMin,userinput_QMax/iMHEA_HMT_01_AREA,iMHEA_HMT_ws);
[iMHEA_HMT_02_MAMFLOW_HRes,iMHEA_HMT_02_MAMFLOW_1yr] = iMHEA_Diversion(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,9),userinput_QMin,userinput_QMax/iMHEA_HMT_02_AREA,iMHEA_HMT_ws);

% Processing C1
[iMHEA_HMT_01_MAMFLOW_1day] = iMHEA_Aggregation(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_01_MAMFLOW_HRes(:,1)*5*60/1000000,1440);
iMHEA_HMT_01_MAMFLOW_1day(:,3:end) = [];
% Processing C2
[iMHEA_HMT_02_MAMFLOW_1day] = iMHEA_Aggregation(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_02_MAMFLOW_HRes(:,1)*5*60/1000000,1440);
iMHEA_HMT_02_MAMFLOW_1day(:,3:end) = [];

disp('Extracting minima and maxima water balance and wet season volumes')
% Matrix iMHEA_HMT_MAMVOL_max/min:
% row index:    catchment number
% column 1:     annual rainfall (hydrological year)
% column 2:     annual streamflow (hydrological year)
% column 3:     potential diversion (umlimited wet season)
% column 4:     effective diversion (wet season & limits)

% Initialise the matrices
iMHEA_HMT_MAMVOL_min = zeros(2,4);
iMHEA_HMT_MAMVOL_max = zeros(2,4);
nM = length(userinput_QMax);
nm = length(userinput_QMin);

% Minima
% Precipitation
iMHEA_HMT_MAMVOL_min(1,1) = min(min(iMHEA_HMT_1yr(:,2:3)));
iMHEA_HMT_MAMVOL_min(2,1) = min(min(iMHEA_HMT_1yr(:,4:5)));
% Streamflow
iMHEA_HMT_MAMVOL_min(1,2) = min(iMHEA_HMT_1yr(:,6));
iMHEA_HMT_MAMVOL_min(2,2) = min(iMHEA_HMT_1yr(:,7));
% Potential diversion
iMHEA_HMT_MAMVOL_min(1,3) = min(iMHEA_HMT_01_MAMFLOW_1yr(:,3));
iMHEA_HMT_MAMVOL_min(2,3) = min(iMHEA_HMT_02_MAMFLOW_1yr(:,3));
% Effective diversion
iMHEA_HMT_MAMVOL_min(1,4) = min(iMHEA_HMT_01_MAMFLOW_1yr(:,nM*nm+2));
iMHEA_HMT_MAMVOL_min(2,4) = min(iMHEA_HMT_02_MAMFLOW_1yr(:,nM*nm+2));

% Maxima
% Precipitation
iMHEA_HMT_MAMVOL_max(1,1) = max(max(iMHEA_HMT_1yr(:,2:3)));
iMHEA_HMT_MAMVOL_max(2,1) = max(max(iMHEA_HMT_1yr(:,4:5)));
% Streamflow
iMHEA_HMT_MAMVOL_max(1,2) = max(iMHEA_HMT_1yr(:,6));
iMHEA_HMT_MAMVOL_max(2,2) = max(iMHEA_HMT_1yr(:,7));
% Potential diversion
iMHEA_HMT_MAMVOL_max(1,3) = max(iMHEA_HMT_01_MAMFLOW_1yr(:,3));
iMHEA_HMT_MAMVOL_max(2,3) = max(iMHEA_HMT_02_MAMFLOW_1yr(:,3));
% Effective diversion
iMHEA_HMT_MAMVOL_max(1,4) = max(iMHEA_HMT_01_MAMFLOW_1yr(:,nM*nm+2));
iMHEA_HMT_MAMVOL_max(2,4) = max(iMHEA_HMT_02_MAMFLOW_1yr(:,nM*nm+2));

% Clear auxiliar variables
clear nM nm

disp('Done')
disp(' ')
%% 8. Estimate residence times from the average tracer concentration
disp('8. Estimating residence times from the average tracer concentration')

disp('Using tracer data and baseflow at daily scale')
% Extract values from the tracer data table
iMHEA_HMT_Tracers_Avg = [cumsum(iMHEA_HMT_Tracers_Pmean{:,4}),iMHEA_HMT_Tracers_Pmean{:,6}];
% Fill any nan data using linear interpolation
iMHEA_HMT_Tracers_Avg(:,2) = interp1(iMHEA_HMT_Tracers_Avg(~isnan(iMHEA_HMT_Tracers_Avg(:,2)),1),iMHEA_HMT_Tracers_Avg(~isnan(iMHEA_HMT_Tracers_Avg(:,2)),2),iMHEA_HMT_Tracers_Avg(:,1));
% Extend the tracer data to deplet at the end of the year
iMHEA_HMT_Tracers_Avg(end+1,1) = 300;
iMHEA_HMT_Tracers_Avg(end,2) = 0;
% Calculate cumulative tracer concentration
iMHEA_HMT_Tracers_Avg(:,3) = cumsum(iMHEA_HMT_Tracers_Avg(:,2));

% Matrix iMHEA_HMT_Tracers_Int_1day:
% column 1:     day index
% column 2:     average tracer concentration interpolated at daily scale
% column 3:     cumulative tracer concentration at daily scale
% column 4:     daily residence time distribution (scaled by spring flow)
% column 5:     daily cumulative residence time distribution
% column 6:     daily spring flow (derived from catchment baseflow)

% Interpolate tracer concentration at 1-day resolution
iMHEA_HMT_Tracers_Int_1day = (iMHEA_HMT_Tracers_Avg(1,1):1:iMHEA_HMT_Tracers_Avg(end,1))';
iMHEA_HMT_Tracers_Int_1day(:,2) = interp1(iMHEA_HMT_Tracers_Avg(:,1),iMHEA_HMT_Tracers_Avg(:,2),iMHEA_HMT_Tracers_Int_1day(:,1),'linear','extrap');
iMHEA_HMT_Tracers_Int_1day(iMHEA_HMT_Tracers_Int_1day(:,2)<0,2) = 0;
iMHEA_HMT_Tracers_Int_1day(:,3) = cumsum(iMHEA_HMT_Tracers_Int_1day(:,2));
iMHEA_HMT_Tracers_Int_1day(:,1) = datenum(iMHEA_HMT_Tracers_Pmean{1,2})+iMHEA_HMT_Tracers_Int_1day(:,1);
% Ponderate relative to the baseflow (proxy for spring flow)
iMHEA_HMT_Tracers_Int_1day(:,6) = iMHEA_HMT_01_HQ_1day_mm(and(iMHEA_HMT_01_HQ_1day_mm(:,1)>=iMHEA_HMT_Tracers_Int_1day(1,1),iMHEA_HMT_01_HQ_1day_mm(:,1)<=iMHEA_HMT_Tracers_Int_1day(end,1)),3);
iMHEA_HMT_Tracers_Int_1day(:,4) = iMHEA_HMT_Tracers_Int_1day(:,2).*iMHEA_HMT_Tracers_Int_1day(:,6);
% Scale the residence times and accumulate them
iMHEA_HMT_Tracers_Int_1day(:,4) = iMHEA_HMT_Tracers_Int_1day(:,4)/sum(iMHEA_HMT_Tracers_Int_1day(:,4));
iMHEA_HMT_Tracers_Int_1day(:,5) = cumsum(iMHEA_HMT_Tracers_Int_1day(:,4));

disp('Interpolating residence times at 5-min scale')
% Matrix iMHEA_HMT_Tracers_Int_5min:
% column 1:     day index (fraction of a day)
% column 2:     average tracer concentration interpolated at 5-min scale
% column 3:     cumulative tracer concentration at 5-min scale
% column 4:     residence time distribution (scaled by spring flow) at 5-min
% column 5:     cumulative residence time distribution at 5-min
% column 6:     spring flow (derived from catchment baseflow) at 5-min

% Interpolate average tracer concentration at 5-min resolution
iMHEA_HMT_Tracers_Int_5min = (iMHEA_HMT_Tracers_Int_1day(1,1):1/(60/5*24):iMHEA_HMT_Tracers_Int_1day(end,1))';
iMHEA_HMT_Tracers_Int_5min(:,2) = interp1(iMHEA_HMT_Tracers_Int_1day(:,1),iMHEA_HMT_Tracers_Int_1day(:,2),iMHEA_HMT_Tracers_Int_5min(:,1),'linear','extrap');
iMHEA_HMT_Tracers_Int_5min(:,3) = cumsum(iMHEA_HMT_Tracers_Int_5min(:,2));
% Interpolate spring flow at 5-min resolution
iMHEA_HMT_Tracers_Int_5min(:,6) = interp1(iMHEA_HMT_Tracers_Int_1day(:,1),iMHEA_HMT_Tracers_Int_1day(:,6),iMHEA_HMT_Tracers_Int_5min(:,1),'linear','extrap');
% Ponderate relative to the baseflow (proxy for spring flow)
iMHEA_HMT_Tracers_Int_5min(:,4) = iMHEA_HMT_Tracers_Int_5min(:,2).*iMHEA_HMT_Tracers_Int_5min(:,6);
% Scale the residence times and accumulate them
iMHEA_HMT_Tracers_Int_5min(:,4) = iMHEA_HMT_Tracers_Int_5min(:,4)/sum(iMHEA_HMT_Tracers_Int_5min(:,4));
iMHEA_HMT_Tracers_Int_5min(:,5) = cumsum(iMHEA_HMT_Tracers_Int_5min(:,4));

disp('Rescale the tracer data')
% Scale the residence times and accumulate them
iMHEA_HMT_Tracers_Avg(:,1) = datenum(iMHEA_HMT_Tracers_Pmean{1,2})+iMHEA_HMT_Tracers_Avg(:,1);
iMHEA_HMT_Tracers_Avg(:,4) = iMHEA_HMT_Tracers_Avg(:,2)/sum(iMHEA_HMT_Tracers_Avg(:,2));
iMHEA_HMT_Tracers_Avg(:,5) = cumsum(iMHEA_HMT_Tracers_Avg(:,4));

disp('Calculating average residence time')
% Matrix iMHEA_HMT_Tracers_MTT:
% column 1: residence time [days]
% row 1:    mean residence time
% row 2:    residence time standard deviation
% row 3:    median residence time

iMHEA_HMT_Tracers_MTT = sum(iMHEA_HMT_Tracers_Int_1day(:,4).*(iMHEA_HMT_Tracers_Int_1day(:,1)-iMHEA_HMT_Tracers_Avg(1,1)));
iMHEA_HMT_Tracers_MTT(2) = sqrt(sum(iMHEA_HMT_Tracers_Int_1day(:,4).*((iMHEA_HMT_Tracers_Int_1day(:,1)-iMHEA_HMT_Tracers_Avg(1,1)-iMHEA_HMT_Tracers_MTT(1)).^2)));
iMHEA_HMT_Tracers_MTT(3) = interp1(iMHEA_HMT_Tracers_Int_1day(20:100,5),iMHEA_HMT_Tracers_Int_1day(20:100,1),0.50,'linear','extrap') - iMHEA_HMT_Tracers_Avg(1,1);

disp('Done')
disp(' ')
%% 9. Modified hydrograph after mamanteo
disp('9. Estimating the modified hydrograph affected by the practice')

% Residence times at 5-min scale
l5mbf = iMHEA_HMT_Tracers_Int_5min(:,4);

% Apply diverted flows, residence times, and recovery rate
[iMHEA_HMT_01_MAMFLOW_HRes(:,3)] = iMHEA_ModifiedFlows(iMHEA_HMT_01_MAMFLOW_HRes(:,2),iMHEA_HMT_01_MAMFLOW_HRes(:,1),l5mbf,userinput_frr);
[iMHEA_HMT_02_MAMFLOW_HRes(:,3)] = iMHEA_ModifiedFlows(iMHEA_HMT_02_MAMFLOW_HRes(:,2),iMHEA_HMT_02_MAMFLOW_HRes(:,1),l5mbf,userinput_frr);
% [iMHEA_HMT_01_MAMFLOW_HRes(:,3)] = iMHEA_ModifiedFlowsBF(iMHEA_HMT_01_MAMFLOW_HRes(:,2),iMHEA_HMT_01_MAMFLOW_HRes(:,1),iMHEA_HMT_Pair_DataHRes(:,10),l5mbf,userinput_frr);
% [iMHEA_HMT_02_MAMFLOW_HRes(:,3)] = iMHEA_ModifiedFlowsBF(iMHEA_HMT_02_MAMFLOW_HRes(:,2),iMHEA_HMT_02_MAMFLOW_HRes(:,1),iMHEA_HMT_Pair_DataHRes(:,11),l5mbf,userinput_frr);

disp('Calculating cumulative water volumes at daily scale')
% Matrix iMHEA_HMT_0x_MAMFLOW_1day
% column 1: date
% column 2: diverted volumes at daily scale [mm]
% column 3: remaning volumes at daily scale[mm]
% column 4: modified hydrograph at daily scale [mm]

% Processing C1
[~,iMHEA_HMT_01_MAMFLOW_1day(:,3)] = iMHEA_Aggregation(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_01_MAMFLOW_HRes(:,2)*5*60/1000000,1440);
[~,iMHEA_HMT_01_MAMFLOW_1day(:,4)] = iMHEA_Aggregation(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_01_MAMFLOW_HRes(:,3)*5*60/1000000,1440);

% Processing C2
[~,iMHEA_HMT_02_MAMFLOW_1day(:,3)] = iMHEA_Aggregation(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_02_MAMFLOW_HRes(:,2)*5*60/1000000,1440);
[~,iMHEA_HMT_02_MAMFLOW_1day(:,4)] = iMHEA_Aggregation(iMHEA_HMT_Pair_DataHRes(:,1),iMHEA_HMT_02_MAMFLOW_HRes(:,3)*5*60/1000000,1440);

disp('Calculating cumulative water volumes at monthly scale')
% Matrix iMHEA_HMT_0x_MAMFLOW_1mon
% column index: month index from Jan-Year1
% column 1:     rainfall at monthly scale [mm]
% column 2:     streamflow at monthly scale[mm]
% column 3:     modified hydrograph at monthly scale[mm]

% Initialise the variable
iMHEA_HMT_01_MAMFLOW_1mon = zeros(48,3);
iMHEA_HMT_02_MAMFLOW_1mon = zeros(48,3);
refmon_hy = month(iMHEA_HMT_hy(1));

% Annual rainfall
[iMHEA_HMT_01_MAMFLOW_1mon(:,1)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,4));
[iMHEA_HMT_02_MAMFLOW_1mon(:,1)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,8));

% Annual discharge
[iMHEA_HMT_01_MAMFLOW_1mon(:,2)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,5)*5*60/1000000);
[iMHEA_HMT_02_MAMFLOW_1mon(:,2)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,9)*5*60/1000000);

% Monthly modified hydrograph [mm]
[iMHEA_HMT_01_MAMFLOW_1mon(:,3)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_01_MAMFLOW_HRes(:,3)*5*60/1000000);
[iMHEA_HMT_02_MAMFLOW_1mon(:,3)] = iMHEA_MonthlyRain(iMHEA_HMT_Pair_datevec,iMHEA_HMT_02_MAMFLOW_HRes(:,3)*5*60/1000000);

% Cut matrices at hydrological years
iMHEA_HMT_01_MAMFLOW_Hmon = iMHEA_HMT_01_MAMFLOW_1mon;
iMHEA_HMT_02_MAMFLOW_Hmon = iMHEA_HMT_02_MAMFLOW_1mon;
iMHEA_HMT_01_MAMFLOW_Hmon([1:refmon_hy-1,refmon_hy+24:end],:) = [];
iMHEA_HMT_02_MAMFLOW_Hmon([1:refmon_hy-1,refmon_hy+24:end],:) = [];

% Aggregate interannual monthly data in a single matrix
iMHEA_HMT_HQ_Hmon = [iMHEA_HMT_01_MAMFLOW_Hmon(1:12,2),iMHEA_HMT_01_MAMFLOW_Hmon(13:24,2),iMHEA_HMT_02_MAMFLOW_Hmon(1:12,2),iMHEA_HMT_02_MAMFLOW_Hmon(13:24,2)]';
iMHEA_HMT_MM_Hmon = [iMHEA_HMT_01_MAMFLOW_Hmon(1:12,3),iMHEA_HMT_01_MAMFLOW_Hmon(13:24,3),iMHEA_HMT_02_MAMFLOW_Hmon(1:12,3),iMHEA_HMT_02_MAMFLOW_Hmon(13:24,3)]';

% Variability ranges
iMHEA_HMT_HQ_Hmonvar = [1:12,12:-1:1];
iMHEA_HMT_HQ_Hmonvar(2,:) = [prctile(iMHEA_HMT_HQ_Hmon,5),prctile(iMHEA_HMT_HQ_Hmon(:,end:-1:1),95)];
iMHEA_HMT_MM_Hmonvar = [1:12,12:-1:1];
iMHEA_HMT_MM_Hmonvar(2,:) = [prctile(iMHEA_HMT_MM_Hmon,5),prctile(iMHEA_HMT_MM_Hmon(:,end:-1:1),95)];

% Summary table of recovered flows and increased dry season flow percentage
% Matrix iMHEA_HMT_MM:
% column 1:     Recovery rate (-)
% column 2:     Recovered volume (m3/yr)
% column 3:     Long-term monthly average dry season increase (%)
% column 4:     Long-term monthly max dry season increase (%)
% column 5:     Maximum monthly mean dry season increase (%)
% column 6:     Average monthly mean dry season increase (%)
% column 7:     Minimum monthly mean dry season increase (%)

iMHEA_HMT_MM = userinput_frr;
auxvol = iMHEA_HMT_HQ_Hmon(:) - iMHEA_HMT_MM_Hmon(:);
iMHEA_HMT_MM(2) = nansum(auxvol)/size(iMHEA_HMT_HQ_Hmon,1)*userinput_frr/(1-userinput_frr);
auxdiff = (iMHEA_HMT_MM_Hmon - iMHEA_HMT_HQ_Hmon)./iMHEA_HMT_HQ_Hmon;
iMHEA_HMT_MM(3) = mean(auxdiff(auxdiff>0))*100;
iMHEA_HMT_MM(4) = max(auxdiff(:))*100;
auxmonth = mean(auxdiff);
iMHEA_HMT_MM(5) = min(auxmonth(auxmonth>0))*100;
iMHEA_HMT_MM(6) = mean(auxmonth(auxmonth>0))*100;
iMHEA_HMT_MM(7) = max(auxmonth)*100;

% Clear auxiliar variables
clear refmon_hy auxdiff auxmonth auxvol

disp('Done')
disp(' ')
%% 10. Plot estimated residence times
disp('10. Ploting estimated residence times')

% Initialise variables
offset_loc = [-0.2 0.1 0.05]; % To locate plot texts and scale axis
aux_loct = interp1(iMHEA_HMT_Tracers_Int_1day(20:100,1),iMHEA_HMT_Tracers_Int_1day(20:100,5),iMHEA_HMT_Tracers_MTT(1)+iMHEA_HMT_Tracers_Avg(1,1),'linear','extrap');

fig_S4 = figure;
set(fig_S4,'Renderer','painters')

% subplot(4,2,2)
subplot(5,1,1)
hold on
plot(datetime(iMHEA_HMT_Tracers_Int_1day(:,1),'ConvertFrom','datenum'),...
    iMHEA_HMT_Tracers_Int_1day(:,6),'-','LineWidth',1,...
    'color',[0.08 0.17 0.55],'DisplayName','Estimated spring flow (baseflow)')
% Plot settings
set(gca,'YGrid','on','TickDir','out','Ylim',[0-offset_loc(3) 1+offset_loc(3)]*1)
auxx = get(gca,'Xlim');
datetick('x','mmm-yy','keeplimits','keepticks')
box off
% Labels
ylabel('mm day^{-1}')
text(offset_loc(1),1+offset_loc(2),'b','Units','normalized','FontWeight','bold','FontSize',16)
% Legend
legend('Location','NorthEast')
legend('boxoff')

% subplot(4,2,4)
subplot(5,1,2)
hold on
% stairs(datetime([iMHEA_HMT_Tracers_Avg(1,1);iMHEA_HMT_Tracers_Avg(:,1)],'ConvertFrom','datenum'),...
%     [iMHEA_HMT_Tracers_Avg(:,2);0],...
%     ':','LineWidth',1,'color',[0.49 0.18 0.56])
a = plot(datetime(iMHEA_HMT_Tracers_Avg(:,1),'ConvertFrom','datenum'),...
    iMHEA_HMT_Tracers_Avg(:,2),...
    'o','LineWidth',1,'color',[0.49 0.18 0.56],...
    'DisplayName','Original tracer data');
b = plot(datetime(iMHEA_HMT_Tracers_Int_1day(:,1),'ConvertFrom','datenum'),...
    iMHEA_HMT_Tracers_Int_1day(:,2),...
    '-','LineWidth',1,'color',[0.85 0.33 0.01],...
    'DisplayName','Interpolation at daily scale');
% Plot settings
set(gca,'YGrid','on','TickDir','out','Ylim',[0-offset_loc(3) 1+offset_loc(3)]*10)
ax17 = get(gca,'Xlim');
datetick('x','mmm-yy','keeplimits','keepticks')
box off
% Labels
ylabel('proxy concentration [ppb]')
text(offset_loc(1),1+offset_loc(2),'c','Units','normalized','FontWeight','bold','FontSize',16)
% Legend
legend('Location','NorthEast')
legend('boxoff')
legend([a,b])
clear a b

% subplot(4,2,6)
subplot(5,1,3)
plot(datetime(iMHEA_HMT_Tracers_Int_1day(:,1),'ConvertFrom','datenum'),...
    iMHEA_HMT_Tracers_Int_1day(:,4),...
    '--','LineWidth',1,'color',[0 0.5 0],...
    'DisplayName','Fractional residence time scaled by spring flow')
% Plot settings
set(gca,'YGrid','on','TickDir','out','Ylim',[0-offset_loc(3) 1+offset_loc(3)]*0.03)
datetick('x','mmm-yy','keeplimits','keepticks')
box off
% Labels
ylabel('fraction [-]')
xlabel('Date')
text(offset_loc(1),1+offset_loc(2),'d','Units','normalized','FontWeight','bold','FontSize',16)
% Legend
legend('Location','NorthEast')
legend('boxoff')

% subplot(4,2,[1,3,5,7])
subplot(5,1,4)
hold on
% plot(iMHEA_HMT_Tracers_Avg(:,1)-iMHEA_HMT_Tracers_Avg(1,1),...
%     iMHEA_HMT_Tracers_Avg(:,5),...
%     'o','LineWidth',1,'color',[0.49 0.18 0.56],...
%     'DisplayName','Original tracer data')
plot(iMHEA_HMT_Tracers_Int_1day(:,1)-iMHEA_HMT_Tracers_Avg(1,1),...
    iMHEA_HMT_Tracers_Int_1day(:,3)/sum(iMHEA_HMT_Tracers_Int_1day(:,2)),...
    '-','LineWidth',1,'color',[0.85 0.33 0.01],...
    'DisplayName','Proxy concentration at daily scale')
plot(iMHEA_HMT_Tracers_Int_1day(:,1)-iMHEA_HMT_Tracers_Avg(1,1),...
    iMHEA_HMT_Tracers_Int_1day(:,5),...
    '--','LineWidth',1,'color',[0 0.5 0],...
    'DisplayName','Residence time scaled by spring flow')
plot(iMHEA_HMT_Tracers_MTT(1),aux_loct,'sr','Markersize',8,...
    'MarkerFaceColor','r','DisplayName','Mean residence time')
plot(iMHEA_HMT_Tracers_MTT(3),0.50,'sk','Markersize',8,...
    'MarkerFaceColor','k','DisplayName','Median residence time')
% Plot settings
set(gca,'YGrid','on','TickDir','out','Ylim',[0-offset_loc(3) 1+offset_loc(3)]*1,...
    'Xlim',datenum(auxx - iMHEA_HMT_Tracers_Avg(1,1)))
box off
% Labels
ylabel('cumulative distribution [-]')
xlabel('residence time [day]')
text(offset_loc(1),1+offset_loc(2)/4,'a','Units','normalized','FontWeight','bold','FontSize',16)
% Legend
legend('Location','SouthEast')
legend('boxoff')
legend('AutoUpdate','off');
% Auxiliar lines
plot(iMHEA_HMT_Tracers_MTT(3)*[1 1],[-.1 0.50],...
    ':k','LineWidth',1)
plot([datenum(auxx(1) - iMHEA_HMT_Tracers_Avg(1,1)) iMHEA_HMT_Tracers_MTT(3)],0.50*[1 1],...
    ':k','LineWidth',1)
plot(iMHEA_HMT_Tracers_MTT(1)*[1 1],[-.1 aux_loct],...
    ':r','LineWidth',1)
plot([datenum(auxx(1) - iMHEA_HMT_Tracers_Avg(1,1)) iMHEA_HMT_Tracers_MTT(1)],aux_loct*[1 1],...
    ':r','LineWidth',1)

% subplot(4,2,8)
subplot(5,1,5)
hold on
plot(datetime(iMHEA_HMT_01_HQ_1day_mm(:,1),'ConvertFrom','datenum'),...
    iMHEA_HMT_01_HQ_1day_mm(:,2),'Color',[0 0.45 0.74],...
    'LineWidth',1,'DisplayName','Discharge without diversion')
plot(datetime(iMHEA_HMT_01_MAMFLOW_1day(:,1),'ConvertFrom','datenum'),...
    iMHEA_HMT_01_MAMFLOW_1day(:,4),'Color',[0 0.5 0],...
    'LineWidth',1,'DisplayName','Discharge with diversion')
% Plot settings
set(gca,'YGrid','on','TickDir','out','Ylim',[0-offset_loc(3) 1+offset_loc(3)]*10,'Xlim',ax17)
datetick('x','mmm-yy','keeplimits','keepticks')
box off
% Labels
ylabel('mm day^{-1}')
xlabel('Date')
text(offset_loc(1),1+offset_loc(2),'e','Units','normalized','FontWeight','bold','FontSize',16)
% Legend
legend('Location','NorthEast')
legend('boxoff')

% Export Figure
disp('Exporting the Figure to pdf')
% Set the paper size [width height]
set(fig_S4,'PaperSize',[20 30]);
% set(fig_S4,'PaperOrientation','landscape')
print(fig_S4,'FigS4_export','-dpdf','-fillpage')

% Clear axuliar variables
clear fig_S4 offset_loc ax17 aux_loct auxx

disp('Done')
disp(' ')
%% 11. Plot results of delayed monthly hydrograph
disp('11. Ploting results of delayed monthly hydrograph')

% Initialise variables
refhyv = datevec(iMHEA_HMT_hy);
refhydt = datetime(refhyv(1,1),refhyv(1,2):12*(refhyv(end,1)-refhyv(1,1))+refhyv(end,2),1);
offset_loc = [-0.2 0.05 0.02]; % To locate plot texts and scale axis
lim1 = 200;

fig_S5 = figure;
set(fig_S5,'Renderer','painters')

% Catchment C1
disp('Ploting delayed hydrograph of catchment C1')
subplot(1,2,1)
hold on
b = bar (refhydt,iMHEA_HMT_01_MAMFLOW_Hmon(:,1),'grouped','barwidth',1);
% plot(refhydt,iMHEA_HMT_01_MAMFLOW_Hmon(:,1),'Linewidth',1,'Color',0.5*[1 1 1],'DisplayName','Catchment C1 rainfall');
plot(refhydt,iMHEA_HMT_01_MAMFLOW_Hmon(:,2),'Linewidth',2,'Color',[0 0.45 0.74],'DisplayName','Original C1 discharge');
plot(refhydt,iMHEA_HMT_01_MAMFLOW_Hmon(:,3),'Linewidth',2,'Color',[0 0.5 0],'DisplayName','Potential effect of infiltration system');
% Plot settings
b(1).FaceColor = 'none';
b(1).EdgeColor = 0.5*[1 1 1];
b(1).DisplayName = 'Catchment C1 rainfall';
set(gca,'Ylim',[-offset_loc(3) 1+offset_loc(3)]*lim1,...
    'Xlim',[refhydt(1)-30 refhydt(end)+30],...
    'YGrid','on',...
    'TickDir','out')
datetick('x','mmm-yy','keeplimits','keepticks')
box off
% Labels
ylabel('mm month^{-1}')
% title('Huamantanga catchment C1')
text(offset_loc(1),1+offset_loc(2),'a','Units','normalized','FontWeight','bold','FontSize',16)
% Legend
legend('AutoUpdate','off');
legend('location','Northwest')
legend('boxoff')

% Plot seasons
plot([iMHEA_HMT_hy(1) iMHEA_HMT_hy(end)],[lim1 lim1],':k','Linewidth',1)
for i = 1:size(iMHEA_HMT_ws,1)
    plot(iMHEA_HMT_ws(i,:),[lim1 lim1],'k','LineWidth',2)
end
text(mean(iMHEA_HMT_ws(1,:)),lim1*(1+offset_loc(2)),'Wet season','HorizontalAlignment','center')% ,'FontSize',14)

% Catchment C2
disp('Ploting delayed hydrograph of catchment C2')
subplot(1,2,2)
hold on
b = bar (refhydt,iMHEA_HMT_02_MAMFLOW_Hmon(:,1),'grouped','barwidth',1);
plot(refhydt,iMHEA_HMT_02_MAMFLOW_Hmon(:,2),'Linewidth',2,'Color',[0 0.45 0.74],'DisplayName','Original C2 discharge');
plot(refhydt,iMHEA_HMT_02_MAMFLOW_Hmon(:,3),'Linewidth',2,'Color',[0 0.5 0],'DisplayName','Potential effect of infiltration system');
% Plot settings
b(1).FaceColor = 'none';
b(1).EdgeColor = 0.5*[1 1 1];
b(1).DisplayName = 'Catchment C2 rainfall';
set(gca,'Ylim',[-offset_loc(3) 1+offset_loc(3)]*lim1,...
    'Xlim',[refhydt(1)-30 refhydt(end)+30],...
    'YGrid','on',...
    'TickDir','out')
datetick('x','mmm-yy','keeplimits','keepticks')
box off
% Labels
ylabel('mm month^{-1}')
% title('Huamantanga catchment C2')
text(offset_loc(1),1+offset_loc(2),'b','Units','normalized','FontWeight','bold','FontSize',16)
% Legend
legend('AutoUpdate','off');
legend('location','Northwest')
legend('boxoff')
% Plot seasons
plot([iMHEA_HMT_hy(1) iMHEA_HMT_hy(end)],[lim1 lim1],':k','Linewidth',1)
for i = 1:size(iMHEA_HMT_ws,1)
    plot(iMHEA_HMT_ws(i,:),[lim1 lim1],'k','LineWidth',2)
end
text(mean(iMHEA_HMT_ws(2:3)),lim1*(1+offset_loc(2)),'Dry season','HorizontalAlignment','center')% ,'FontSize',14)

% Export Figure
disp('Exporting the Figure to pdf')
% Set the paper size [width height]
set(fig_S5,'PaperSize',[30 15]);
% set(fig_S5,'PaperOrientation','landscape')
print(fig_S5,'FigS5_export','-dpdf','-fillpage')

% Clear auxiliar variables
clear fig_S5 refhyv refhydt offset_loc lim1 i b

disp('Done')
disp(' ')
%% 12. Figure: Rainfall, runoff, and tracer concentrations
disp('12. Creating Figure of hydrological monitoring and tracer data')

% Initialise variables
offset_loc = [-0.1 0.2 0.04];
% Extract data for the analysis periods
refhyv = datevec(iMHEA_HMT_hy);
refhydt = datetime(refhyv(1,1),refhyv(1,2):12*(refhyv(end,1)-refhyv(1,1))+refhyv(end,2),1);

fig_3 = figure;
set(fig_3,'Renderer','painters')

% Tracer concentrations
disp('Ploting tracer concentrations')
subplot(4,1,4)
% area([iMHEA_HMT_Tracers_Pmean{9,2};iMHEA_HMT_Tracers_Pmean{9,3}],[20 20],'facecolor',0.9*[1 1 1],'edgecolor',0.9*[1 1 1])
hold on
c = plot(iMHEA_HMT_Tracers_P1{:,2},...
    iMHEA_HMT_Tracers_P1{:,5},...
    '-o','DisplayName','Spring 1',...
    'Color',[0.64 0.08 0.18],'Linewidth',1);
d = plot(iMHEA_HMT_Tracers_P2{:,2},...
    iMHEA_HMT_Tracers_P2{:,5},...
    '--o','DisplayName','Spring 2',...
    'Color',[0.93 0.69 0.13],'Linewidth',1);
e = plot(iMHEA_HMT_Tracers_P3{:,2},...
    iMHEA_HMT_Tracers_P3{:,5},...
    '-.o','DisplayName','Spring 3',...
    'Color',[0.49 0.18 0.56],'Linewidth',1);
f = plot(iMHEA_HMT_Tracers_P4{:,2},...
    iMHEA_HMT_Tracers_P4{:,5},...
    ':o','DisplayName','Spring 4',...
    'Color',[0 0.5 0],'Linewidth',1);
% Plot settings
set(gca,'XTick',refhydt(1:2:end),'YLim',[-offset_loc(3) 1+offset_loc(3)]*20,...
    'Xlim',[refhydt(1)-3 refhydt(16)+3],...
    'TickDir','out') % ,'FontSize',16) % ,'xticklabelrotation',90)
datetick('x','mmm-yy','keepticks','keeplimits')
box off
% Labels
text(offset_loc(1),1+offset_loc(2),'c','Units','normalized','FontWeight','bold','FontSize',16)
% text(0.025,0.875,'    Tracer concentration','Units','normalized','FontSize',14)
ylabel('ppb')
xlabel('Date')
% Legend
legend('AutoUpdate','off')
legend('location','NorthWest')
legend('boxoff')
legend([c,d,e,f],'Spring 1 tracer concentration',...
    'Spring 2 tracer concentration',...
    'Spring 3 tracer concentration',...
    'Spring 4 tracer concentration')
% legend('off')
clear c d e f

% Plot seasons
disp('Ploting seasons')
plot([iMHEA_HMT_hy(1) iMHEA_HMT_hy(end)],[20 20],':k','Linewidth',1)
for i = 1:size(iMHEA_HMT_ws,1)
    plot(iMHEA_HMT_ws(i,:),[20 20],'k','LineWidth',2)
end
% text(mean(iMHEA_HMT_ws(1,:)),22,'Wet season','HorizontalAlignment','center')% ,'FontSize',14)
% text(mean(iMHEA_HMT_ws(2:3)),22,'Dry season','HorizontalAlignment','center')% ,'FontSize',14)

% Discharge
disp('Ploting discharge')
subplot(4,1,2:3)
hold on
a = area(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,5),...
    'FaceColor',[0.3 0.75 0.93],'EdgeColor','none');
plot(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,5),...
    'Color',[0.3 0.75 0.93],'Linewidth',1);
% Baseflow
f = plot(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,10),...
    'Color',[0.08 0.17 0.55],'Linewidth',1);
% Diversion limits
plot([iMHEA_HMT_hy(1)-30 iMHEA_HMT_hy(end)+30],userinput_QMin(end)*[1 1],'--r','Linewidth',1); % Min
plot([iMHEA_HMT_hy(1)-30 iMHEA_HMT_hy(end)+30],userinput_QMax(end)*[1 1]/iMHEA_HMT_01_AREA+userinput_QMin(end),'-r','Linewidth',1); % Max

% Plot settings
set(gca,'XTick',refhydt(1:2:end),'YLim',[-offset_loc(3)/2 1+offset_loc(3)/2]*100,...
    'Xlim',[refhydt(1)-3 refhydt(16)+3],...
    'TickDir','out') % ,'FontSize',16)
datetick('x','mmm-yy','keepticks','keeplimits')
% set(gca,'Xticklabel',[])
box off
% Labels
text(offset_loc(1),1+offset_loc(2)/2,'b','Units','normalized','FontWeight','bold','FontSize',16)
% text(0.025,0.95,'    Discharge','Units','normalized') % ,'FontSize',14)
text(mean(iMHEA_HMT_ws(2:3)),9,'Minimum flow left in the stream','HorizontalAlignment','center')% ,'FontSize',14)
text(mean(iMHEA_HMT_ws(2:3)),45,'Maximum diversion canal capacity','HorizontalAlignment','center')% ,'FontSize',14)
ylabel('l s^{-1} km^{-2}')
% Legend
legend('AutoUpdate','off');
legend('location','Northeast')
legend('boxoff')
legend([a,f],'Catchment outlet discharge','Baseflow')
% legend('off')

% Peak labels
disp('Ploting peak labels')
% Identify flows over 300 l/s/km2
peaks = find(iMHEA_HMT_Pair_DataHRes(:,5)>300);
% Consider only those peaks within the shown dates
limit = find(iMHEA_HMT_Pair_DataHRes(:,1)>=datenum(iMHEA_HMT_hy(1,1)),1);
peaks(peaks<=limit) = [];
limit = find(iMHEA_HMT_Pair_DataHRes(:,1)>=datenum(iMHEA_HMT_hy(1,2)),1);
peaks(peaks>=limit) = [];
peaks(:,2) = iMHEA_HMT_Pair_DataHRes(peaks(:,1),5);
% Initialise text
peaktext = 0;
pos1 = false;
pos2 = true;
% Iterate through them
for i = 1:length(peaks)
    % Calculate max peak for printing
    peaktext = max(peaktext,iMHEA_HMT_Pair_DataHRes(peaks(i),5));
    % Check if they correspond to the same peak
    if i==length(peaks) || peaks(i)+1~=peaks(i+1)
        % Print peak label
        text(iMHEA_HMT_Pair_datevec(peaks(i)),...
            107.5+3.5*pos1,...
            num2str(ceil(peaktext)),...
            'HorizontalAlignment','center') % ,'FontSize',14)
        peaktext = 0;
        % Iterate their position
        pos1 = ~pos1;
        if pos2 && pos1
            pos1 = -pos1;
            pos2 = ~pos2;
        elseif pos1
            pos2 = ~pos2;
        end
    end
end
clear peaks peaktext limit pos i

% %-- For legend
% b = plot([iMHEA_HMT_hy(1)-30 iMHEA_HMT_hy(end)+30],userinput_QMin(end)*[1 1],'--r','Linewidth',1); % Min
% c = plot([iMHEA_HMT_hy(1)-30 iMHEA_HMT_hy(end)+30],userinput_QMax(end)*[1 1]/iMHEA_HMT_01_AREA+userinput_QMin(end),'-r','Linewidth',1); % Max
% cS = plot(iMHEA_HMT_Tracers_P1{:,2}-9999,...
%     iMHEA_HMT_Tracers_P1{:,5},...
%     '-o','DisplayName','Spring 1',...
%     'Color',[0.64 0.08 0.18],'Linewidth',1);
% dS = plot(iMHEA_HMT_Tracers_P2{:,2}-9999,...
%     iMHEA_HMT_Tracers_P2{:,5},...
%     '--o','DisplayName','Spring 2',...
%     'Color',[0.93 0.69 0.13],'Linewidth',1);
% eS = plot(iMHEA_HMT_Tracers_P3{:,2}-9999,...
%     iMHEA_HMT_Tracers_P3{:,5},...
%     '-.o','DisplayName','Spring 3',...
%     'Color',[0.49 0.18 0.56],'Linewidth',1);
% fS = plot(iMHEA_HMT_Tracers_P4{:,2}-9999,...
%     iMHEA_HMT_Tracers_P4{:,5},...
%     ':o','DisplayName','Spring 4',...
%     'Color',[0 0.5 0],'Linewidth',1);
% R = area(datetime(iMHEA_HMT_Pair_DataHRes(:,1)-9999,'ConvertFrom','datenum'),iMHEA_HMT_Pair_DataHRes(:,4),...
%     'FaceColor',0.5*[1 1 1],'EdgeColor','none');
% legend('AutoUpdate','off');
% legend('location','Northeast')
% legend('boxoff')
% legend([R,a,f,b,c,cS,dS,eS,fS],'Average catchment rainfall',...
%     'Catchment outlet discharge','Baseflow',...
%     'Minimum flow left in the stream','Maximum diversion canal capacity',...
%     'Spring 1 tracer concentration','Spring 2 tracer concentration',...
%     'Spring 3 tracer concentration','Spring 4 tracer concentration')
% clear cS dS eS fS R
% %--
clear a f b c

% Rainfall
disp('Ploting rainfall')
subplot(4,1,1)
hold on
a = area(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,4),...
    'FaceColor',0.5*[1 1 1],'EdgeColor','none');
plot(iMHEA_HMT_Pair_datevec,iMHEA_HMT_Pair_DataHRes(:,4),...
    'Color',0.5*[1 1 1],'Linewidth',1);
% Plot settings
set(gca,'XTick',refhydt(1:2:end),'YLim',[-offset_loc(3) 1+offset_loc(3)]*4,...    
    'Xlim',[refhydt(1)-3 refhydt(16)+3],...
    'TickDir','out') % ,'FontSize',16)
datetick('x','mmm-yy','keepticks','keeplimits')
% set(gca,'Xticklabel',[])
box off
% Labels
text(offset_loc(1),1+offset_loc(2),'a','Units','normalized','FontWeight','bold' ,'FontSize',16)
% text(0.025,0.875,'    Rainfall','Units','normalized') % ,'FontSize',14)
ylabel('mm 5min^{-1}')
% Legend
legend('AutoUpdate','off');
legend('location','Northeast')
legend('boxoff')
legend(a,'Average catchment rainfall')
% legend('off')
clear a

% Plot seasons
disp('Ploting seasons')
plot([iMHEA_HMT_hy(1) iMHEA_HMT_hy(end)],[4 4],':k','Linewidth',1)
for i = 1:size(iMHEA_HMT_ws,1)
    plot(iMHEA_HMT_ws(i,:),[4 4],'k','LineWidth',2)
end
text(mean(iMHEA_HMT_ws(1,:)),4.5,'Wet season','HorizontalAlignment','center')% ,'FontSize',14)
text(mean(iMHEA_HMT_ws(2:3)),4.5,'Dry season','HorizontalAlignment','center')% ,'FontSize',14)

% Export Figure
disp('Exporting the Figure to pdf')
% Set the paper size [width height]
set(fig_3,'PaperSize',[25 25]);
% set(fig_3,'PaperOrientation','landscape')
print(fig_3,'Fig3_export','-dpdf','-fillpage')

% Clear auxiliar variables offset_loc
clear fig_3 offset_loc i pos1 pos2 refhydt refhyv

disp('Done')
disp(' ')
%% 13. Estimate regional effects for Rimac
disp('13. Estimating regional effects for Rimac')

% Residence times at daily scale
l1dbf = iMHEA_HMT_Tracers_Int_1day(:,4);

% Matrix ANA_CHO_HQ_19612018_date
% column 1: datetime [number format]

% Matrix ANA_CHO_HQ_XXXXYYYY_1day
% column 1: daily streamflow between years XXXX and YYYY [m3/s]
% column 2: daily discharge volumes between years XXXX and YYYY [mm]

% Matrix ANA_CHO_HQ_XXXXYYYY_Hmon
% columns 1-12: monthly discharge using XXXX-YYYY hydrological years [mm]

% Natural flow regime 1921 - 1960
refyrs_hy = year(userinput_ax3);
% Reorganise in hydrological years
ANA_CHO_HQ_19211960_Hmon = iMHEA_HydroYear(ANA_CHO_HQ_Mmon,ANA_CHO_HQ_1yr(:,1),userinput_ax3);

% Extract daily data
ANA_CHO_HQ_19211960_1day = ANA_CHO_HQ_1day(and(year(ANA_CHO_HQ_DDATE)>=refyrs_hy(1),year(ANA_CHO_HQ_DDATE)<=refyrs_hy(2)),3);
ANA_CHO_HQ_19211960_1day(:,2) = ANA_CHO_HQ_19211960_1day*1e9/(ANA_CHO_AREA_station*1e6);

% Altered flow regime with grey infrastructure 1961 - present
refyrs_hy = year(userinput_ax4);
% Reorganise in hydrological years
[ANA_CHO_HQ_19612018_Hmon,ANA_CHO_HQ_19612018_Hyr] = iMHEA_HydroYear(ANA_CHO_HQ_Mmon,ANA_CHO_HQ_1yr(:,1),userinput_ax4);
% Extract daily data
ANA_CHO_HQ_19612018_1day = ANA_CHO_HQ_1day(and(year(ANA_CHO_HQ_DDATE)>=refyrs_hy(1),year(ANA_CHO_HQ_DDATE)<=refyrs_hy(2)),3);
ANA_CHO_HQ_19612018_DDATE = ANA_CHO_HQ_DDATE(and(year(ANA_CHO_HQ_DDATE)>=refyrs_hy(1),year(ANA_CHO_HQ_DDATE)<=refyrs_hy(2)));
ANA_CHO_HQ_19612018_1day(:,2) = ANA_CHO_HQ_19612018_1day*1e9/(ANA_CHO_AREA_station*1e6);

% Matrix ANA_CHO_MAMFLOW_XXXXYYYY_1day
% column 1: wet season discharge at daily scale [m3/s]
% column 2: diverted volumes at daily scale [m3/s]
% column 3: remaning volumes at daily scale [m3/s]
% column 4: modified hydrograph at daily scale [m3/s]
% column 5: modified hydrograph at daily scale [mm]

% Potential diversion to mamanteo systems 1961 - present
refmon_ws = month(iMHEA_HMT_ws(1,:));
n = refyrs_hy(2)-refyrs_hy(1);
% Define wet seasons
ANA_CHO_ws_19602018 = [datetime([(refyrs_hy(1):refyrs_hy(2)-1)',refmon_ws(1)*ones(n,1),ones(n,1)]),datetime([(refyrs_hy(1)+1:refyrs_hy(2))',refmon_ws(2)*ones(n,1),30*ones(n,1)])];

% Extract flow during the wet seasons
ANA_CHO_MAMFLOW_19612018_1day = zeros(size(ANA_CHO_HQ_19612018_1day(:,1)));
for k = 1:n
    % Fill with data only between wet season limits
    ANA_CHO_MAMFLOW_19612018_1day(and(ANA_CHO_HQ_19612018_DDATE>ANA_CHO_ws_19602018(k,1),ANA_CHO_HQ_19612018_DDATE<ANA_CHO_ws_19602018(k,2))) = ANA_CHO_HQ_19612018_1day(and(ANA_CHO_HQ_19612018_DDATE>ANA_CHO_ws_19602018(k,1),ANA_CHO_HQ_19612018_DDATE<ANA_CHO_ws_19602018(k,2)));
end
% Calculate flow diversion percentage using precipitation climatology in the highlands
% divert_fdv = userinput_fdv*ANA_CHO_TRMM_4000{15,7}/(nansum(ANA_CHO_MAMFLOW_19612018_1day)/n);
divert_fdv = userinput_fdv*ANA_CHO_TRMM_4000{15,7}/ANA_CHO_TRMM_station{15,7};
% Calculate actual diverted flow
ANA_CHO_MAMFLOW_19612018_1day(:,2) = ANA_CHO_MAMFLOW_19612018_1day(:,1)*divert_fdv;

% Flows remaining in the stream
ANA_CHO_MAMFLOW_19612018_1day(:,3) = ANA_CHO_HQ_19612018_1day(:,1) - ANA_CHO_MAMFLOW_19612018_1day(:,2);
% Apply diverted flows, residence times, and recovery rate
[ANA_CHO_MAMFLOW_19612018_1day(:,4)] = iMHEA_ModifiedFlows(ANA_CHO_MAMFLOW_19612018_1day(:,3),ANA_CHO_MAMFLOW_19612018_1day(:,2),l1dbf,userinput_frr);
% [ANA_CHO_MAMFLOW_19612018_1day(:,4)] = iMHEA_ModifiedFlowsBF(ANA_CHO_MAMFLOW_19612018_1day(:,3),ANA_CHO_MAMFLOW_19612018_1day(:,2),ANA_CHO_HQ_1day(:,5),l1dbf,userinput_frr);
ANA_CHO_MAMFLOW_19612018_1day(:,5) = ANA_CHO_MAMFLOW_19612018_1day(:,4)*1e9/(ANA_CHO_AREA_station*1e6);

% Monthly modified hydrograph
[~,ANA_CHO_MM_19612018_1yr,~,~,ANA_CHO_MM_19612018_Mmon] = iMHEA_MonthlyRain(ANA_CHO_HQ_19612018_DDATE,ANA_CHO_MAMFLOW_19612018_1day(:,4));
% Reorganise in hydrological years
ANA_CHO_MM_19612018_Hmon = iMHEA_HydroYear(ANA_CHO_MM_19612018_Mmon,ANA_CHO_MM_19612018_1yr(:,1),userinput_ax4);

% Variability ranges
ANA_CHO_HQ_19211960_Hmonvar = [1:12,12:-1:1];
ANA_CHO_HQ_19211960_Hmonvar(2,:) = [prctile(userinput_frl*ANA_CHO_HQ_19211960_Hmon,5),prctile(userinput_frl*ANA_CHO_HQ_19211960_Hmon(:,end:-1:1),90)];
ANA_CHO_HQ_19612018_Hmonvar = [1:12,12:-1:1];
ANA_CHO_HQ_19612018_Hmonvar(2,:) = [prctile(userinput_frl*ANA_CHO_HQ_19612018_Hmon,5),prctile(userinput_frl*ANA_CHO_HQ_19612018_Hmon(:,end:-1:1),90)];
ANA_CHO_MM_19612018_Hmonvar = [1:12,12:-1:1];
ANA_CHO_MM_19612018_Hmonvar(2,:) = [prctile(userinput_frl*ANA_CHO_MM_19612018_Hmon,5),prctile(userinput_frl*ANA_CHO_MM_19612018_Hmon(:,end:-1:1),90)];

% Summary table of original, diverted, and modified flows
disp('Summary table of original, diverted, and modified flows')
% Matrix ANA_CHO_MAMVOL_mean:
% column 1:     mean annual river flow (hydrological year)
% column 2:     wet season river flow (hydrological year)
% column 3:     mean potential diversion (using userinput_fdv % diversion)
% column 4:     mean remaming flow (approximately column 1 - column 2)
% column 5:     mean modified flow (using userinput_frr recovery rate)
% row 1:        units in m3/s
% row 2:        units in million m3/yr
% row 3:        units in mm/yr

% million m3/day
aux_mean = nanmean(ANA_CHO_MAMFLOW_19612018_1day);
ANA_CHO_MAMVOL_mean(1) = nanmean(ANA_CHO_HQ_19612018_1day(:,1));
ANA_CHO_MAMVOL_mean(1,2:5) = aux_mean(1:4);
% million m3/yr
aux_sum = nansum(ANA_CHO_MAMFLOW_19612018_1day);
ANA_CHO_MAMVOL_mean(2,1) = nansum(ANA_CHO_HQ_19612018_1day(:,1))/n;
ANA_CHO_MAMVOL_mean(2,2:5) = aux_sum(1:4)/n;
% mm/yr
aux_factor = 1000/(ANA_CHO_AREA_station/1e6)/1e6/n;
ANA_CHO_MAMVOL_mean(3,1) = nansum(ANA_CHO_HQ_19612018_1day(:,2))/n;
ANA_CHO_MAMVOL_mean(3,2:4) = aux_sum(1:3)*aux_factor;
ANA_CHO_MAMVOL_mean(3,5) = aux_sum(4)/n;

% Clear auxiliar variables n
clear k n refmon_hy refmon_ws refyrs_hy aux_mean aux_sum aux_factor

disp('Done')
disp(' ')
%% 14. Modified flow drought analysis
disp('14. Performing drought analysis for grey and green+grey flows')

drought_aggr_Q = 1;
% aggreg = number 0 or 1 specyfing further method details as follows:
%           1: smooth series with 30-day moving average 
%           0: use original series (to evaluate Q only)

% Streamflow drought indices compared against 1920-1940 threshold
[ANA_CHO_HQ_DIndices_19211960_vs_19211960] = iMHEA_DroughtIndices(ANA_CHO_HQ_DDATE(and(ANA_CHO_HQ_DDATE>=userinput_ax3(1),ANA_CHO_HQ_DDATE<=userinput_ax3(2))),ANA_CHO_HQ_1day(and(ANA_CHO_HQ_DDATE>=userinput_ax3(1),ANA_CHO_HQ_DDATE<=userinput_ax3(2)),4),ANA_CHO_HQ_DThr_19211960,drought_aggr_Q,1);
[ANA_CHO_HQ_DIndices_19612018_vs_19211960] = iMHEA_DroughtIndices(ANA_CHO_HQ_DDATE(and(ANA_CHO_HQ_DDATE>=userinput_ax4(1),ANA_CHO_HQ_DDATE<=userinput_ax4(2))),ANA_CHO_HQ_1day(and(ANA_CHO_HQ_DDATE>=userinput_ax4(1),ANA_CHO_HQ_DDATE<=userinput_ax4(2)),4),ANA_CHO_HQ_DThr_19211960,drought_aggr_Q,1);
[ANA_CHO_MM_DIndices_19612018_vs_19211960] = iMHEA_DroughtIndices(ANA_CHO_HQ_19612018_DDATE(and(ANA_CHO_HQ_19612018_DDATE>=userinput_ax4(1),ANA_CHO_HQ_19612018_DDATE<=userinput_ax4(2))),ANA_CHO_MAMFLOW_19612018_1day(and(ANA_CHO_HQ_19612018_DDATE>=userinput_ax4(1),ANA_CHO_HQ_19612018_DDATE<=userinput_ax4(2)),5),ANA_CHO_HQ_DThr_19211960,drought_aggr_Q,1);

% Streamflow drought indices against their own thresholds
[ANA_CHO_HQ_DIndices_19612018,ANA_CHO_HQ_DThr_19612018] = iMHEA_Drought(ANA_CHO_HQ_DDATE(and(ANA_CHO_HQ_DDATE>=userinput_ax4(1),ANA_CHO_HQ_DDATE<=userinput_ax4(2))),ANA_CHO_HQ_1day(and(ANA_CHO_HQ_DDATE>=userinput_ax4(1),ANA_CHO_HQ_DDATE<=userinput_ax4(2)),4),drought_method,drought_aggr_Q,1);
[ANA_CHO_MM_DIndices_19211960,ANA_CHO_MM_DThr_19211960] = iMHEA_Drought(ANA_CHO_HQ_19612018_DDATE(and(ANA_CHO_HQ_19612018_DDATE>=userinput_ax4(1),ANA_CHO_HQ_19612018_DDATE<=userinput_ax4(2))),ANA_CHO_MAMFLOW_19612018_1day(and(ANA_CHO_HQ_19612018_DDATE>=userinput_ax4(1),ANA_CHO_HQ_19612018_DDATE<=userinput_ax4(2)),5),drought_method,drought_aggr_Q,1);

% Streamflow drought indices compared against 1960-2018 threshold
[ANA_CHO_MM_DIndices_19612018_vs_19612018] = iMHEA_DroughtIndices(ANA_CHO_HQ_19612018_DDATE(and(ANA_CHO_HQ_19612018_DDATE>=userinput_ax4(1),ANA_CHO_HQ_19612018_DDATE<=userinput_ax4(2))),ANA_CHO_MAMFLOW_19612018_1day(and(ANA_CHO_HQ_19612018_DDATE>=userinput_ax4(1),ANA_CHO_HQ_19612018_DDATE<=userinput_ax4(2)),5),ANA_CHO_HQ_DThr_19612018,drought_aggr_Q,1);

disp('Done')
disp(' ')
%% 15. Flow Duration Curve
disp('15. Extracting flow duration curves')

% Natural flow regime 1921 - 1960
disp('Calculating flow duration curves for 1921 - 1960')
% Define hydrological years seasons
refyrs_hy = year(userinput_ax3);
refmon_hy = month(userinput_ax3(1,:));
n = refyrs_hy(2)-refyrs_hy(1);
ANA_CHO_hy_19201960 = [datetime([(refyrs_hy(1):refyrs_hy(2)-1)',refmon_hy(1)*ones(n,1),ones(n,1)]),datetime([(refyrs_hy(1)+1:refyrs_hy(2))',refmon_hy(2)*ones(n,1),31*ones(n,1)])];
% Calculate annual flow duration curves (FDC)
[ANA_CHO_HQ_19211960_1day_FDC] = iMHEA_AnnualFDC(ANA_CHO_HQ_DDATE,ANA_CHO_HQ_1day(:,4),ANA_CHO_hy_19201960);

% Altered flow regime with grey and green infrastructure 1961 - present
disp('Calculating flow duration curves for 1961 - present')
% Define hydrological years seasons
refyrs_hy = year(userinput_ax4);
refmon_hy = month(userinput_ax4(1,:));
n = refyrs_hy(2)-refyrs_hy(1);
ANA_CHO_hy_19602018 = [datetime([(refyrs_hy(1):refyrs_hy(2)-1)',refmon_hy(1)*ones(n,1),ones(n,1)]),datetime([(refyrs_hy(1)+1:refyrs_hy(2))',refmon_hy(2)*ones(n,1),31*ones(n,1)])];
% Calculate annual flow duration curves (FDC)
[ANA_CHO_HQ_19612018_1day_FDC] = iMHEA_AnnualFDC(ANA_CHO_HQ_DDATE,ANA_CHO_HQ_1day(:,4),ANA_CHO_hy_19602018);
[ANA_CHO_MM_19612018_1day_FDC] = iMHEA_AnnualFDC(ANA_CHO_HQ_19612018_DDATE,ANA_CHO_MAMFLOW_19612018_1day(:,5),ANA_CHO_hy_19602018);

% Drought and demand thresholds
disp('Calculating flow duration curves for thresholds')
[~,~,~,~,ANA_CHO_HQ_DThr_19211960_FDC] = iMHEA_FDC(ANA_CHO_HQ_DThr_19211960);
[~,~,~,~,FT_RIMAC_HQ_dem_FDC] = iMHEA_FDC(FT_RIMAC_HQ_dem*1e9/(ANA_CHO_AREA_station*1e6)./eomday(2018,[9:12,1:8]));

% Variability ranges
ANA_CHO_HQ_19211960_1day_FDCvar = [ANA_CHO_HQ_19211960_1day_FDC(1:end,1)',ANA_CHO_HQ_19211960_1day_FDC(end:-1:1,1)'];
ANA_CHO_HQ_19211960_1day_FDCvar(2,1:100) = prctile(ANA_CHO_HQ_19211960_1day_FDC(:,2:end)',5);
ANA_CHO_HQ_19211960_1day_FDCvar(2,200:-1:101) = prctile(ANA_CHO_HQ_19211960_1day_FDC(:,2:end)',90);
ANA_CHO_HQ_19612018_1day_FDCvar = [ANA_CHO_HQ_19612018_1day_FDC(1:end,1)',ANA_CHO_HQ_19612018_1day_FDC(end:-1:1,1)'];
ANA_CHO_HQ_19612018_1day_FDCvar(2,1:100) = prctile(ANA_CHO_HQ_19612018_1day_FDC(:,2:end)',5);
ANA_CHO_HQ_19612018_1day_FDCvar(2,200:-1:101) = prctile(ANA_CHO_HQ_19612018_1day_FDC(:,2:end)',90);
ANA_CHO_MM_19612018_1day_FDCvar = [ANA_CHO_MM_19612018_1day_FDC(1:end,1)',ANA_CHO_MM_19612018_1day_FDC(end:-1:1,1)'];
ANA_CHO_MM_19612018_1day_FDCvar(2,1:100) = prctile(ANA_CHO_MM_19612018_1day_FDC(:,2:end)',5);
ANA_CHO_MM_19612018_1day_FDCvar(2,200:-1:101) = prctile(ANA_CHO_MM_19612018_1day_FDC(:,2:end)',90);

% Clear auxiliar variables
clear n refmon_hy refmon_hy refyrs_hy

disp('Done')
disp(' ')
%% 16. Figure: Regional effects
disp('16. Creating Figure of regional effects')

% Initialise variables
offset_loc = [-0.1 0.1 0.02];
lim2 = 350;

fig_4 = figure;
set(fig_4,'Renderer','painters')

% Supply - Demand plot for discharge in Lima
disp('Ploting monthly regional hydrograph')
hold on
% Variability ranges
patch(ANA_CHO_HQ_19211960_Hmonvar(1,:),ANA_CHO_HQ_19211960_Hmonvar(2,:),[0 0.45 0.74],'Edgecolor','none','FaceAlpha',.1)
patch(ANA_CHO_HQ_19612018_Hmonvar(1,:),ANA_CHO_HQ_19612018_Hmonvar(2,:),'r','Edgecolor','none','FaceAlpha',.1)
patch(ANA_CHO_MM_19612018_Hmonvar(1,:),ANA_CHO_MM_19612018_Hmonvar(2,:),[0 0.5 0],'Edgecolor','none','FaceAlpha',.1)
% Mean volumes
nat = plot(userinput_frl*mean(ANA_CHO_HQ_19211960_Hmon),'Linewidth',1,'color',[0 0.45 0.74]);
alt = plot(userinput_frl*mean(ANA_CHO_HQ_19612018_Hmon),'Linewidth',1,'color','r');
grn = plot(userinput_frl*mean(ANA_CHO_MM_19612018_Hmon),'Linewidth',1,'color',[0 0.5 0]);
% D = plot(userinput_frl*ANA_CHO_QDThrmon_19211960([refmon_hy:12,1:refmon_hy-1])*ANA_CHO_AREA/1e3,':r','Linewidth',1);
d = plot(FT_RIMAC_HQ_dem,'--k','Linewidth',1);% ,'color',[0.85 0.33 0.1]);
% Plot settings
set(gca,'Xlim',[0.5 12.5],...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*lim2,...
    'XTick',1:2:12,...
    'XTickLabel',MonthLabels(1:2:end),...
    'TickDir','out') % ,'FontSize',16)
box off
% Labels
% text(offset_loc(1)*2,1+offset_loc(2)/2,'c','Units','normalized','FontWeight','bold' ,'FontSize',16)
ylabel('million m^{3} month^{-1}')
% Legend
legend('AutoUpdate','off');
legend('location','NorthWest')
legend('boxoff')
legend([nat,alt,grn,d],'Natural flow regime in Rimac',...
    'Artificial hydrological regulation (grey infrastructure)',...
    'Potential effect of infiltration system (green+grey infrastructure)',...
    'Lima city''s water demand')
clear nat alt grn d D
% Plot seasons
plot([1 12],[lim2 lim2],':k','Linewidth',1)
plot([4 8],[lim2 lim2],'k','LineWidth',2)
text(6,lim2+10,'Wet season','HorizontalAlignment','center')% ,'FontSize',14)
text(10,lim2+10,'Dry season','HorizontalAlignment','center')% ,'FontSize',14)

% Export Figure
disp('Exporting the Figure to pdf')
% Set the paper size [width height]
set(fig_4,'PaperSize',[15 15]);
% set(fig_4,'PaperOrientation','landscape')
print(fig_4,'Fig4_export','-dpdf','-fillpage')

% Clear auxiliar variables
clear fig_4 offset_loc lim2 refmon_hy

disp('Done')
disp(' ')
%% 17. Sensitivity analysis for the recovery rate
disp('17. Sensitivity analysis for the recovery rate')

% Residence times at 5-min scale
l1dbf = iMHEA_HMT_Tracers_Int_1day(:,4);

% Sensitivity analysis of the recovery rate +/- 0.10 and optimistic 0.70
userinput_frs = [userinput_frr-0.1 userinput_frr userinput_frr+0.1 0.70];

% Number of days during the wet season
refyrs_hy = year(userinput_ax4);
n = refyrs_hy(2)-refyrs_hy(1);

% Initialise variables
ANA_CHO_SA_19612018_1day = nan(length(ANA_CHO_MAMFLOW_19612018_1day(:,3)),4);
ANA_CHO_SA_19612018_1yr = cell(1,4);
ANA_CHO_SA_19612018_Mmon = cell(1,4);
ANA_CHO_SA_19612018_Hmon = cell(1,4);
ANA_CHO_SA = nan(length(userinput_frs),6);
ANA_CHO_SA_19612018_1day_FDC = cell(1,4);

% Apply diverted flows, residence times, and recovery rate
for j = 1:length(userinput_frs)
    disp('Recovery rate used [fraction]:')
    disp(userinput_frs(j))
    [ANA_CHO_SA_19612018_1day(:,j)] = iMHEA_ModifiedFlows(ANA_CHO_MAMFLOW_19612018_1day(:,3),ANA_CHO_MAMFLOW_19612018_1day(:,2),l1dbf,userinput_frs(j));
    % [ANA_CHO_SA_19612018_1day(:,j)] = iMHEA_ModifiedFlowsBF(ANA_CHO_MAMFLOW_19612018_1day(:,3),ANA_CHO_MAMFLOW_19612018_1day(:,2),ANA_CHO_HQ_1day(:,5),l1dbf,userinput_frs(j));
    % Monthly modified hydrograph
    [~,ANA_CHO_SA_19612018_1yr{j},~,~,ANA_CHO_SA_19612018_Mmon{j}] = iMHEA_MonthlyRain(ANA_CHO_HQ_19612018_DDATE,ANA_CHO_SA_19612018_1day(:,j));
    % Reorganise in hydrological years
    ANA_CHO_SA_19612018_Hmon{j} = iMHEA_HydroYear(ANA_CHO_SA_19612018_Mmon{j},ANA_CHO_SA_19612018_1yr{j}(:,1),userinput_ax4);

    % Summary table of recovered flows and increased dry season flow percentage
    % Matrix ANA_CHO_SA:
    % column 1:     Recovery rate (-)
    % column 2:     Recovered volume (m3/yr)
    % column 3:     Long-term monthly average dry season increase (%)
    % column 4:     Long-term monthly max dry season increase (%)
    % column 5:     Maximum monthly mean dry season increase (%)
    % column 6:     Average monthly mean dry season increase (%)
    % column 7:     Minimum monthly mean dry season increase (%)
    
    ANA_CHO_SA(j,1) = userinput_frs(j);    
    auxvol = ANA_CHO_SA_19612018_1day(:,j) - ANA_CHO_MAMFLOW_19612018_1day(:,3);
    ANA_CHO_SA(j,2) = nansum(auxvol)/n;
    auxdiff = (ANA_CHO_SA_19612018_Hmon{j} - ANA_CHO_HQ_19612018_Hmon)./ANA_CHO_HQ_19612018_Hmon;
    ANA_CHO_SA(j,3) = mean(auxdiff(auxdiff>0))*100;
    ANA_CHO_SA(j,4) = max(auxdiff(:))*100;
    auxmonth = mean(auxdiff);
    ANA_CHO_SA(j,5) = min(auxmonth(auxmonth>0))*100;
    ANA_CHO_SA(j,6) = mean(auxmonth(auxmonth>0))*100;
    ANA_CHO_SA(j,7) = max(auxmonth)*100;
    % Calculate flow duration curves (FDC)
    ANA_CHO_SA_19612018_1day_FDC{j} = iMHEA_AnnualFDC(ANA_CHO_HQ_19612018_DDATE,ANA_CHO_SA_19612018_1day(:,j)*1e9/(ANA_CHO_AREA_station*1e6),ANA_CHO_hy_19602018);
end

% Clear auxiliar variables
clear j refyrs_hy n auxdiff auxmonth auxvol

disp('Done')
disp(' ')
%% 18. Figure: Extended regional effects
disp('18. Creating Figure of extended regional effects')

% Initialise variables
offset_loc = [-0.1 0.1 0.02];
lim1 = 100;
lim2 = 350;
% Extract data for the analysis periods
refmon_hy = month(userinput_ax4(1));

fig_S7 = figure;
set(fig_S7,'Renderer','painters')

% Monthly effects in Huamantanga
disp('Ploting monthly local hydrograph')
subplot(2,2,1)
hold on
% Variability ranges
patch(iMHEA_HMT_HQ_Hmonvar(1,:),iMHEA_HMT_HQ_Hmonvar(2,:),[0 0.45 0.74],'Edgecolor','none','FaceAlpha',.1)
patch(iMHEA_HMT_MM_Hmonvar(1,:),iMHEA_HMT_MM_Hmonvar(2,:),[0 0.5 0],'Edgecolor','none','FaceAlpha',.1)
% Mean volumes
b = plot(1:12,mean(iMHEA_HMT_HQ_Hmon),'Linewidth',1,'color',[0 0.45 0.74]);
c = plot(1:12,mean(iMHEA_HMT_MM_Hmon),'Linewidth',1,'color',[0 0.5 0]);
% Plot settings
set(gca,'Xlim',[0.5 12.5],...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*lim1,...
    'XTick',1:2:12,...
    'XTickLabel',MonthLabels(1:2:end),...
    'TickDir','out') % ,'FontSize',16)
% set(gca,'Xticklabel',[])
box off
% Labels
text(offset_loc(1)*2,1+offset_loc(2)/2,'a','Units','normalized','FontWeight','bold' ,'FontSize',16)
ylabel('mm month^{-1}') % ylabel('mm month^{-1}')
% Legend
legend('AutoUpdate','off');
legend('location','NorthWest')
legend('boxoff')
legend([b,c],'Original average discharge in Huamantanga',...
    'Potential effect of infiltration system (local effects)')
clear b c
% Plot seasons
plot([1 12],[lim1 lim1],':k','Linewidth',1)
plot([4 8],[lim1 lim1],'k','LineWidth',2)
text(6,lim1+5,'Wet season','HorizontalAlignment','center')% ,'FontSize',14)
text(10,lim1+5,'Dry season','HorizontalAlignment','center')% ,'FontSize',14)

% Supply - Demand plot for discharge in Lima + Grey infrastructure
disp('Ploting monthly regional hydrograph')
subplot(2,2,2)
hold on
% Variability ranges
patch(ANA_CHO_HQ_19211960_Hmonvar(1,:),ANA_CHO_HQ_19211960_Hmonvar(2,:),[0 0.45 0.74],'Edgecolor','none','FaceAlpha',.1)
patch(ANA_CHO_HQ_19612018_Hmonvar(1,:),ANA_CHO_HQ_19612018_Hmonvar(2,:),'r','Edgecolor','none','FaceAlpha',.1)
% patch(ANA_CHO_MM_19612018_Hmonvar(1,:),ANA_CHO_MM_19612018_Hmonvar(2,:),[0 0.5 0],'Edgecolor','none','FaceAlpha',.1)
% Mean volumes
nat = plot(userinput_frl*mean(ANA_CHO_HQ_19211960_Hmon),'Linewidth',1,'color',[0 0.45 0.74]);
alt = plot(userinput_frl*mean(ANA_CHO_HQ_19612018_Hmon),'Linewidth',1,'color','r');
% grn = plot(userinput_frl*mean(ANA_CHO_MM_19612018),'Linewidth',1,'color',[0 0.5 0]);
D = plot(userinput_frl*ANA_CHO_HQ_DThrmon_19211960([refmon_hy:12,1:refmon_hy-1])*ANA_CHO_AREA_station/1e3,':k','Linewidth',1);
d = plot(FT_RIMAC_HQ_dem,'--','Linewidth',2,'color',[0.85 0.33 0.1]);
% Plot settings
set(gca,'Xlim',[0.5 12.5],...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*lim2,...
    'XTick',1:2:12,...
    'XTickLabel',MonthLabels(1:2:end),...
    'TickDir','out') % ,'FontSize',16)
box off
% Labels
text(offset_loc(1)*2,1+offset_loc(2)/2,'b','Units','normalized','FontWeight','bold' ,'FontSize',16)
ylabel('million m^{3} month^{-1}')
% Legend
legend('AutoUpdate','off');
legend('location','NorthWest')
legend('boxoff')
legend([nat,alt,D,d],'Natural flow regime in Rimac',...
    'Artificial hydrological regulation (grey infrastructure)',...
    'Rimac river flow drought threshold',...
    'Lima city''s water demand')
clear nat alt grn d D
% Plot seasons
plot([1 12],[lim2 lim2],':k','Linewidth',1)
plot([4 8],[lim2 lim2],'k','LineWidth',2)

% Supply - Demand plot for discharge in Lima + Green infrastructure
disp('Ploting monthly regional hydrograph')
subplot(2,2,3)
hold on
% Variability ranges
patch(ANA_CHO_HQ_19211960_Hmonvar(1,:),ANA_CHO_HQ_19211960_Hmonvar(2,:),[0 0.45 0.74],'Edgecolor','none','FaceAlpha',.1)
% patch(ANA_CHO_HQ_19612018_Hmonvar(1,:),ANA_CHO_HQ_19612018_Hmonvar(2,:),'r','Edgecolor','none','FaceAlpha',.1)
patch(ANA_CHO_MM_19612018_Hmonvar(1,:),ANA_CHO_MM_19612018_Hmonvar(2,:),[0 0.5 0],'Edgecolor','none','FaceAlpha',.1)
% Mean volumes
nat = plot(userinput_frl*mean(ANA_CHO_HQ_19211960_Hmon),'Linewidth',1,'color',[0 0.45 0.74]);
% alt = plot(userinput_frl*mean(ANA_CHO_HV_19612018),'Linewidth',1,'color','r');
grn = plot(userinput_frl*mean(ANA_CHO_MM_19612018_Hmon),'Linewidth',1,'color',[0 0.5 0]);
D = plot(userinput_frl*ANA_CHO_HQ_DThrmon_19211960([refmon_hy:12,1:refmon_hy-1])*ANA_CHO_AREA_station/1e3,':k','Linewidth',1);
d = plot(FT_RIMAC_HQ_dem,'--','Linewidth',2,'color',[0.85 0.33 0.1]);
% Plot settings
set(gca,'Xlim',[0.5 12.5],...
    'Ylim',[-offset_loc(3) 1+offset_loc(3)]*lim2,...
    'XTick',1:2:12,...
    'XTickLabel',MonthLabels(1:2:end),...
    'TickDir','out') % ,'FontSize',16)
box off
% Labels
text(offset_loc(1)*2,1+offset_loc(2)/2,'c','Units','normalized','FontWeight','bold' ,'FontSize',16)
ylabel('million m^{3} month^{-1}')
% Legend
legend('AutoUpdate','off');
legend('location','NorthWest')
legend('boxoff')
legend([nat,grn,D,d],'Natural flow regime in Rimac',...
    'Potential effect of infiltration system (green+grey infrastructure)',...
    'Rimac river flow drought threshold',...
    'Lima city''s water demand')
clear nat alt grn d D
% Plot seasons
plot([1 12],[lim2 lim2],':k','Linewidth',1)
plot([4 8],[lim2 lim2],'k','LineWidth',2)
% Add sensitivity analyis
for j = 1:length(userinput_frs)
    plot(userinput_frl*mean(ANA_CHO_SA_19612018_Hmon{j}),'Linewidth',1,'color',userinput_frs(j)*[.75 1 .75])
    plot([10 11],(150+150*userinput_frs(j))*[1 1],'color',userinput_frs(j)*[.75 1 .75])
    text(11,150+150*userinput_frs(j),[' ',num2str(userinput_frs(j))],'FontSize',8)
end
    text(10,150+150*(userinput_frs(j)+0.1),'Recovery rate','FontSize',8)
    
% Flow Duration Curves
disp('Ploting flow duration curves')
subplot(2,2,4)
hold on
% Variability ranges
patch(ANA_CHO_HQ_19211960_1day_FDCvar(1,:),ANA_CHO_HQ_19211960_1day_FDCvar(2,:),[0 0.45 0.74],'Edgecolor','none','FaceAlpha',.1)
patch(ANA_CHO_HQ_19612018_1day_FDCvar(1,:),ANA_CHO_HQ_19612018_1day_FDCvar(2,:),'r','Edgecolor','none','FaceAlpha',.1)
patch(ANA_CHO_MM_19612018_1day_FDCvar(1,:),ANA_CHO_MM_19612018_1day_FDCvar(2,:),[0 0.5 0],'Edgecolor','none','FaceAlpha',.1)
% Mean FDCs
nat = semilogy(ANA_CHO_HQ_19211960_1day_FDC(:,1),mean(userinput_frl*ANA_CHO_HQ_19211960_1day_FDC(:,2:end),2),'Linewidth',1,'color',[0 0.45 0.74]);
alt = semilogy(ANA_CHO_HQ_19612018_1day_FDC(:,1),mean(userinput_frl*ANA_CHO_HQ_19612018_1day_FDC(:,2:end),2),'Linewidth',1,'color','r');
grn = semilogy(ANA_CHO_MM_19612018_1day_FDC(:,1),mean(userinput_frl*ANA_CHO_MM_19612018_1day_FDC(:,2:end),2),'Linewidth',1,'color',[0 0.5 0]);
D = semilogy(ANA_CHO_HQ_DThr_19211960_FDC(:,1),userinput_frl*ANA_CHO_HQ_DThr_19211960_FDC(:,2),':k','Linewidth',1);
d = semilogy(FT_RIMAC_HQ_dem_FDC(:,1),FT_RIMAC_HQ_dem_FDC(:,2),'--','Linewidth',2,'color',[0.85 0.33 0.1]);
% Plot settings
set(gca,'Xlim',[-1 101],'Yscale','log',...
    'Ylim',[8e-2 1.2e1],...
    'TickDir','out') %'Ylim',[-offset_loc(3) 1+offset_loc(3)]*lim2,'FontSize',16)
box off
% Labels
text(offset_loc(1)*2,1+offset_loc(2)/2,'d','Units','normalized','FontWeight','bold' ,'FontSize',16)
ylabel('mm day^{-1}')
xlabel('% exceedance')
% Legend
legend('AutoUpdate','off');
legend('location','SouthWest')
legend('boxoff')
legend([nat,alt,grn,D,d],'Natural flow regime in Rimac',...
    'Artificial hydrological regulation (grey infrastructure)',...
    'Potential effect of infiltration system (green+grey infrastructure)',...
    'Rimac river flow drought threshold',...
    'Lima city''s water demand')
clear nat alt grn D d
for j = 1:length(userinput_frs)
    semilogy(ANA_CHO_SA_19612018_1day_FDC{j}(:,1),mean(userinput_frl*ANA_CHO_SA_19612018_1day_FDC{j}(:,2:end),2),'Linewidth',1,'color',userinput_frs(j)*[.75 1 .75])
%     plot([80 90],(1+j*0.5)*[1 1],'color',userinput_frs(j)*[.75 1 .75])
%     text(90,1+j*0.5,[' ',num2str(userinput_frs(j))],'FontSize',8)
    plot([80 90],(exp((userinput_frs(j)-0.15)*2))*[1 1],'color',userinput_frs(j)*[.75 1 .75])
    text(90,exp((userinput_frs(j)-0.15)*2),[' ',num2str(userinput_frs(j))],'FontSize',8)
end
    text(80,exp((userinput_frs(j)-0.05)*2),'Recovery rate','FontSize',8)

% Export Figure
disp('Exporting the Figure to pdf')
% Set the paper size [width height]
set(fig_S7,'PaperSize',[30 30]);
% set(fig_S7,'PaperOrientation','landscape')
print(fig_S7,'FigS7_export','-dpdf','-fillpage')

% Clear auxiliar variables
clear fig_S7 offset_loc refhyv refhydt refmon_hy refyrs lim1 lim2 j

disp('Done')
disp(' ')
%% 19. Dynamic plot (animation) at daily scale
disp('19. Generating dynamic plot (animation) at daily scale')

% Plot Rimac daily time series
h = figure;
filename = 'testAnimated_1day.gif';

hold on
b = plot(ANA_CHO_HQ_19612018_DDATE,ANA_CHO_HQ_19612018_1day(:,1),'Linewidth',1,'color',[0 0.45 0.74]);
c = plot(ANA_CHO_HQ_19612018_DDATE(1),ANA_CHO_MAMFLOW_19612018_1day(1,4),'Linewidth',1,'color',[0 0.5 0]);
d = plot([ANA_CHO_HQ_19612018_DDATE(1) ANA_CHO_HQ_19612018_DDATE(1)],[0 20],':k','Linewidth',1);
% Plot settings
set(gca,'Xlim',[ANA_CHO_HQ_19612018_DDATE(1) ANA_CHO_HQ_19612018_DDATE(366)],...
    'Ylim',[0 20],...
    'TickDir','out')
box off
% Labels
ylabel('m^{3} s^{-1}')
xlabel('Date')
legend('AutoUpdate','off');
legend('location','NorthWest')
legend('boxoff')
legend([b,c],'Discharge without diversion','Discharge with diversion')
drawnow

% Capture the plot as an image 
frame = getframe(h); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 

% Starting animation
for i = 2:15:183
    set(b,'Visible','off')
    b = plot(ANA_CHO_HQ_19612018_DDATE(i-1:end),ANA_CHO_HQ_19612018_1day(i-1:end,1),'Linewidth',1,'color',[0 0.45 0.74]);
    set(b,'Visible','on')
    set(c,'Visible','off')
    c = plot(ANA_CHO_HQ_19612018_DDATE(1:i),ANA_CHO_MAMFLOW_19612018_1day(1:i,4),'Linewidth',1,'color',[0 0.5 0]);
    set(c,'Visible','on')
    set(d,'Visible','off')
    d = plot([ANA_CHO_HQ_19612018_DDATE(i) ANA_CHO_HQ_19612018_DDATE(i)],[0 20],':k');
    set(d,'Visible','on')
    drawnow
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end
% Dynamic animation
for i = 183:15:length(ANA_CHO_HQ_19612018_DDATE)-183
    set(b,'Visible','off')
    b = plot(ANA_CHO_HQ_19612018_DDATE(i-1:end),ANA_CHO_HQ_19612018_1day(i-1:end,1),'Linewidth',1,'color',[0 0.45 0.74]);
    set(b,'Visible','on')
    set(c,'Visible','off')
    c = plot(ANA_CHO_HQ_19612018_DDATE(1:i),ANA_CHO_MAMFLOW_19612018_1day(1:i,4),'Linewidth',1,'color',[0 0.5 0]);
    set(c,'Visible','on')
    set(d,'Visible','off')
    d = plot([ANA_CHO_HQ_19612018_DDATE(i) ANA_CHO_HQ_19612018_DDATE(i)],[0 20],':k');
    set(d,'Visible','on')
    set(gca,'Xlim',[ANA_CHO_HQ_19612018_DDATE(1)+i-183 ANA_CHO_HQ_19612018_DDATE(366)+i-183])
    drawnow
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end
% Finishing animation
for i = length(ANA_CHO_HQ_19612018_DDATE)-183:15:length(ANA_CHO_HQ_19612018_DDATE)
    set(b,'Visible','off')
    b = plot(ANA_CHO_HQ_19612018_DDATE(i-1:end),ANA_CHO_HQ_19612018_1day(i-1:end,1),'Linewidth',1,'color',[0 0.45 0.74]);
    set(b,'Visible','on')
    set(c,'Visible','off')
    c = plot(ANA_CHO_HQ_19612018_DDATE(1:i),ANA_CHO_MAMFLOW_19612018_1day(1:i,4),'Linewidth',1,'color',[0 0.5 0]);
    set(c,'Visible','on')
    set(d,'Visible','off')
    d = plot([ANA_CHO_HQ_19612018_DDATE(i) ANA_CHO_HQ_19612018_DDATE(i)],[0 20],':k');
    set(d,'Visible','on')
    drawnow
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end

% Clear auxiliar variables
clear b c d  i h filename frame im imind cm auxyr auxmo

disp('Done')
disp(' ')
%% 20. Dynamic plot (animation) at monthly scale
disp('20. Generating dynamic plot (animation) at monthly scale')

% Initialise variables
n = length(ANA_CHO_HQ_19612018_Hyr(:,1));
auxyr = repmat(ANA_CHO_HQ_19612018_Hyr(:,1)',12,1);
auxmo = repmat((1:12)',n,1);
ANA_CHO_HQ_19612018_MDATE = datetime([auxyr(:),auxmo-4,15*ones(n*12,1)]);
ANA_CHO_HQ_19612018_1mo = ANA_CHO_HQ_19612018_Hmon';
ANA_CHO_HQ_19612018_1mo = ANA_CHO_HQ_19612018_1mo(:);
ANA_CHO_MM_19612018_1mo = ANA_CHO_MM_19612018_Hmon';
ANA_CHO_MM_19612018_1mo = ANA_CHO_MM_19612018_1mo(:);

% Plot Rimac daily time series
h = figure;
filename = 'testAnimated_1mo.gif';

hold on
b = plot(ANA_CHO_HQ_19612018_MDATE,ANA_CHO_HQ_19612018_1mo(:,1),'Linewidth',1,'color',[0 0.45 0.74]);
c = plot(ANA_CHO_HQ_19612018_MDATE(1),ANA_CHO_MM_19612018_1mo(1,1),'Linewidth',1,'color',[0 0.5 0]);
d = plot([ANA_CHO_HQ_19612018_MDATE(1) ANA_CHO_HQ_19612018_MDATE(1)],[0 400],':k','Linewidth',1);
% Plot settings
set(gca,'Xlim',[ANA_CHO_HQ_19612018_MDATE(1) ANA_CHO_HQ_19612018_MDATE(36)],...
    'Ylim',[0 400],...
    'TickDir','out')
box off
% Labels
ylabel('mm month^{-1}')
xlabel('Date')
title('Potential effect of infiltration systems')
%Legend
legend('AutoUpdate','off');
legend('location','NorthWest')
legend('boxoff')
legend([b,c],'Discharge without diversion','Discharge with diversion')
drawnow

% Capture the plot as an image 
frame = getframe(h); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 

% Starting animation
for i = 2:18
    set(b,'Visible','off')
    b = plot(ANA_CHO_HQ_19612018_MDATE(i-1:end),ANA_CHO_HQ_19612018_1mo(i-1:end,1),'Linewidth',1,'color',[0 0.45 0.74]);
    set(b,'Visible','on')
    set(c,'Visible','off')
    c = plot(ANA_CHO_HQ_19612018_MDATE(1:i),ANA_CHO_MM_19612018_1mo(1:i,1),'Linewidth',1,'color',[0 0.5 0]);
    set(c,'Visible','on')
    set(d,'Visible','off')
    d = plot([ANA_CHO_HQ_19612018_MDATE(i) ANA_CHO_HQ_19612018_MDATE(i)],[0 400],':k');
    set(d,'Visible','on')
    drawnow
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end
% Dynamic animation
for i = 18:length(ANA_CHO_HQ_19612018_MDATE)-18
    set(b,'Visible','off')
    b = plot(ANA_CHO_HQ_19612018_MDATE(i-1:end),ANA_CHO_HQ_19612018_1mo(i-1:end,1),'Linewidth',1,'color',[0 0.45 0.74]);
    set(b,'Visible','on')
    set(c,'Visible','off')
    c = plot(ANA_CHO_HQ_19612018_MDATE(1:i),ANA_CHO_MM_19612018_1mo(1:i,1),'Linewidth',1,'color',[0 0.5 0]);
    set(c,'Visible','on')
    set(d,'Visible','off')
    d = plot([ANA_CHO_HQ_19612018_MDATE(i) ANA_CHO_HQ_19612018_MDATE(i)],[0 400],':k');
    set(d,'Visible','on')
    set(gca,'Xlim',[ANA_CHO_HQ_19612018_MDATE(1)+(365/12)*(i-18) ANA_CHO_HQ_19612018_MDATE(36)+(365/12)*(i-18)])
    drawnow
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end
% Finishing animation
for i = length(ANA_CHO_HQ_19612018_MDATE)-18:length(ANA_CHO_HQ_19612018_MDATE)
    set(b,'Visible','off')
    b = plot(ANA_CHO_HQ_19612018_MDATE(i-1:end),ANA_CHO_HQ_19612018_1mo(i-1:end,1),'Linewidth',1,'color',[0 0.45 0.74]);
    set(b,'Visible','on')
    set(c,'Visible','off')
    c = plot(ANA_CHO_HQ_19612018_MDATE(1:i),ANA_CHO_MM_19612018_1mo(1:i,1),'Linewidth',1,'color',[0 0.5 0]);
    set(c,'Visible','on')
    set(d,'Visible','off')
    d = plot([ANA_CHO_HQ_19612018_MDATE(i) ANA_CHO_HQ_19612018_MDATE(i)],[0 400],':k');
    set(d,'Visible','on')
    drawnow
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end

% Clear auxiliar variables
clear b c d  i h filename frame im imind cm auxyr auxmo

disp('Done')
disp(' ')
%% 21. Save results
disp('21. Saving results from the analysis')

save iMHEA_MAMANTEO_results

disp('Done')
disp(' ')
