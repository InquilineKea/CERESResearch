
load LatWeights.mat
%zonal weights from http://ceres.larc.nasa.gov/data/zone_weights_lou.txt
ncdisp('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc');
rlutcs = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','rlutcs');
rsutcs = ncread('rsutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc', 'rsutcs');
rsdt = ncread('rsdt_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','rsdt');
rsut = ncread('rsut_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','rsut');
rlut = ncread('rlut_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','rlut');
SW = rsdt-rsut;
LW = -rlut;
SWCF = rsutcs-rsut; %what about shortwave down forcing? since this is what clouds should reflect...
%rsutcs < rsut. less shortwave up in clearsky? because less shortwave is
%reflected in clearsky..
LWCF = rlutcs-rlut; 
%rlutcs > rlut. more longwave up in clearsky
netcs = rsdt - rsutcs - rlutcs; %isnt this just clear sky??
net = rsdt - rsut - rlut;
lat = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');
lon = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lon');
time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));
Net = permute(net,[2 1 3]);
SWCF = permute(SWCF,[2 1 3]);
LWCF = permute(LWCF,[2 1 3]);
SW = permute(SW,[2 1 3]);
LW = permute(LW,[2 1 3]);
SWclear = permute(rsdt-rsutcs,[2 1 3]);
LWclear = permute(-rlutcs,[2 1 3]);
netclear = permute(netcs,[2 1 3]);
TotalCloudForcing = SWCF + LWCF;
clearvars net rsut rlut rsutcs rlutcs

NetClimatologySubtracted = GlobalValuesMinusClimatology(Net,'Net');
SWClimatologySubtracted = GlobalValuesMinusClimatology(SW,'SW');
LWClimatologySubtracted = GlobalValuesMinusClimatology(LW,'LW');
SWCFClimatologySubtracted = GlobalValuesMinusClimatology(SWCF,'SWCF');
LWCFClimatologySubtracted = GlobalValuesMinusClimatology(LWCF,'LWCF');
netclearClimatologySubtracted = GlobalValuesMinusClimatology(netclear,'netclear');
LWclearClimatologySubtracted = GlobalValuesMinusClimatology(LWclear,'LWclear');
SWclearClimatologySubtracted = GlobalValuesMinusClimatology(SWclear,'SWclear');
TotalCloudForcingClimatologySubtracted = GlobalValuesMinusClimatology(TotalCloudForcing,'TotalCloudForcing');

%%%%%%%%%%%%%%%%%%%%%

load ClimatologySubtractedVars.mat

%%%%%%%%%%%%%%%%%%%%%%%what about specific months only?
%then indices should be changed...

%SPRING
SpringIndices = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(1:3,1,13);
%WINTER
WinterIndices = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(10:12,1,13);
%AUTUMN
AutumnIndices = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(7:9,1,13);
%SUMMER
SummerIndices = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(4:6,1,13);

ensoTime = ncread('enso.cdf','T');
NINO12 = ncread('enso.cdf','NINO12');
NINO12 = NINO12(find(ensoTime > 482 & ensoTime < 638));
NINO3 = ncread('enso.cdf','NINO3');
NINO3 = NINO3(find(ensoTime > 482 & ensoTime < 638));
NINO34 = ncread('enso.cdf','NINO34');
NINO34 = NINO34(find(ensoTime > 482 & ensoTime < 638));
NINO4 = ncread('enso.cdf','NINO4');
NINO4 = NINO4(find(ensoTime > 482 & ensoTime < 638));
SAOTime = ncread('sao.cdf','T');
SAO = ncread('sao.cdf','anomaly');
SAO = SAO(find(SAOTime > 482 & SAOTime < 638));

NAM = dlmread('NAM.ascii');
NAO = dlmread('NAO.ascii');
SAM = dlmread('SAM.ascii');
PNA = dlmread('PNA.ascii');
T = dlmread('monthly.land_ocean.90S.90N.df_1901-2000mean.dat');
%NAM(find(NAM(:,1) >= 2000),3)
%find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2)
NAM = NAM(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2),3);
NAO = NAO(find(NAO(:,1) == 2000 & NAO(:,2)==3):find(NAO(:,1) == 2013 & NAO(:,2)==2),3);
SAM = SAM(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2),3);
PNA = PNA(find(PNA(:,1) == 2000 & PNA(:,2)==3):find(PNA(:,1) == 2013 & PNA(:,2)==2),3);
T = T(find(T(:,1) == 2000 & T(:,2)==3):find(T(:,1) == 2013 & T(:,2)==2),3);

MonthlyFluxDepartures = struct('net',NetClimatologySubtracted,'SW',SWClimatologySubtracted,'LW',LWClimatologySubtracted,...
    'SWCF',SWCFClimatologySubtracted,'LWCF',LWCFClimatologySubtracted,'netclear',netclearClimatologySubtracted,'LWclear',...
    LWclearClimatologySubtracted,'SWclear',SWclearClimatologySubtracted,'TotalCloudForcing',TotalCloudForcingClimatologySubtracted,...
    'Precip',PrecipClimatologySubtracted);
fields = fieldnames(MonthlyFluxDepartures);
GlobalCorrMap(MonthlyFluxDepartures.(fields{1}),MonthlyFluxDepartures.(fields{3}))

for i = 1:length(fields)
    for j = i+1:length(fields)
        GlobalCorrMapVarNameInput(MonthlyFluxDepartures.(fields{i}),MonthlyFluxDepartures.(fields{j}),fields{i},fields{j});
    end
end

MonthlyIndicesTimeSeries = struct('NAM',NAM,'NAO',NAO,'NINO34',NINO34,'SAM',SAM,'PNA',PNA,'Temp',T);

HemisphericFluxesClimSubtracted(NetClimatologySubtracted,0,90,'Net') %do I really expect monthly-subtracted to be different...?
HemisphericFluxesClimSubtracted(NetClimatologySubtracted,0,30,'Net') %do I really expect monthly-subtracted to be different...?
HemisphericFluxesClimSubtracted(NetClimatologySubtracted,30,90,'Net') %do I really expect monthly-subtracted to be different...?

[NetNHMonthAllLats, NetSHMonthAllLats, NetNHMinusSHMonthAllLats] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,0,90,'Net');
[Net0_15N_Month, Net0_15S_Month, Net0_15Dif_Month] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,0,14.5,'Net');
[Net15_30N_Month, Net15_30S_Month, Net15_30Dif_Month] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,14.5,30,'Net');
[Net30_49N_Month, Net30_49S_Month, Net30_49Dif_Month] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,30,48.5,'Net');
[Net49_90N_Month, Net49_90S_Month, Net49_90Dif_Month] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,48.5,90,'Net');

plot([corr(NetNHMinusSHMonthAllLats',Net49_90S_Month') corr(NetNHMinusSHMonthAllLats',Net30_49S_Month') ...
    corr(NetNHMinusSHMonthAllLats',Net15_30S_Month') corr(NetNHMinusSHMonthAllLats',Net0_15S_Month') ...
    corr(NetNHMinusSHMonthAllLats',Net0_15N_Month') corr(NetNHMinusSHMonthAllLats',Net15_30N_Month')...
    corr(NetNHMinusSHMonthAllLats',Net30_49N_Month') corr(NetNHMinusSHMonthAllLats',Net49_90N_Month')],'rs','MarkerFaceColor','g','markersize', 20)
set(gca,'FontSize',20)
title('Correlation Coefficients between each Latitudinal Band and NH-SH net flux')
xlabel('Latitudinal Band')
ylabel('Correlation Coefficient')
grid on;
set(gca,'XTickLabel',{'90S-49S','49S-30S','30S-15S','15S-0S','0N-15N','15N-30N','30N-49N','49N-90N'})
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['NetLatitudinalNHMinusSHCorr','.png']);
hold off;

Corr1DTimeSeriesMap(NetClimatologySubtracted,Net0_15N_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net0_15S_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net0_15Dif_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net15_30N_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net15_30S_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net15_30Dif_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net30_49N_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net30_49S_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net30_49Dif_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net49_90N_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net49_90S_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,Net49_90Dif_Month);
Corr1DTimeSeriesMap(NetClimatologySubtracted,NetNHMinusSHMonthAllLats);
Corr1DTimeSeriesMap(NetClimatologySubtracted,NetNHMonthAllLats);
Corr1DTimeSeriesMap(NetClimatologySubtracted,NetSHMonthAllLats);

[SWNHMonthAllLats, SWSHMonthAllLats, SWNHMinusSHMonthAllLats] = HemisphericFluxesClimSubtracted(SWClimatologySubtracted,0,90,'SW');
[SW0_15N_Month, SW0_15S_Month, SW0_15Dif_Month] =HemisphericFluxesClimSubtracted(SWClimatologySubtracted,0,14.5,'SW');
[SW15_30N_Month, SW15_30S_Month, SW15_30Dif_Month] = HemisphericFluxesClimSubtracted(SWClimatologySubtracted,14.5,30,'SW');
[SW30_49N_Month, SW30_49S_Month, SW30_49Dif_Month] = HemisphericFluxesClimSubtracted(SWClimatologySubtracted,30,48.5,'SW');
[SW49_90N_Month, SW49_90S_Month, SW49_90Dif_Month] = HemisphericFluxesClimSubtracted(SWClimatologySubtracted,48.5,90,'SW');

plot([corr(SWNHMinusSHMonthAllLats',SW49_90S_Month') corr(SWNHMinusSHMonthAllLats',SW30_49S_Month') ...
    corr(SWNHMinusSHMonthAllLats',SW15_30S_Month') corr(SWNHMinusSHMonthAllLats',SW0_15S_Month') ...
    corr(SWNHMinusSHMonthAllLats',SW0_15N_Month') corr(SWNHMinusSHMonthAllLats',SW15_30N_Month')...
    corr(SWNHMinusSHMonthAllLats',SW30_49N_Month') corr(SWNHMinusSHMonthAllLats',SW49_90N_Month')],'rs','MarkerFaceColor','g','markersize', 20)
set(gca,'FontSize',20)
title('Correlation Coefficients between each Latitudinal Band and NH-SH SW flux')
xlabel('Latitudinal Band')
ylabel('Correlation Coefficient')
grid on;
set(gca,'XTickLabel',{'90S-49S','49S-30S','30S-15S','15S-0S','0N-15N','15N-30N','30N-49N','49N-90N'})
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['SWLatitudinalNHMinusSHCorr','.png']);
hold off;

Corr1DTimeSeriesMap(SWClimatologySubtracted,SW0_15N_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW0_15S_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW0_15Dif_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW15_30N_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW15_30S_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW15_30Dif_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW30_49N_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW30_49S_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW30_49Dif_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW49_90N_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW49_90S_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SW49_90Dif_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SWNHMinusSHMonthAllLats);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SWNHMonthAllLats);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SWSHMonthAllLats);

[LWNHMonthAllLats, LWSHMonthAllLats, LWNHMinusSHMonthAllLats] = HemisphericFluxesClimSubtracted(LWClimatologySubtracted,0,90,'LW');
[LW0_15N_Month, LW0_15S_Month, LW0_15Dif_Month] =HemisphericFluxesClimSubtracted(LWClimatologySubtracted,0,14.5,'LW');
[LW15_30N_Month, LW15_30S_Month, LW15_30Dif_Month] = HemisphericFluxesClimSubtracted(LWClimatologySubtracted,14.5,30,'LW');
[LW30_49N_Month, LW30_49S_Month, LW30_49Dif_Month] = HemisphericFluxesClimSubtracted(LWClimatologySubtracted,30,48.5,'LW');
[LW49_90N_Month, LW49_90S_Month, LW49_90Dif_Month] = HemisphericFluxesClimSubtracted(LWClimatologySubtracted,48.5,90,'LW');

plot([corr(LWNHMinusSHMonthAllLats',LW49_90S_Month') corr(LWNHMinusSHMonthAllLats',LW30_49S_Month') ...
    corr(LWNHMinusSHMonthAllLats',LW15_30S_Month') corr(LWNHMinusSHMonthAllLats',LW0_15S_Month') ...
    corr(LWNHMinusSHMonthAllLats',LW0_15N_Month') corr(LWNHMinusSHMonthAllLats',LW15_30N_Month')...
    corr(LWNHMinusSHMonthAllLats',LW30_49N_Month') corr(LWNHMinusSHMonthAllLats',LW49_90N_Month')],'rs','MarkerFaceColor','g','markersize', 20)
set(gca,'FontSize',20)
title('Correlation Coefficients between each Latitudinal Band and NH-SH LW flux')
xlabel('Latitudinal Band')
ylabel('Correlation Coefficient')
grid on;
set(gca,'XTickLabel',{'90S-49S','49S-30S','30S-15S','15S-0S','0N-15N','15N-30N','30N-49N','49N-90N'})
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['LWLatitudinalNHMinusSHCorr','.png']);
hold off;

plot([corr(SWNHMinusSHMonthAllLats',LW49_90S_Month') corr(SWNHMinusSHMonthAllLats',LW30_49S_Month') ...
    corr(SWNHMinusSHMonthAllLats',LW15_30S_Month') corr(SWNHMinusSHMonthAllLats',LW0_15S_Month') ...
    corr(SWNHMinusSHMonthAllLats',LW0_15N_Month') corr(SWNHMinusSHMonthAllLats',LW15_30N_Month')...
    corr(SWNHMinusSHMonthAllLats',LW30_49N_Month') corr(SWNHMinusSHMonthAllLats',LW49_90N_Month')],'rs','MarkerFaceColor','g','markersize', 20)
set(gca,'FontSize',20)
title('Correlation Coefficients between each Latitudinal Band in SW and NH-SH LW flux')
xlabel('Latitudinal Band')
ylabel('Correlation Coefficient')
grid on;
set(gca,'XTickLabel',{'90S-49S','49S-30S','30S-15S','15S-0S','0N-15N','15N-30N','30N-49N','49N-90N'})
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['SWBandLWGlobalLatitudinalNHMinusSHCorr','.png']);
hold off;

plot([corr(NetNHMinusSHMonthAllLats',LW49_90S_Month') corr(NetNHMinusSHMonthAllLats',LW30_49S_Month') ...
    corr(NetNHMinusSHMonthAllLats',LW15_30S_Month') corr(NetNHMinusSHMonthAllLats',LW0_15S_Month') ...
    corr(NetNHMinusSHMonthAllLats',LW0_15N_Month') corr(NetNHMinusSHMonthAllLats',LW15_30N_Month')...
    corr(NetNHMinusSHMonthAllLats',LW30_49N_Month') corr(NetNHMinusSHMonthAllLats',LW49_90N_Month')],'rs','MarkerFaceColor','g','markersize', 20)
set(gca,'FontSize',20)
title('Correlation Coefficients between each Latitudinal Band in LW and NH-SH Net flux')
xlabel('Latitudinal Band')
ylabel('Correlation Coefficient')
grid on;
set(gca,'XTickLabel',{'90S-49S','49S-30S','30S-15S','15S-0S','0N-15N','15N-30N','30N-49N','49N-90N'})
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['LWBandNetGlobalLatitudinalNHMinusSHCorr','.png']);
hold off;

plot([corr(NetNHMinusSHMonthAllLats',SW49_90S_Month') corr(NetNHMinusSHMonthAllLats',SW30_49S_Month') ...
    corr(NetNHMinusSHMonthAllLats',SW15_30S_Month') corr(NetNHMinusSHMonthAllLats',SW0_15S_Month') ...
    corr(NetNHMinusSHMonthAllLats',SW0_15N_Month') corr(NetNHMinusSHMonthAllLats',SW15_30N_Month')...
    corr(NetNHMinusSHMonthAllLats',SW30_49N_Month') corr(NetNHMinusSHMonthAllLats',SW49_90N_Month')],'rs','MarkerFaceColor','g','markersize', 20)
set(gca,'FontSize',20)
title('Correlation Coefficients between each Latitudinal Band in SW and NH-SH Net flux')
xlabel('Latitudinal Band')
ylabel('Correlation Coefficient')
grid on;
set(gca,'XTickLabel',{'90S-49S','49S-30S','30S-15S','15S-0S','0N-15N','15N-30N','30N-49N','49N-90N'})
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['SWBandNetGlobalLatitudinalNHMinusSHCorr','.png']);
hold off;



Corr1DTimeSeriesMap(SWClimatologySubtracted,LW0_15N_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW0_15S_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW0_15Dif_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW15_30N_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW15_30S_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW15_30Dif_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW30_49N_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW30_49S_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW30_49Dif_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW49_90N_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW49_90S_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LW49_90Dif_Month);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LWNHMinusSHMonthAllLats);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LWNHMonthAllLats);
Corr1DTimeSeriesMap(SWClimatologySubtracted,LWSHMonthAllLats);

Corr1DTimeSeriesMap(LWClimatologySubtracted,LW0_15N_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW0_15S_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW0_15Dif_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW15_30N_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW15_30S_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW15_30Dif_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW30_49N_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW30_49S_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW30_49Dif_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW49_90N_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW49_90S_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LW49_90Dif_Month);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LWNHMinusSHMonthAllLats);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LWNHMonthAllLats);
Corr1DTimeSeriesMap(LWClimatologySubtracted,LWSHMonthAllLats);

[netclearNHMonthAllLats, netclearSHMonthAllLats, netclearNHMinusSHMonthAllLats] = HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,0,90,'netclear')
[netclear0_15N_Month, netclear0_15S_Month, netclear0_15Dif_Month] =HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,0,14.5,'netclear');
[netclear15_30N_Month, netclear15_30S_Month, netclear15_30Dif_Month] = HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,14.5,30,'netclear');
[netclear30_49N_Month, netclear30_49S_Month, netclear30_49Dif_Month] = HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,30,48.5,'netclear');
[netclear49_90N_Month, netclear49_90S_Month, netclear49_90Dif_Month] = HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,48.5,90,'netclear');

[SWCFNHMonthAllLats, SWCFSHMonthAllLats, SWCFNHMinusSHMonthAllLats] = HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,0,90,'SWCF')
[SWCF0_15N_Month, SWCF0_15S_Month, SWCF0_15Dif_Month] =HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,0,14.5,'SWCF');
[SWCF15_30N_Month, SWCF15_30S_Month, SWCF15_30Dif_Month] = HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,14.5,30,'SWCF');
[SWCF30_49N_Month, SWCF30_49S_Month, SWCF30_49Dif_Month] = HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,30,48.5,'SWCF');
[SWCF49_90N_Month, SWCF49_90S_Month, SWCF49_90Dif_Month] = HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,48.5,90,'SWCF');

[LWCFNHMonthAllLats, LWCFSHMonthAllLats, LWCFNHMinusSHMonthAllLats] = HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,0,90,'LWCF')
[LWCF0_15N_Month, LWCF0_15S_Month, LWCF0_15Dif_Month] =HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,0,14.5,'LWCF');
[LWCF15_30N_Month, LWCF15_30S_Month, LWCF15_30Dif_Month] = HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,14.5,30,'LWCF');
[LWCF30_49N_Month, LWCF30_49S_Month, LWCF30_49Dif_Month] = HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,30,48.5,'LWCF');
[LWCF49_90N_Month, LWCF49_90S_Month, LWCF49_90Dif_Month] = HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,48.5,90,'LWCF');

[SWclearNHMonthAllLats, SWclearSHMonthAllLats, SWclearNHMinusSHMonthAllLats] = HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,0,90,'SWclear')
[SWclear0_15N_Month, SWclear0_15S_Month, SWclear0_15Dif_Month] =HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,0,14.5,'SWclear');
[SWclear15_30N_Month, SWclear15_30S_Month, SWclear15_30Dif_Month] = HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,14.5,30,'SWclear');
[SWclear30_49N_Month, SWclear30_49S_Month, SWclear30_49Dif_Month] = HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,30,48.5,'SWclear');
[SWclear49_90N_Month, SWclear49_90S_Month, SWclear49_90Dif_Month] = HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,48.5,90,'SWclear');

[LWclearNHMonthAllLats, LWclearSHMonthAllLats, LWclearNHMinusSHMonthAllLats] = HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,0,90,'LWclear')
[LWclear0_15N_Month, LWclear0_15S_Month, LWclear0_15Dif_Month] =HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,0,14.5,'LWclear');
[LWclear15_30N_Month, LWclear15_30S_Month, LWclear15_30Dif_Month] = HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,14.5,30,'LWclear');
[LWclear30_49N_Month, LWclear30_49S_Month, LWclear30_49Dif_Month] = HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,30,48.5,'LWclear');
[LWclear49_90N_Month, LWclear49_90S_Month, LWclear49_90Dif_Month] = HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,48.5,90,'LWclear');

[TotalCloudForcingNHMonthAllLats, TotalCloudForcingSHMonthAllLats, TotalCloudForcingNHMinusSHMonthAllLats] = HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,0,90,'TotalCloudForcing')
[TotalCloudForcing0_15N_Month, TotalCloudForcing0_15S_Month, TotalCloudForcing0_15Dif_Month] =HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,0,14.5,'TotalCloudForcing');
[TotalCloudForcing15_30N_Month, TotalCloudForcing15_30S_Month, TotalCloudForcing15_30Dif_Month] = HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,14.5,30,'TotalCloudForcing');
[TotalCloudForcing30_49N_Month, TotalCloudForcing30_49S_Month, TotalCloudForcing30_49Dif_Month] = HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,30,48.5,'TotalCloudForcing');
[TotalCloudForcing49_90N_Month, TotalCloudForcing49_90S_Month, TotalCloudForcing49_90Dif_Month] = HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,48.5,90,'TotalCloudForcing');

load LatitudinalDifFluxTimeSeries.mat

%FluxTimeSeries.net= struct('Lat49_90',{},'Lat30_49',{})
FluxTimeSeries.net.Lat49_90 = struct('N',Net49_90N_Month,'S',Net49_90S_Month,'Dif',Net49_90Dif_Month);
FluxTimeSeries.net.Lat30_49 = struct('N',Net30_49N_Month,'S',Net30_49S_Month,'Dif',Net30_49Dif_Month);
FluxTimeSeries.net.Lat15_30 = struct('N',Net15_30N_Month,'S',Net15_30S_Month,'Dif',Net15_30Dif_Month);
FluxTimeSeries.net.Lat0_15 = struct('N',Net0_15N_Month,'S',Net0_15S_Month,'Dif',Net0_15Dif_Month);
FluxTimeSeries.net.AllLats = struct('N',NetNHMonthAllLats,'S',NetSHMonthAllLats,'Dif',NetNHMinusSHMonthAllLats);
FluxTimeSeries.LW.Lat49_90 = struct('N',LW49_90N_Month,'S',LW49_90S_Month,'Dif',LW49_90Dif_Month);
FluxTimeSeries.LW.Lat30_49 = struct('N',LW30_49N_Month,'S',LW30_49S_Month,'Dif',LW30_49Dif_Month);
FluxTimeSeries.LW.Lat15_30 = struct('N',LW15_30N_Month,'S',LW15_30S_Month,'Dif',LW15_30Dif_Month);
FluxTimeSeries.LW.Lat0_15 = struct('N',LW0_15N_Month,'S',LW0_15S_Month,'Dif',LW0_15Dif_Month);
FluxTimeSeries.LW.AllLats = struct('N',LWNHMonthAllLats,'S',LWSHMonthAllLats,'Dif',LWNHMinusSHMonthAllLats);
FluxTimeSeries.SW.Lat49_90 = struct('N',SW49_90N_Month,'S',SW49_90S_Month,'Dif',SW49_90Dif_Month);
FluxTimeSeries.SW.Lat30_49 = struct('N',SW30_49N_Month,'S',SW30_49S_Month,'Dif',SW30_49Dif_Month);
FluxTimeSeries.SW.Lat15_30 = struct('N',SW15_30N_Month,'S',SW15_30S_Month,'Dif',SW15_30Dif_Month);
FluxTimeSeries.SW.Lat0_15 = struct('N',SW0_15N_Month,'S',SW0_15S_Month,'Dif',SW0_15Dif_Month);
FluxTimeSeries.SW.AllLats = struct('N',SWNHMonthAllLats,'S',SWSHMonthAllLats,'Dif',SWNHMinusSHMonthAllLats);
FluxTimeSeries.LWCF.Lat49_90 = struct('N',LWCF49_90N_Month,'S',LWCF49_90S_Month,'Dif',LWCF49_90Dif_Month);
FluxTimeSeries.LWCF.Lat30_49 = struct('N',LWCF30_49N_Month,'S',LWCF30_49S_Month,'Dif',LWCF30_49Dif_Month);
FluxTimeSeries.LWCF.Lat15_30 = struct('N',LWCF15_30N_Month,'S',LWCF15_30S_Month,'Dif',LWCF15_30Dif_Month);
FluxTimeSeries.LWCF.Lat0_15 = struct('N',LWCF0_15N_Month,'S',LWCF0_15S_Month,'Dif',LWCF0_15Dif_Month);
FluxTimeSeries.LWCF.AllLats = struct('N',LWCFNHMonthAllLats,'S',LWCFSHMonthAllLats,'Dif',LWCFNHMinusSHMonthAllLats);
FluxTimeSeries.SWCF.Lat49_90 = struct('N',SWCF49_90N_Month,'S',SWCF49_90S_Month,'Dif',SWCF49_90Dif_Month);
FluxTimeSeries.SWCF.Lat30_49 = struct('N',SWCF30_49N_Month,'S',SWCF30_49S_Month,'Dif',SWCF30_49Dif_Month);
FluxTimeSeries.SWCF.Lat15_30 = struct('N',SWCF15_30N_Month,'S',SWCF15_30S_Month,'Dif',SWCF15_30Dif_Month);
FluxTimeSeries.SWCF.Lat0_15 = struct('N',SWCF0_15N_Month,'S',SWCF0_15S_Month,'Dif',SWCF0_15Dif_Month);
FluxTimeSeries.SWCF.AllLats = struct('N',SWCFNHMonthAllLats,'S',SWCFSHMonthAllLats,'Dif',SWCFNHMinusSHMonthAllLats);
FluxTimeSeries.SWclear.Lat49_90 = struct('N',SWclear49_90N_Month,'S',SWclear49_90S_Month,'Dif',SWclear49_90Dif_Month);
FluxTimeSeries.SWclear.Lat30_49 = struct('N',SWclear30_49N_Month,'S',SWclear30_49S_Month,'Dif',SWclear30_49Dif_Month);
FluxTimeSeries.SWclear.Lat15_30 = struct('N',SWclear15_30N_Month,'S',SWclear15_30S_Month,'Dif',SWclear15_30Dif_Month);
FluxTimeSeries.SWclear.Lat0_15 = struct('N',SWclear0_15N_Month,'S',SWclear0_15S_Month,'Dif',SWclear0_15Dif_Month);
FluxTimeSeries.SWclear.AllLats = struct('N',SWclearNHMonthAllLats,'S',SWclearSHMonthAllLats,'Dif',SWclearNHMinusSHMonthAllLats);
FluxTimeSeries.LWclear.Lat49_90 = struct('N',LWclear49_90N_Month,'S',LWclear49_90S_Month,'Dif',LWclear49_90Dif_Month);
FluxTimeSeries.LWclear.Lat30_49 = struct('N',LWclear30_49N_Month,'S',LWclear30_49S_Month,'Dif',LWclear30_49Dif_Month);
FluxTimeSeries.LWclear.Lat15_30 = struct('N',LWclear15_30N_Month,'S',LWclear15_30S_Month,'Dif',LWclear15_30Dif_Month);
FluxTimeSeries.LWclear.Lat0_15 = struct('N',LWclear0_15N_Month,'S',LWclear0_15S_Month,'Dif',LWclear0_15Dif_Month);
FluxTimeSeries.LWclear.AllLats = struct('N',LWclearNHMonthAllLats,'S',LWclearSHMonthAllLats,'Dif',LWclearNHMinusSHMonthAllLats);
FluxTimeSeries.TotalCloudForcing.Lat49_90 = struct('N',TotalCloudForcing49_90N_Month,'S',TotalCloudForcing49_90S_Month,'Dif',TotalCloudForcing49_90Dif_Month);
FluxTimeSeries.TotalCloudForcing.Lat30_49 = struct('N',TotalCloudForcing30_49N_Month,'S',TotalCloudForcing30_49S_Month,'Dif',TotalCloudForcing30_49Dif_Month);
FluxTimeSeries.TotalCloudForcing.Lat15_30 = struct('N',TotalCloudForcing15_30N_Month,'S',TotalCloudForcing15_30S_Month,'Dif',TotalCloudForcing15_30Dif_Month);
FluxTimeSeries.TotalCloudForcing.Lat0_15 = struct('N',TotalCloudForcing0_15N_Month,'S',TotalCloudForcing0_15S_Month,'Dif',TotalCloudForcing0_15Dif_Month);
FluxTimeSeries.TotalCloudForcing.AllLats = struct('N',TotalCloudForcingNHMonthAllLats,'S',TotalCloudForcingSHMonthAllLats,'Dif',TotalCloudForcingNHMinusSHMonthAllLats);
%FluxTimeSeries.net.Lat49_90.N
FluxTimeSeries.(fields{1}).Lat49_90.Dif

%orderfields to order them all


%accuracy check - should be close to 0
%max(max(mean(LWClimatologySubtracted,3)))
%mean(max(mean(LWClimatologySubtracted,3)))

LWSWCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*SWClimatologySubtracted,3),(std(LWClimatologySubtracted,0,3).*std(LWClimatologySubtracted,0,3)));
LWSWCorrMap = bsxfun(@rdivide,155^2*sum(LWClimatologySubtracted.*SWClimatologySubtracted,3),(std(LWClimatologySubtracted,0,3).*std(LWClimatologySubtracted,0,3)));
LWSWCorrMap = bsxfun(@rdivide,mean(LWClimatologySubtracted.*SWClimatologySubtracted,3),(std(LWClimatologySubtracted,0,3).*std(LWClimatologySubtracted,0,3)));
%LWSWCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*SWClimatologySubtracted,3),(moment(LWClimatologySubtracted,2,3).*moment(LWClimatologySubtracted,2,3)));


LWSWCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*SWClimatologySubtracted,3),(sqrt(sum(LWClimatologySubtracted.^2,3)).*sqrt(sum(SWClimatologySubtracted.^2,3))));
LWNetCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*NetClimatologySubtracted,3),(sqrt(sum(LWClimatologySubtracted.^2,3)).*sqrt(sum(NetClimatologySubtracted.^2,3))));
SWNetCorrMap = bsxfun(@rdivide,sum(SWClimatologySubtracted.*NetClimatologySubtracted,3),(sqrt(sum(SWClimatologySubtracted.^2,3)).*sqrt(sum(NetClimatologySubtracted.^2,3))));
LW_LWClearCorrMap = bsxfun(@rdivide,sum(LWclearClimatologySubtracted.*LWClimatologySubtracted,3),(sqrt(sum(LWclearClimatologySubtracted.^2,3)).*sqrt(sum(LWClimatologySubtracted.^2,3))));
LW_LWCFCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*LWCFClimatologySubtracted,3),(sqrt(sum(LWCFClimatologySubtracted.^2,3)).*sqrt(sum(LWClimatologySubtracted.^2,3))));
LWclear_SWclearCorrMap = bsxfun(@rdivide,sum(LWclearClimatologySubtracted.*SWclearClimatologySubtracted,3),(sqrt(sum(LWclearClimatologySubtracted.^2,3)).*sqrt(sum(SWclearClimatologySubtracted.^2,3))));

%9*8/2 = 36 possible pairwise comparisons

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(LWSWCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(LWNetCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(SWNetCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(LW_LWClearCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(LW_LWCFCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(LWclear_SWclearCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar

GlobalCorrMap(LW,SW)
GlobalCorrMap(LWClimatologySubtracted,SWClimatologySubtracted)
GlobalCorrMap(NetClimatologySubtracted,SWClimatologySubtracted)
GlobalCorrMap(LWclearClimatologySubtracted,SWClimatologySubtracted)
GlobalCorrMap(SWclearClimatologySubtracted,LWClimatologySubtracted);
GlobalCorrMap(SWCFClimatologySubtracted,LWCFClimatologySubtracted);
GlobalCorrMap(SWCFClimatologySubtracted,LWClimatologySubtracted);
GlobalCorrMap(NetClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(SWclearClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(LWclearClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(SWclearClimatologySubtracted,SWClimatologySubtracted);
GlobalCorrMap(LWClimatologySubtracted,LWCFClimatologySubtracted);
GlobalCorrMap(SWCFClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(LWCFClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(NetClimatologySubtracted,netclearClimatologySubtracted);
GlobalCorrMap(TotalCloudForcingClimatologySubtracted,netclearClimatologySubtracted);

GlobalRegressionMap(LWClimatologySubtracted,SWClimatologySubtracted)
GlobalRegressionMap(LWclearClimatologySubtracted,SWClimatologySubtracted)
GlobalRegressionMap(SWclearClimatologySubtracted,LWClimatologySubtracted);
GlobalRegressionMap(SWCFClimatologySubtracted,LWCFClimatologySubtracted);
GlobalRegressionMap(SWCFClimatologySubtracted,LWClimatologySubtracted);
GlobalRegressionMap(NetClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalRegressionMap(SWclearClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalRegressionMap(LWclearClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalRegressionMap(SWclearClimatologySubtracted,SWClimatologySubtracted);
GlobalRegressionMap(LWClimatologySubtracted,LWCFClimatologySubtracted);
GlobalRegressionMap(SWCFClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalRegressionMap(LWCFClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalRegressionMap(NetClimatologySubtracted,netclearClimatologySubtracted);
GlobalRegressionMap(TotalCloudForcingClimatologySubtracted,netclearClimatologySubtracted);

GlobalCovMap(NetClimatologySubtracted,LWClimatologySubtracted);
GlobalCovMap(NetClimatologySubtracted,SWClimatologySubtracted);
GlobalCovMap(LWClimatologySubtracted,SWClimatologySubtracted)
GlobalCovMap(LWclearClimatologySubtracted,SWClimatologySubtracted)
GlobalCovMap(SWclearClimatologySubtracted,LWClimatologySubtracted);
GlobalCovMap(SWCFClimatologySubtracted,LWCFClimatologySubtracted);
GlobalCovMap(SWCFClimatologySubtracted,LWClimatologySubtracted);
GlobalCovMap(NetClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCovMap(SWclearClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCovMap(LWclearClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCovMap(SWclearClimatologySubtracted,SWClimatologySubtracted);
GlobalCovMap(LWClimatologySubtracted,LWCFClimatologySubtracted);
GlobalCovMap(SWCFClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCovMap(LWCFClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCovMap(NetClimatologySubtracted,netclearClimatologySubtracted);
GlobalCovMap(TotalCloudForcingClimatologySubtracted,netclearClimatologySubtracted);

Corr1DTimeSeriesMap(NetClimatologySubtracted,NINO12);
Corr1DTimeSeriesMap(NetClimatologySubtracted,NINO3);
Corr1DTimeSeriesMap(NetClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(NetClimatologySubtracted,NINO4);
Corr1DTimeSeriesMap(NetClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(NetClimatologySubtracted,T);
Corr1DTimeSeriesMap(SWClimatologySubtracted,NINO12);
Corr1DTimeSeriesMap(SWClimatologySubtracted,NINO3);
Corr1DTimeSeriesMap(SWClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(SWClimatologySubtracted,NINO4);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(LWClimatologySubtracted,NINO12);
Corr1DTimeSeriesMap(LWClimatologySubtracted,NINO3);
Corr1DTimeSeriesMap(LWClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(LWClimatologySubtracted,NINO4);
Corr1DTimeSeriesMap(LWClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,NINO12);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,NINO3);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,NINO4);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,T);
Corr1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NINO12);
Corr1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NINO3);
Corr1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NINO4);
Corr1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(SWclearClimatologySubtracted,NINO12);
Corr1DTimeSeriesMap(SWclearClimatologySubtracted,NINO3);
Corr1DTimeSeriesMap(SWclearClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(SWclearClimatologySubtracted,NINO4);
Corr1DTimeSeriesMap(SWclearClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,NINO12);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,NINO3);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,NINO4);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,T);
Corr1DTimeSeriesMap(LWCFClimatologySubtracted,NINO12);
Corr1DTimeSeriesMap(LWCFClimatologySubtracted,NINO3);
Corr1DTimeSeriesMap(LWCFClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(LWCFClimatologySubtracted,NINO4);
Corr1DTimeSeriesMap(LWCFClimatologySubtracted,SAO);

Regress1DTimeSeriesMap(NetClimatologySubtracted,NINO12);
Regress1DTimeSeriesMap(NetClimatologySubtracted,NINO3);
Regress1DTimeSeriesMap(NetClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(NetClimatologySubtracted,NINO4);
Regress1DTimeSeriesMap(NetClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(NetClimatologySubtracted,T);
Regress1DTimeSeriesMap(SWClimatologySubtracted,NINO12);
Regress1DTimeSeriesMap(SWClimatologySubtracted,NINO3);
Regress1DTimeSeriesMap(SWClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(SWClimatologySubtracted,NINO4);
Regress1DTimeSeriesMap(SWClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(SWClimatologySubtracted,T);
Regress1DTimeSeriesMap(LWClimatologySubtracted,NINO12);
Regress1DTimeSeriesMap(LWClimatologySubtracted,NINO3);
Regress1DTimeSeriesMap(LWClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(LWClimatologySubtracted,NINO4);
Regress1DTimeSeriesMap(LWClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,NINO12);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,NINO3);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,NINO4);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,T);
Regress1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NINO12);
Regress1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NINO3);
Regress1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NINO4);
Regress1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(SWclearClimatologySubtracted,NINO12);
Regress1DTimeSeriesMap(SWclearClimatologySubtracted,NINO3);
Regress1DTimeSeriesMap(SWclearClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(SWclearClimatologySubtracted,NINO4);
Regress1DTimeSeriesMap(SWclearClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,NINO12);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,NINO3);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,NINO4);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,T);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,NINO12);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,NINO3);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,NINO4);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,T);

Corr1DTimeSeriesMap(NetClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(NetClimatologySubtracted,NAO);

subplot(2,2,1)
Corr1DTimeSeriesMap(NetClimatologySubtracted(:,:,WinterIndices),NAO(WinterIndices));
subplot(2,2,2)
Corr1DTimeSeriesMap(NetClimatologySubtracted(:,:,SpringIndices),NAO(SpringIndices));
subplot(2,2,3)
Corr1DTimeSeriesMap(NetClimatologySubtracted(:,:,SummerIndices),NAO(SummerIndices));
subplot(2,2,4)
Corr1DTimeSeriesMap(NetClimatologySubtracted(:,:,AutumnIndices),NAO(AutumnIndices));
title(subplot(2,2,1),'Net-NAO Winter Crrln')
title(subplot(2,2,2),'Net-NAO Spring Crrln')
title(subplot(2,2,3),'Net-NAO Summer Crrln')
title(subplot(2,2,4),'Net-NAO Autumn Crrln')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['AllSeasonTimeSeriesCorr_','NetClimatologySubtracted','-','NAO','.png']);

PlotSeasonalCorrRegMaps1DTimeSeries(NetClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(LWClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(SWClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(LWclearClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(SWclearClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(LWCFClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(SWCFClimatologySubtracted,NAO)

PlotSeasonalCorrRegMaps1DTimeSeries(NetClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(LWClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(SWClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(LWclearClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(SWclearClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(LWCFClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(SWCFClimatologySubtracted,NAM)

Corr1DTimeSeriesMap(NetClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(NetClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(NetClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(NetClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(NetClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(NetClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(SWClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(SWClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(SWClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(SWClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(SWClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(SWClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(SWClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(SWClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(LWClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(LWClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(LWClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(LWClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(LWClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(LWClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(LWClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(LWClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(TotalCloudForcingClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(LWclearClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(LWclearClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(SWclearClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(SWclearClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(SWclearClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(SWclearClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(SWclearClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(SWclearClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(SWclearClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(SWclearClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(SWCFClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(SWCFClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(LWCFClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(LWCFClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(LWCFClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(LWCFClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(LWCFClimatologySubtracted,PNA);


% LinearModel.fit(squeeze(LWclearClimatologySubtracted(1,1,:)),NINO12)
% 
% GlobalCorrMap(LWClimatologySubtracted,LW)
% GlobalCorrMap(LW,LW);


corr(NAM,NetNHMinusSHMonthAllLats')
corr(SAM,NetNHMinusSHMonthAllLats')
corr(NAO,NetNHMinusSHMonthAllLats')
corr(NINO34,NetNHMinusSHMonthAllLats')

corr(NAM,Net0_15Dif_Month')
corr(SAM,Net0_15Dif_Month')%-0.15706
corr(NAO,Net0_15Dif_Month')
corr(NINO34,Net0_15Dif_Month')%0.4631

corr(NAM,Net15_30Dif_Month')
corr(SAM,Net15_30Dif_Month')
corr(NAO,Net15_30Dif_Month')
corr(NINO34,Net15_30Dif_Month')

corr(NAM,Net30_49Dif_Month')
corr(SAM,Net30_49Dif_Month')% 0.31789
corr(NAO,Net30_49Dif_Month')
corr(NINO34,Net30_49Dif_Month')% -0.19975

corr(NAM,Net49_90Dif_Month')
corr(SAM,Net49_90Dif_Month')%-0.17849
corr(NAO,Net49_90Dif_Month')
corr(NAO(WinterIndices),Net49_90Dif_Month(WinterIndices)') %.189
corr(NINO34,Net49_90Dif_Month')

corr(NAM,SWNHMinusSHMonthAllLats')%-0.10823
corr(SAM,SWNHMinusSHMonthAllLats')%0.20108
corr(NAO,SWNHMinusSHMonthAllLats')
corr(NINO34,SWNHMinusSHMonthAllLats')

corr(NAM,SW0_15Dif_Month')%-0.19588
corr(SAM,SW0_15Dif_Month')%-0.14271
corr(NAO,SW0_15Dif_Month')
corr(NINO34,SW0_15Dif_Month')%0.27466

corr(NAM,SW15_30Dif_Month')% -0.21443
corr(SAM,SW15_30Dif_Month')%  0.17455
corr(NAO,SW15_30Dif_Month')
corr(NINO34,SW15_30Dif_Month')

corr(NAM,SW30_49Dif_Month')%0.21936
corr(SAM,SW30_49Dif_Month')%0.27816
corr(NAO,SW30_49Dif_Month')%0.13408
corr(NINO34,SW30_49Dif_Month')%-0.20794

corr(NAM,SW49_90Dif_Month')
corr(SAM,SW49_90Dif_Month')%0.15879
corr(NAO,SW49_90Dif_Month')%-0.19248
corr(NAO(WinterIndices),SW49_90Dif_Month(WinterIndices)')%-0.228
corr(NINO34,SW49_90Dif_Month')%-0.19958

corr(NAM,LWNHMinusSHMonthAllLats')
corr(SAM,LWNHMinusSHMonthAllLats')%-0.18163
corr(NAO,LWNHMinusSHMonthAllLats')
corr(NINO34,LWNHMinusSHMonthAllLats')

corr(NAM,LW0_15Dif_Month')% 0.14101
corr(SAM,LW0_15Dif_Month')
corr(NAO,LW0_15Dif_Month')
corr(NINO34,LW0_15Dif_Month')

corr(NAM,LW15_30Dif_Month')% 0.16164
corr(SAM,LW15_30Dif_Month')
corr(NAO,LW15_30Dif_Month')
corr(NINO34,LW15_30Dif_Month')

corr(NAM,LW30_49Dif_Month')%-0.36864
corr(SAM,LW30_49Dif_Month')
corr(NAO,LW30_49Dif_Month')%-0.15669
corr(NINO34,LW30_49Dif_Month')

corr(NAM,LW49_90Dif_Month')
corr(SAM,LW49_90Dif_Month') %-.372
corr(NAO,LW49_90Dif_Month') %.205
corr(NINO34,LW49_90Dif_Month') %.152
%%%%%%%%%
corr(NAM,LWclearNHMinusSHMonthAllLats')
corr(SAM,LWclearNHMinusSHMonthAllLats') %-.22
corr(NAO,LWclearNHMinusSHMonthAllLats')
corr(NINO34,LWclearNHMinusSHMonthAllLats') %.14

corr(NAM,LWclear0_15Dif_Month')%.168
corr(SAM,LWclear0_15Dif_Month')%.11
corr(NAO,LWclear0_15Dif_Month')
corr(NINO34,LWclear0_15Dif_Month')%

corr(NAM,LWclear15_30Dif_Month')%.165
corr(SAM,LWclear15_30Dif_Month')
corr(NAO,LWclear15_30Dif_Month')
corr(NINO34,LWclear15_30Dif_Month')

corr(NAM,LWclear30_49Dif_Month')%-.225
corr(SAM,LWclear30_49Dif_Month')% 
corr(NAO,LWclear30_49Dif_Month')
corr(NINO34,LWclear30_49Dif_Month')% 

corr(NAM,LWclear49_90Dif_Month')%-.16
corr(SAM,LWclear49_90Dif_Month')%-0.33214
corr(NAO,LWclear49_90Dif_Month')
corr(NINO34,LWclear49_90Dif_Month') %.17
%%
corr(NAM,LWCFNHMinusSHMonthAllLats') %.18
corr(SAM,LWCFNHMinusSHMonthAllLats')
corr(NAO,LWCFNHMinusSHMonthAllLats')
corr(NINO34,LWCFNHMinusSHMonthAllLats')

corr(NAM,LWCF0_15Dif_Month')
corr(SAM,LWCF0_15Dif_Month')
corr(NAO,LWCF0_15Dif_Month')
corr(NINO34,LWCF0_15Dif_Month')

corr(NAM,LWCF15_30Dif_Month')
corr(SAM,LWCF15_30Dif_Month')
corr(NAO,LWCF15_30Dif_Month')
corr(NINO34,LWCF15_30Dif_Month')

corr(NAM,LWCF30_49Dif_Month')%-0.35308
corr(SAM,LWCF30_49Dif_Month')
corr(NAO,LWCF30_49Dif_Month')%-0.18171
corr(NINO34,LWCF30_49Dif_Month')% 0.08036

corr(NAM,LWCF49_90Dif_Month')%0.37794
corr(SAM,LWCF49_90Dif_Month')%-0.15386
corr(NAO,LWCF49_90Dif_Month')%0.27265
corr(NINO34,LWCF49_90Dif_Month')
%%

corr(NAM,SWCFNHMinusSHMonthAllLats')
corr(SAM,SWCFNHMinusSHMonthAllLats') %.112
corr(NAO,SWCFNHMinusSHMonthAllLats')
corr(NINO34,SWCFNHMinusSHMonthAllLats') %.15

corr(NAM,SWCF0_15Dif_Month')%-0.16869
corr(SAM,SWCF0_15Dif_Month')%-0.13631
corr(NAO,SWCF0_15Dif_Month')
corr(NINO34,SWCF0_15Dif_Month')%0.2972

corr(NAM,SWCF15_30Dif_Month')%-0.16364
corr(SAM,SWCF15_30Dif_Month')%0.1758
corr(NAO,SWCF15_30Dif_Month')
corr(NINO34,SWCF15_30Dif_Month')

corr(NAM,SWCF30_49Dif_Month')%0.29425
corr(SAM,SWCF30_49Dif_Month')% 0.2868
corr(NAO,SWCF30_49Dif_Month')%0.1543
corr(NINO34,SWCF30_49Dif_Month')% -0.19975

corr(NAM,SWCF49_90Dif_Month')
corr(SAM,SWCF49_90Dif_Month')
corr(NAO,SWCF49_90Dif_Month')
corr(NINO34,SWCF49_90Dif_Month')% 0.1836
%%
corr(NAM,SWclearNHMinusSHMonthAllLats') %-.10
corr(SAM,SWclearNHMinusSHMonthAllLats') %.122
corr(NAO,SWclearNHMinusSHMonthAllLats')%-.118
corr(NINO34,SWclearNHMinusSHMonthAllLats') %-.25

corr(NAM,SWclear0_15Dif_Month')%-0.18936
corr(SAM,SWclear0_15Dif_Month')%
corr(NAO,SWclear0_15Dif_Month')
corr(NINO34,SWclear0_15Dif_Month')%-.15

corr(NAM,SWclear15_30Dif_Month')%-.17
corr(SAM,SWclear15_30Dif_Month')
corr(NAO,SWclear15_30Dif_Month')
corr(NINO34,SWclear15_30Dif_Month')

corr(NAM,SWclear30_49Dif_Month')
corr(SAM,SWclear30_49Dif_Month')% 
corr(NAO,SWclear30_49Dif_Month')
corr(NAO(WinterIndices),SWclear30_49Dif_Month(WinterIndices)')
corr(NINO34,SWclear30_49Dif_Month')% 

corr(NAM,SWclear49_90Dif_Month')
corr(SAM,SWclear49_90Dif_Month')%.18
corr(NAO,SWclear49_90Dif_Month')%-.17
corr(NINO34,SWclear49_90Dif_Month')%-0.30672
%%
corr(NAM,TotalCloudForcingNHMinusSHMonthAllLats')
corr(SAM,TotalCloudForcingNHMinusSHMonthAllLats')
corr(NAO,TotalCloudForcingNHMinusSHMonthAllLats') %.14
corr(NINO34,TotalCloudForcingNHMinusSHMonthAllLats') %.255

corr(NAM,netclearNHMinusSHMonthAllLats') %-.13
corr(SAM,netclearNHMinusSHMonthAllLats')
corr(NAO,netclearNHMinusSHMonthAllLats')
corr(NINO34,netclearNHMinusSHMonthAllLats') %-.17
corr(NAM,netclear0_15Dif_Month')
corr(SAM,netclear0_15Dif_Month')%
corr(NAO,netclear0_15Dif_Month')
corr(NINO34,netclear0_15Dif_Month')%

corr(NAM,netclear15_30Dif_Month')
corr(SAM,netclear15_30Dif_Month')
corr(NAO,netclear15_30Dif_Month')
corr(NINO34,netclear15_30Dif_Month')

corr(NAM,netclear30_49Dif_Month')%-.25
corr(SAM,netclear30_49Dif_Month')%
corr(NAO,netclear30_49Dif_Month')
corr(NINO34,netclear30_49Dif_Month')% 

corr(NAM,netclear49_90Dif_Month')%_.16
corr(SAM,netclear49_90Dif_Month')%
corr(NAO,netclear49_90Dif_Month')%-.16
corr(NINO34,netclear49_90Dif_Month')%-.22


%%%%%%%%%%%%%
corr(LW49_90Dif_Month,SW49_90Dif_Month)
corr(LW49_90Dif_Month',SW49_90Dif_Month')
corr(LW49_90Dif_Month',Net49_90Dif_Month')

corr(NAM,SAM) %.18
corr(NAM,NAO) %.65
corr(NAM,NINO34) %-.15
corr(NAM,PNA)%-.17
corr(NAO,SAM)
corr(NAO,NINO34)%-.17
corr(NAO,PNA)
corr(NINO34,PNA)
corr(NINO34,SAM)%-.26
corr(SAM,PNA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trend data? ONLY VALID IF YOU DO PROPER AREA-WEIGHTING OF GRID CELLS
%but actually this was already done in NetNHMonthAllLats
mean(diff(NetNHMonthAllLats))
mean(diff(NetSHMonthAllLats))

%but in the north...

% mean(diff(Net49_90N_Month))
% 
% ans =
% 
%       0.01099
% 
% 0.01099*156
% 
% ans =
% 
%        1.7144

% mean(diff(LWNHMonthAllLats))
% 
% ans =
% 
%      0.011813
% 
% mean(diff(SWNHMonthAllLats))
% 
% ans =
% 
%     -0.016836
% 

%%%%%

%average slope of each variable over time for each grid cell

load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(mean(diff(NetClimatologySubtracted,1,3),3),geoidrefvec,'DisplayType','texturemap');colorbar

load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(mean(diff(LWclearClimatologySubtracted,1,3),3),geoidrefvec,'DisplayType','texturemap');colorbar

% corr(NINO34,T)
% 
% ans =
% 
%       0.41754
% 
% corr(NAM,T)
% 
% ans =
% 
%         0.129
% 
% corr(SAM,T)
% 
% ans =
% 
%     -0.028051
% 
% corr(NAO,T)
% 
% ans =
% 
%      -0.10703

