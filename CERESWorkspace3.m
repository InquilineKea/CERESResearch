clear all

precip = ncread('precip.mon.mean.nc','precip');
lat = ncread('precip.mon.mean.nc','lat');
long = ncread('precip.mon.mean.nc','lon');
preciptime = ncread('precip.mon.mean.nc','time');
precip = permute(precip,[2 1 3]);
precip = precip(:,:,find(preciptime== 73108):find(preciptime== 77828));
temp = ncread('t2m.nc','t2m');
temptime = ncread('t2m.nc','time');
tempLats= ncread('t2m.nc','latitude');
tempLongs = ncread('t2m.nc','longitude');
temp = temp(:,:,find(temptime==878016):find(temptime==991296));
temp = permute(temp,[2 1 3]);
seaice = ncread('seaice.nc','ci');
seaice = seaice(:,:,find(temptime==878016):find(temptime==991296));
seaice = permute(seaice,[2 1 3]);
albedo= ncread('albedo.nc','al');
albedo = albedo(:,:,find(temptime==878016):find(temptime==991296));
albedo = permute(albedo,[2 1 3]);
EastHeatFlux = ncread('albedo.nc','p69.162');
EastHeatFlux = EastHeatFlux(:,:,find(temptime==878016):find(temptime==991296));
EastHeatFlux = permute(EastHeatFlux,[2 1 3]);
NorthHeatFlux = ncread('albedo.nc','p70.162');
NorthHeatFlux = NorthHeatFlux(:,:,find(temptime==878016):find(temptime==991296));
NorthHeatFlux = permute(NorthHeatFlux,[2 1 3]);


for depth=1:size(precip,3)
  precip(:,:,depth) = flipud(precip(:,:,depth));
end

E=zeros(180,360,156);
for depth=1:size(temp,3)
  temp(:,:,depth) = flipud(temp(:,:,depth));
  E(:,:,depth)=imresize(temp(:,:,depth),[180 360]);
end
temp = E;
clear E;
GlobePlot(temp,1);

E=zeros(180,360,156);
for depth=1:size(seaice,3)
  seaice(:,:,depth) = flipud(seaice(:,:,depth));
  E(:,:,depth)=imresize(seaice(:,:,depth),[180 360]);
end
seaice = E;
clear E;

E=zeros(180,360,156);
for depth=1:size(albedo,3)
  albedo(:,:,depth) = flipud(albedo(:,:,depth));
  E(:,:,depth)=imresize(albedo(:,:,depth),[180 360]);
end
albedo = E;
clear E;
E=zeros(180,360,156);
for depth=1:size(EastHeatFlux,3)
  EastHeatFlux(:,:,depth) = flipud(EastHeatFlux(:,:,depth));
  E(:,:,depth)=imresize(EastHeatFlux(:,:,depth),[180 360]);
end
EastHeatFlux= E;
clear E;

for depth=1:size(NorthHeatFlux,3)
  NorthHeatFlux(:,:,depth) = flipud(NorthHeatFlux(:,:,depth));
  E(:,:,depth)=imresize(NorthHeatFlux(:,:,depth),[180 360]);
end
NorthHeatFlux= E;
clear E;


% [NewX,NewY,NewZ] = meshgrid(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lon'),ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat'),find(temptime==878016):find(temptime==991296));
% [OldX,OldY,OldZ] = meshgrid(double(ncread('t2m.nc','longitude')),double(ncread('t2m.nc','latitude')),find(temptime==878016):find(temptime==991296));
% temp = interp3(OldX,OldY,OldZ,temp,NewX,NewY,NewZ,'linear',mean(temp(:,end-1,1)));
TempClimatologySubtracted = GlobalValuesMinusClimatology(temp,'temp');
PrecipClimatologySubtracted = GlobalValuesMinusClimatology(precip,'precip');
SeaIceClimatologySubtracted = GlobalValuesMinusClimatology(seaice,'seaice');
AlbedoClimatologySubtracted = GlobalValuesMinusClimatology(albedo,'albedo');
EastHeatFluxClimatologySubtracted= GlobalValuesMinusClimatology(EastHeatFlux,'EastHeatFlux');
NorthHeatFluxClimatologySubtracted= GlobalValuesMinusClimatology(NorthHeatFlux,'NorthHeatFlux');

save('ClimatologySubtractedVars.mat','PrecipClimatologySubtracted','-append')
save('ClimatologySubtractedVars.mat','AlbedoClimatologySubtracted','-append')
%save('ClimatologySubtractedVars.mat','SeaIceClimatologySubtracted','-append')
%save('ClimatologySubtractedVars.mat','TempClimatologySubtracted','-append')
%load('ClimatologySubtractedVars.mat','PrecipClimatologySubtracted');

load geoid;
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(NewY(:,:,1),geoidrefvec,'DisplayType','texturemap');colorbar


load MonthlyDeparturesStructure.mat
MonthlyDeparturesStructure.Temp = TempClimatologySubtracted;
save MonthlyDeparturesStructure.mat

load('ClimatologySubtractedVars','TempClimatologySubtracted');
load('ClimatologySubtractedVars','NetClimatologySubtracted');
load('ClimatologySubtractedVars','LWCFClimatologySubtracted');
load('ClimatologySubtractedVars','SWCFClimatologySubtracted');
load('ClimatologySubtractedVars','LWclearClimatologySubtracted');
load('ClimatologySubtractedVars','SWclearClimatologySubtracted');
load('ClimatologySubtractedVars','LWClimatologySubtracted');
load('ClimatologySubtractedVars','SWClimatologySubtracted');

PlotIndexThreshold(TempClimatologySubtracted,NINO34,90,'Temperature','NINO34')
PlotIndexThresholdLower(TempClimatologySubtracted,NINO34,10,'Temperature','NINO34')

PlotIndexThreshold(NetClimatologySubtracted,NAM,90,'Net','NAM')
PlotIndexThreshold(NetClimatologySubtracted,SAM,90,'Net','SAM')
PlotIndexThreshold(NetClimatologySubtracted,NAO,90,'Net','NAO')
PlotIndexThresholdLower(NetClimatologySubtracted,NAM,10,'Temperature','NINO34')

PlotIndexThreshold(NetClimatologySubtracted,NAM,90,'Net','NAM')
PlotIndexThreshold(NetClimatologySubtracted,SAM,90,'Net','SAM')
PlotIndexThreshold(NetClimatologySubtracted,NAO,90,'Net','NAO')

PlotIndexThreshold(LWCFClimatologySubtracted,NAM,90,'LWCF','NAM')
PlotIndexThreshold(LWCFClimatologySubtracted,SAM,90,'LWCF','SAM')
PlotIndexThreshold(LWCFClimatologySubtracted,NAO,90,'LWCF','NAO')

PlotIndexThreshold(SWCFClimatologySubtracted,NAM,90,'SWCF','NAM')
PlotIndexThreshold(SWCFClimatologySubtracted,SAM,90,'SWCF','SAM')
PlotIndexThreshold(SWCFClimatologySubtracted,NAO,90,'SWCF','NAO')
PlotIndexThreshold(LWclearClimatologySubtracted,NAM,90,'LWclear','NAM')
PlotIndexThreshold(LWclearClimatologySubtracted,SAM,90,'LWclear','SAM')
PlotIndexThreshold(LWclearClimatologySubtracted,NAO,90,'LWclear','NAO')
PlotIndexThreshold(SWclearClimatologySubtracted,NAM,90,'SWclear','NAM')
PlotIndexThreshold(SWclearClimatologySubtracted,SAM,90,'SWclear','SAM')
PlotIndexThreshold(SWclearClimatologySubtracted,NAO,90,'SWclear','NAO')
PlotIndexThreshold(SWClimatologySubtracted,NAM,90,'SW','NAM')
PlotIndexThreshold(SWClimatologySubtracted,SAM,90,'SW','SAM')
PlotIndexThreshold(SWClimatologySubtracted,NAO,90,'SW','NAO')
PlotIndexThreshold(LWClimatologySubtracted,NAM,90,'LW','NAM')
PlotIndexThreshold(LWClimatologySubtracted,SAM,90,'LW','SAM')
PlotIndexThreshold(LWClimatologySubtracted,NAO,90,'LW','NAO')

PlotIndexThreshold(SWClimatologySubtracted,NAM.*NINO34,90,'SW','NAM_NINO34')
PlotIndexThreshold(SWClimatologySubtracted,NAM+NINO34,90,'SW','NAM_NINO34')

find(NINO34 < prctile(NINO34,90) & NINO34 > prctile(NINO34,10) & NAM > prctile(NAM, 90))

plot(NAO)
hold on;
plot(NINO34,'r')
plot(NAM,'g')
plot(SAM,'m')
corr(NAM,NINO34)
corr(NAO,NINO34)
corr(SAM,NINO34)

plot(NAO.*NINO34)

[B,BIN] = regress(NAM,NINO34)
[B,BIN] = regress(NINO34,NAM)
[B,BIN] =regress(NAO,NINO34)
[B,BIN] =regress(SAM,NINO34)
[B,BIN] =regress(NAO,NAM)
[B,BIN] =regress(NAM,NAO)
[B,BIN] =regress(SAM,NAM)
[B,BIN] =regress(SAM,NAO)

GlobalCorrMap(ClimateSubtractStruct.Temp,ClimateSubtractStruct.Precip)

% 
% load geoid
% set(gca,'FontSize',20)
% ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
% geoshow(resizem(precip(:,:,1)',2.5),geoidrefvec,'DisplayType','texturemap');colorbar
% 
% set(gca,'FontSize',20)
% ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
% geoshow(resizem(mean(precip,3)',2.5),geoidrefvec,'DisplayType','texturemap');colorbar

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

Indices.NAM = NAM';
Indices.SAM = SAM';
Indices.NAO = NAO';
Indices.NINO34 = NINO34';
Indices.PNA = PNA';
Indices.SAO = double(SAO');
Indices.T = T';

save('Indices.mat','Indices')

hold off;

Corr1DTimeSeriesMap(PrecipClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(PrecipClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(PrecipClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(PrecipClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(PrecipClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(PrecipClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(PrecipClimatologySubtracted,T);

Regress1DTimeSeriesMap(PrecipClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(PrecipClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(PrecipClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(PrecipClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(PrecipClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(PrecipClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(PrecipClimatologySubtracted,T);

%load ClimatologySubtractedVars.mat
%save('ClimatologySubtractedVars.mat','PrecipClimatologySubtracted','-append')

GlobalCorrMap(LWClimatologySubtracted,PrecipClimatologySubtracted)
GlobalCorrMap(LWclearClimatologySubtracted,PrecipClimatologySubtracted)
GlobalCorrMap(SWClimatologySubtracted,PrecipClimatologySubtracted)
GlobalCorrMap(SWclearClimatologySubtracted,PrecipClimatologySubtracted)
GlobalCorrMap(SWCFClimatologySubtracted,PrecipClimatologySubtracted)
GlobalCorrMap(LWCFClimatologySubtracted,PrecipClimatologySubtracted)
GlobalCorrMap(TotalCloudForcingClimatologySubtracted,PrecipClimatologySubtracted)
GlobalCorrMap(NetClimatologySubtracted,PrecipClimatologySubtracted)
GlobalCorrMap(netclearClimatologySubtracted,PrecipClimatologySubtracted)
GlobalCorrMap(TempClimatologySubtracted,PrecipClimatologySubtracted)


GlobalRegressionMap(PrecipClimatologySubtracted,LWClimatologySubtracted)
%GlobalRegressionMap(LWClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(PrecipClimatologySubtracted,LWclearClimatologySubtracted)
GlobalRegressionMap(PrecipClimatologySubtracted,SWClimatologySubtracted)
GlobalRegressionMap(PrecipClimatologySubtracted,SWCFClimatologySubtracted)
GlobalRegressionMap(PrecipClimatologySubtracted,LWCFClimatologySubtracted)
GlobalRegressionMap(PrecipClimatologySubtracted,TotalCloudForcingClimatologySubtracted)
GlobalRegressionMap(PrecipClimatologySubtracted,NetClimatologySubtracted)
GlobalRegressionMap(PrecipClimatologySubtracted,netclearClimatologySubtracted)
GlobalRegressionMap(PrecipClimatologySubtracted,SWclearClimatologySubtracted)
GlobalRegressionMap(PrecipClimatologySubtracted,TempClimatologySubtracted)

GlobalRegressionMap(LWClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(LWclearClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(SWClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(SWclearClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(SWCFClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(LWCFClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(TotalCloudForcingClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(NetClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(netclearClimatologySubtracted,PrecipClimatologySubtracted)
GlobalRegressionMap(TempClimatologySubtracted,PrecipClimatologySubtracted)


%%%%%%%%%%%%%
PlotSeasonalCorrRegMaps1DTimeSeries(PrecipClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(PrecipClimatologySubtracted,SAM)
PlotSeasonalCorrRegMaps1DTimeSeries(PrecipClimatologySubtracted,NINO34)
PlotSeasonalCorrRegMaps1DTimeSeries(PrecipClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(PrecipClimatologySubtracted,T)

%%%%

Corr1DTimeSeriesMap(TempClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(TempClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(TempClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(TempClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(TempClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(TempClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(TempClimatologySubtracted,T);

Regress1DTimeSeriesMap(TempClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(TempClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(TempClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(TempClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(TempClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(TempClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(TempClimatologySubtracted,T);

%load ClimatologySubtractedVars.mat
%save('ClimatologySubtractedVars.mat','TempClimatologySubtracted','-append')

GlobalCorrMap(LWClimatologySubtracted,TempClimatologySubtracted)
GlobalCorrMap(LWclearClimatologySubtracted,TempClimatologySubtracted)
GlobalCorrMap(SWClimatologySubtracted,TempClimatologySubtracted)
GlobalCorrMap(SWclearClimatologySubtracted,TempClimatologySubtracted)
GlobalCorrMap(SWCFClimatologySubtracted,TempClimatologySubtracted)
GlobalCorrMap(LWCFClimatologySubtracted,TempClimatologySubtracted)
GlobalCorrMap(TotalCloudForcingClimatologySubtracted,TempClimatologySubtracted)
GlobalCorrMap(NetClimatologySubtracted,TempClimatologySubtracted)
GlobalCorrMap(netclearClimatologySubtracted,TempClimatologySubtracted)

GlobalRegressionMap(TempClimatologySubtracted,LWClimatologySubtracted)
%GlobalRegressionMap(LWClimatologySubtracted,TempClimatologySubtracted)
GlobalRegressionMap(TempClimatologySubtracted,LWclearClimatologySubtracted)
GlobalRegressionMap(TempClimatologySubtracted,SWClimatologySubtracted)
GlobalRegressionMap(TempClimatologySubtracted,SWCFClimatologySubtracted)
GlobalRegressionMap(TempClimatologySubtracted,LWCFClimatologySubtracted)
GlobalRegressionMap(TempClimatologySubtracted,TotalCloudForcingClimatologySubtracted)
GlobalRegressionMap(TempClimatologySubtracted,NetClimatologySubtracted)
GlobalRegressionMap(TempClimatologySubtracted,netclearClimatologySubtracted)
GlobalRegressionMap(TempClimatologySubtracted,SWclearClimatologySubtracted)

GlobalRegressionMap(LWClimatologySubtracted,TempClimatologySubtracted)
GlobalRegressionMap(LWclearClimatologySubtracted,TempClimatologySubtracted)
GlobalRegressionMap(SWClimatologySubtracted,TempClimatologySubtracted)
GlobalRegressionMap(SWclearClimatologySubtracted,TempClimatologySubtracted)
GlobalRegressionMap(SWCFClimatologySubtracted,TempClimatologySubtracted)
GlobalRegressionMap(LWCFClimatologySubtracted,TempClimatologySubtracted)
GlobalRegressionMap(TotalCloudForcingClimatologySubtracted,TempClimatologySubtracted)
GlobalRegressionMap(NetClimatologySubtracted,TempClimatologySubtracted)
GlobalRegressionMap(netclearClimatologySubtracted,TempClimatologySubtracted)

%%%%%%%%%%%%%
PlotSeasonalCorrRegMaps1DTimeSeries(TempClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(TempClimatologySubtracted,SAM)
PlotSeasonalCorrRegMaps1DTimeSeries(TempClimatologySubtracted,NINO34)
PlotSeasonalCorrRegMaps1DTimeSeries(TempClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(TempClimatologySubtracted,T)

%%%%%

Corr1DTimeSeriesMap(SeaIceClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(SeaIceClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(SeaIceClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(SeaIceClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(SeaIceClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(SeaIceClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(SeaIceClimatologySubtracted,T);

Regress1DTimeSeriesMap(SeaIceClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(SeaIceClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(SeaIceClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(SeaIceClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(SeaIceClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(SeaIceClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(SeaIceClimatologySubtracted,T);

%load ClimatologySubtractedVars.mat
%save('ClimatologySubtractedVars.mat','SeaIceClimatologySubtracted','-append')

GlobalCorrMap(LWClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalCorrMap(LWclearClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalCorrMap(SWClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalCorrMap(SWclearClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalCorrMap(SWCFClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalCorrMap(LWCFClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalCorrMap(TotalCloudForcingClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalCorrMap(NetClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalCorrMap(netclearClimatologySubtracted,SeaIceClimatologySubtracted)

GlobalRegressionMap(SeaIceClimatologySubtracted,LWClimatologySubtracted)
%GlobalRegressionMap(LWClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalRegressionMap(SeaIceClimatologySubtracted,LWclearClimatologySubtracted)
GlobalRegressionMap(SeaIceClimatologySubtracted,SWClimatologySubtracted)
GlobalRegressionMap(SeaIceClimatologySubtracted,SWCFClimatologySubtracted)
GlobalRegressionMap(SeaIceClimatologySubtracted,LWCFClimatologySubtracted)
GlobalRegressionMap(SeaIceClimatologySubtracted,TotalCloudForcingClimatologySubtracted)
GlobalRegressionMap(SeaIceClimatologySubtracted,NetClimatologySubtracted)
GlobalRegressionMap(SeaIceClimatologySubtracted,netclearClimatologySubtracted)
GlobalRegressionMap(SeaIceClimatologySubtracted,SWclearClimatologySubtracted)

GlobalRegressionMap(LWClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalRegressionMap(LWclearClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalRegressionMap(SWClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalRegressionMap(SWclearClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalRegressionMap(SWCFClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalRegressionMap(LWCFClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalRegressionMap(TotalCloudForcingClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalRegressionMap(NetClimatologySubtracted,SeaIceClimatologySubtracted)
GlobalRegressionMap(netclearClimatologySubtracted,SeaIceClimatologySubtracted)

%%%%%%%%%%%%%
PlotSeasonalCorrRegMaps1DTimeSeries(SeaIceClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(SeaIceClimatologySubtracted,SAM)
PlotSeasonalCorrRegMaps1DTimeSeries(SeaIceClimatologySubtracted,NINO34)
PlotSeasonalCorrRegMaps1DTimeSeries(SeaIceClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(SeaIceClimatologySubtracted,T)

%%%%%%%%%%%
%%%%%%%%%%%
Corr1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,T);

Regress1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(EastHeatFluxClimatologySubtracted,T);

%load ClimatologySubtractedVars.mat
%save('ClimatologySubtractedVars.mat','EastHeatFluxClimatologySubtracted','-append')

GlobalCorrMap(LWClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalCorrMap(LWclearClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalCorrMap(SWClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalCorrMap(SWclearClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalCorrMap(SWCFClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalCorrMap(LWCFClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalCorrMap(TotalCloudForcingClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalCorrMap(NetClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalCorrMap(netclearClimatologySubtracted,EastHeatFluxClimatologySubtracted)

GlobalRegressionMap(EastHeatFluxClimatologySubtracted,LWClimatologySubtracted)
%GlobalRegressionMap(LWClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalRegressionMap(EastHeatFluxClimatologySubtracted,LWclearClimatologySubtracted)
GlobalRegressionMap(EastHeatFluxClimatologySubtracted,SWClimatologySubtracted)
GlobalRegressionMap(EastHeatFluxClimatologySubtracted,SWCFClimatologySubtracted)
GlobalRegressionMap(EastHeatFluxClimatologySubtracted,LWCFClimatologySubtracted)
GlobalRegressionMap(EastHeatFluxClimatologySubtracted,TotalCloudForcingClimatologySubtracted)
GlobalRegressionMap(EastHeatFluxClimatologySubtracted,NetClimatologySubtracted)
GlobalRegressionMap(EastHeatFluxClimatologySubtracted,netclearClimatologySubtracted)
GlobalRegressionMap(EastHeatFluxClimatologySubtracted,SWclearClimatologySubtracted)

GlobalRegressionMap(LWClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalRegressionMap(LWclearClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalRegressionMap(SWClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalRegressionMap(SWclearClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalRegressionMap(SWCFClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalRegressionMap(LWCFClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalRegressionMap(TotalCloudForcingClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalRegressionMap(NetClimatologySubtracted,EastHeatFluxClimatologySubtracted)
GlobalRegressionMap(netclearClimatologySubtracted,EastHeatFluxClimatologySubtracted)

%%%%%%%%%%%%%
PlotSeasonalCorrRegMaps1DTimeSeries(EastHeatFluxClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(EastHeatFluxClimatologySubtracted,SAM)
PlotSeasonalCorrRegMaps1DTimeSeries(EastHeatFluxClimatologySubtracted,NINO34)
PlotSeasonalCorrRegMaps1DTimeSeries(EastHeatFluxClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(EastHeatFluxClimatologySubtracted,T)

%%%%
hold off;
Corr1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,NAM);
Corr1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,SAM);
Corr1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,NAO);
Corr1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,PNA);
Corr1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,NINO34);
Corr1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,SAO);
Corr1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,T);

Regress1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,NAM);
Regress1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,SAM);
Regress1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,NAO);
Regress1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,PNA);
Regress1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,NINO34);
Regress1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,SAO);
Regress1DTimeSeriesMap(NorthHeatFluxClimatologySubtracted,T);

%load ClimatologySubtractedVars.mat
%save('ClimatologySubtractedVars.mat','NorthHeatFluxClimatologySubtracted','-append')

GlobalCorrMap(LWClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalCorrMap(LWclearClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalCorrMap(SWClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalCorrMap(SWclearClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalCorrMap(SWCFClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalCorrMap(LWCFClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalCorrMap(TotalCloudForcingClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalCorrMap(NetClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalCorrMap(netclearClimatologySubtracted,NorthHeatFluxClimatologySubtracted)

GlobalRegressionMap(NorthHeatFluxClimatologySubtracted,LWClimatologySubtracted)
%GlobalRegressionMap(LWClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalRegressionMap(NorthHeatFluxClimatologySubtracted,LWclearClimatologySubtracted)
GlobalRegressionMap(NorthHeatFluxClimatologySubtracted,SWClimatologySubtracted)
GlobalRegressionMap(NorthHeatFluxClimatologySubtracted,SWCFClimatologySubtracted)
GlobalRegressionMap(NorthHeatFluxClimatologySubtracted,LWCFClimatologySubtracted)
GlobalRegressionMap(NorthHeatFluxClimatologySubtracted,TotalCloudForcingClimatologySubtracted)
GlobalRegressionMap(NorthHeatFluxClimatologySubtracted,NetClimatologySubtracted)
GlobalRegressionMap(NorthHeatFluxClimatologySubtracted,netclearClimatologySubtracted)
GlobalRegressionMap(NorthHeatFluxClimatologySubtracted,SWclearClimatologySubtracted)

GlobalRegressionMap(LWClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalRegressionMap(LWclearClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalRegressionMap(SWClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalRegressionMap(SWclearClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalRegressionMap(SWCFClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalRegressionMap(LWCFClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalRegressionMap(TotalCloudForcingClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalRegressionMap(NetClimatologySubtracted,NorthHeatFluxClimatologySubtracted)
GlobalRegressionMap(netclearClimatologySubtracted,NorthHeatFluxClimatologySubtracted)

%%%%%%%%%%%%%
PlotSeasonalCorrRegMaps1DTimeSeries(NorthHeatFluxClimatologySubtracted,NAM)
PlotSeasonalCorrRegMaps1DTimeSeries(NorthHeatFluxClimatologySubtracted,SAM)
PlotSeasonalCorrRegMaps1DTimeSeries(NorthHeatFluxClimatologySubtracted,NINO34)
PlotSeasonalCorrRegMaps1DTimeSeries(NorthHeatFluxClimatologySubtracted,NAO)
PlotSeasonalCorrRegMaps1DTimeSeries(NorthHeatFluxClimatologySubtracted,T)


