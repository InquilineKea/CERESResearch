
% load(fullfile('D:\2013-Research\MonthlyFluxDeparturesStruct.Mat'),'Net');
load MonthlyFluxDeparturesStruct.mat
FluxNames = fieldnames(MonthlyFluxDepartures);
load LatitudinalDifFluxTimeSeriesAsField.mat
load Indices.mat
Indices = rmfield(Indices,'PNA');
Indices = rmfield(Indices,'SAO');
Indices = rmfield(Indices,'T');
IndexNames = fieldnames(Indices);
MonthlyFluxDepartures = rmfield(MonthlyFluxDepartures,'Clouds');
MonthlyFluxDepartures = rmfield(MonthlyFluxDepartures,'Precip');
FluxNames = fieldnames(MonthlyFluxDepartures);

NormIndices.NAM = Indices.NAM/std(Indices.NAM);
NormIndices.SAM= Indices.SAM/std(Indices.SAM);
NormIndices.NAO = Indices.NAO/std(Indices.NAO);
NormIndices.NINO34= Indices.NINO34/std(Indices.NINO34);

TimeSeries = Indices.NAM;
Flux1 = MonthlyFluxDepartures.Net;
i=1;
[a,a,a,MaxBeta(i)] = Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.NAM,FluxNames{i},'NAM')

% clear MonthlyFluxDepartures
MaxBeta = zeros(1,4);

% subplot(1,2,1)
% [a,a,a,MaxBeta(1)] = Regress1DTimeSeriesMapWithCorrContours(Flux1,Indices.NAM,'Net','NAM');
% subplot(1,2,2)
% [a,a,a,MaxBeta(2)] = Regress1DTimeSeriesMapWithCorrContours(MonthlyFluxDepartures.(FluxNames{2}),Indices.NAM,FluxNames{2},'NAM');
% caxis([-2 2])

[a,a,a,MaxBeta(i)] = Regress1DTimeSeriesMapWithCorrContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.NAM,FluxNames{i},'NAM')

for i=1:8
   [a,a,a,MaxBeta(i)] = Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.NAM,FluxNames{i},'NAM')
end
for i=1:8
   [a,a,a,MaxBeta(i)] = Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.SAM,FluxNames{i},'SAM')
end
for i=1:8
   [a,a,a,MaxBeta(i)] = Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.NINO34,FluxNames{i},'NINO34');
end


for i=1:4
   subplot(2,2,i)
   [a,a,a,MaxBeta(i)] = Regress1DTimeSeriesMapWithCorrContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.NAM,FluxNames{i},'NAM');
end
 subplots = get(gcf,'Children');
 AllPositions= get(subplots,'Position');

 for i=1:length(subplots) % for each subplot
    caxis(subplots(i),[-max(MaxBeta),max(MaxBeta)]); % set the clim
 end
set(subplots(1),'Position',get(subplots(1),'OuterPosition'))
set(subplots(2),'Position',get(subplots(2),'OuterPosition'))
set(subplots(3),'Position',get(subplots(3),'OuterPosition'))
set(subplots(4),'Position',get(subplots(4),'OuterPosition'))
set(subplots(4),'Title','test')
MakeLizMap
colormap(lizmap)

suptitle('Net/SW/LW/SWCF Regressed on NAM')
subplots = get(gcf,'Children');
AllPositions= get(subplots,'Position');
h=colorbar;
set(subplots(1),'Position',AllPositions{1})

% pos2 = get(subplots(2),'OuterPosition');
% pos4 = get(subplots(4),'OuterPosition');
% pos2(2) = pos4(2)-pow2(4);
% set(subplots(2),'Position',pos2)

 %%%
%  pos1 = get(subplots(8),'Position');
%  OuterPos1 = get(subplots(1),'OuterPosition');
%  pos3 = get(subplots(4),'Position');
% %  pos4= get(subplots(4),'Position');
% pos3(2) = pos1(2) - pos3(4);
% pos4(2) = pos3(2);
% % set(subplots(1),'Position',OuterPos1)
% % set(subplots(2),'Position',OuterPos1)
% % set(subplots(5),'Position',OuterPos1)
% % set(subplots(4),'Position',OuterPos1)
% 
% set(subplots(4),'Position',pos3)
% set(subplots(8),'Position',pos4)
% set(subplots(3),'Position',get(subplots(1),'Position'))
% 
% TestPos1 = pos1;
% TestPos1(2) = 0;
% set(subplots(1),'Position',TestPos1)
%  
% set(subplots(5),'Position',pos1)


 set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['AllFlux','-RegressedOn-','NAM','CorrContour.png']);

load LatWeights.mat
for i=1:length(FluxNames)
%contourf(squeeze(mean(MonthlyFluxDepartures.Net,2)));colorbar
    [X,Y]=meshgrid(1:156,sind(LatWeights(:,1)));
contourf(X,Y,squeeze(mean(bsxfun(@times,MonthlyFluxDepartures.(FluxNames{i}),LatWeights(:,2)),2)));colorbar
grid on;
set(gca,'FontSize',20)
xlabel('Time')
ylabel('Latitude')
set(gca,'xtick',11:12:156)
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
    set(gca,'ytick',sind((0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{i}),2),1))/180:179.5)-90))
    set(gca,'yticklabel',num2cell(-90:10:90))
title(['Mean Deviation of Latitudinal Climatology-Subtracted-',FluxNames{i},' (Watts/m^2)'])
caxis([-max(max(abs(squeeze(mean(bsxfun(@times,MonthlyFluxDepartures.(FluxNames{i}),LatWeights(:,2)),2))))) max(max(abs(squeeze(mean(bsxfun(@times,MonthlyFluxDepartures.(FluxNames{i}),LatWeights(:,2)),2)))))])
    MakeLizMap
    colormap(lizmap)
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['ContourWeightedLatTime',FluxNames{i},'.png']);
end

[NetNHMonthAllLats, NetSHMonthAllLats, NetNHMinusSHMonthAllLats,NHPlusSHFlux] = HemisphericFluxesClimSubtracted(MonthlyFluxDepartures.Net,0,90,'Net');
for i=1:length(FluxNames)
    [blah,blah,blah,NHPlusSHFlux.(FluxNames{i})] = HemisphericFluxesClimSubtracted(MonthlyFluxDepartures.(FluxNames{i}),0,90,FluxNames{i});
end

plot(NHPlusSHFlux)

% NHPlusSHFlux=sum(squeeze(bsxfun(@times,sum(MonthlyFluxDepartures.Net,2),LatWeights(:,2))),1);
% NHPlusSHFlux=sum(squeeze(sum(bsxfun(@times,MonthlyFluxDepartures.Net,LatWeights(:,2)),2)),1);
% AllLatsNetFlux = squeeze(sum(bsxfun(@times,MonthlyFluxDepartures.Net,LatWeights(:,2)),2));
% 
% Difdif=sum(AllLatsNetFlux(1:end/2,:),1)-sum(AllLatsNetFlux(end/2+2:end,:),1);
% plot(Difdif)
% 
% %plot(squeeze(bsxfun(@times,sum(MonthlyFluxDepartures.Net,2),LatWeights(:,2))))
% subplot(2,1,1)
% plot(mean(squeeze(bsxfun(@times,sum(MonthlyFluxDepartures.Net,2),LatWeights(:,2)))))
% subplot(2,1,2)
% plot(NHPlusSHFlux)
% mean(NHPlusSHFlux)
% NHPlusSHFlux=squeeze(bsxfun(@times,MonthlyFluxDepartures.Net,LatWeights(:,2)));
% contourf(squeeze(mean(NHPlusSHFlux,2)))

FluxNames = fieldnames(MonthlyFluxDepartures);

for i=1:length(FluxNames)
Regress1DTimeSeriesMap(MonthlyFluxDepartures.(FluxNames{i}),FluxTimeSeries.Net.AllLats.N,FluxNames{i},'NH Total Net Flux')
Regress1DTimeSeriesMap(MonthlyFluxDepartures.(FluxNames{i}),FluxTimeSeries.Net.AllLats.S,FluxNames{i},'SH Total Net Flux')
Regress1DTimeSeriesMap(MonthlyFluxDepartures.(FluxNames{i}),NHPlusSHFlux,FluxNames{i},'NH+SH Total Net Flux')
end

for i=1:length(FluxNames)
Regress1DTimeSeriesMap(MonthlyFluxDepartures.Net,NHPlusSHFlux.(FluxNames{i}),'Net',['NH+SH of ',FluxNames{i}])
end

load LatitudinalDifFluxTimeSeriesAsField.mat

% FluxLatitude.SW = [FluxTimeSeries.SW.Lat49_90.S;FluxTimeSeries.SW.Lat30_49.S;...
%     FluxTimeSeries.SW.Lat15_30.S;FluxTimeSeries.SW.Lat0_15.S;FluxTimeSeries.SW.Lat0_15.N;
%     FluxTimeSeries.SW.Lat15_30.N;FluxTimeSeries.SW.Lat30_49.N;FluxTimeSeries.SW.Lat49_90.N;...
%     FluxTimeSeries.SW.AllLats.S;FluxTimeSeries.SW.AllLats.N;FluxTimeSeries.SW.Lat49_90.Dif;...
%     FluxTimeSeries.SW.Lat30_49.Dif;FluxTimeSeries.SW.Lat15_30.Dif;FluxTimeSeries.SW.Lat0_15.Dif];

fields = fieldnames(FluxTimeSeries);
% FluxTimeSeries.(fields{1})
% FluxTimeSeries.(fields{1}).Lat49_90.Dif
% 
for i=1:length(fields)
   FluxLatitude.(fields{i}) = [FluxTimeSeries.(fields{i}).Lat49_90.S;FluxTimeSeries.(fields{i}).Lat30_49.S;...
    FluxTimeSeries.(fields{i}).Lat15_30.S;FluxTimeSeries.(fields{i}).Lat0_15.S;FluxTimeSeries.(fields{i}).Lat0_15.N;
    FluxTimeSeries.(fields{i}).Lat15_30.N;FluxTimeSeries.(fields{i}).Lat30_49.N;FluxTimeSeries.(fields{i}).Lat49_90.N;...
    FluxTimeSeries.(fields{i}).AllLats.S;FluxTimeSeries.(fields{i}).AllLats.N;FluxTimeSeries.(fields{i}).Lat49_90.Dif;...
    FluxTimeSeries.(fields{i}).Lat30_49.Dif;FluxTimeSeries.(fields{i}).Lat15_30.Dif;FluxTimeSeries.(fields{i}).Lat0_15.Dif];
end


LatitudeBand = {'90-49S','49-30S','30-15S','15-0S','0-15N','15-30N','30-49N','49-90N','SH','NH','49-90Dif','30-49Dif','15-30Dif','0-15Dif'};
NHPlusSHFlux = rmfield(NHPlusSHFlux,'NetClear');

for i=6:length(FluxNames)
    FluxMatrix = FluxLatitude.(FluxNames{i});
    for j=1:14
    Regress1DTimeSeriesMap(MonthlyFluxDepartures.Net,FluxMatrix(j,:),'Net',[FluxNames{i},LatitudeBand{j}])
    end
end




%%%%%%

ensoTime = ncread('enso.cdf','T');
% NINO12 = ncread('enso.cdf','NINO12');
% NINO12 = NINO12(find(ensoTime > 482 & ensoTime < 638));
% NINO3 = ncread('enso.cdf','NINO3');
% NINO3 = NINO3(find(ensoTime > 482 & ensoTime < 638));
FullIndices.NINO34 = ncread('enso.cdf','NINO34')';
% NINO34 = NINO34(find(ensoTime > 482 & ensoTime < 638));
% NINO4 = ncread('enso.cdf','NINO4');
% NINO4 = NINO4(find(ensoTime > 482 & ensoTime < 638));
% SAOTime = ncread('sao.cdf','T');
% SAO = ncread('sao.cdf','anomaly');
% SAO = SAO(find(SAOTime > 482 & SAOTime < 638));
% 
FullIndices.NAM = dlmread('NAM.ascii');
FullIndices.NAM = FullIndices.NAM(:,3)';
FullIndices.NAO = dlmread('NAO.ascii');
FullIndices.NAO = FullIndices.NAO(:,3)';
FullIndices.SAM = dlmread('SAM.ascii');
FullIndices.SAM = FullIndices.SAM(:,3)';
% PNA = dlmread('PNA.ascii');
% T = dlmread('monthly.land_ocean.90S.90N.df_1901-2000mean.dat');
%NAM(find(NAM(:,1) >= 2000),3)
%find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2)
% NAM = NAM(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2),3);
% NAO = NAO(find(NAO(:,1) == 2000 & NAO(:,2)==3):find(NAO(:,1) == 2013 & NAO(:,2)==2),3);
% SAM = SAM(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2),3);
% PNA = PNA(find(PNA(:,1) == 2000 & PNA(:,2)==3):find(PNA(:,1) == 2013 & PNA(:,2)==2),3);
% T = T(find(T(:,1) == 2000 & T(:,2)==3):find(T(:,1) == 2013 & T(:,2)==2),3);

THadCRUT.T = ncread('HadCRUT.4.2.0.0.median.nc','temperature_anomaly');
THadCRUT.Time = ncread('HadCRUT.4.2.0.0.median.nc','time');
THadCRUT.T = permute(THadCRUT.T,[2 1 3]);
THadCRUT.T = circshift(THadCRUT.T,[0 36 0]);
%timeunits         = 'days since 1850-1-1 00:00:00'
size(circshift(THadCRUT.T,[0 36 0]))

% Regress1DTimeSeriesMapWithCorrContours(THadCRUT.T(:,:,end-length(FullIndices.NAM)+1:end),FullIndices.NAM,'HadCRUT-Temp','NAM')
Regress1DTimeSeriesMap(THadCRUT.T(:,:,end-length(FullIndices.NAM)+1:end),FullIndices.NAM,'HadCRUT-Temp','NAM')
Regress1DTimeSeriesMap(THadCRUT.T(:,:,end-length(FullIndices.SAM)+1:end),FullIndices.SAM,'HadCRUT-Temp','SAM');
Regress1DTimeSeriesMap(THadCRUT.T(:,:,end-length(FullIndices.NAO)+1:end),FullIndices.NAO,'HadCRUT-Temp','NAO');
Regress1DTimeSeriesMap(THadCRUT.T(:,:,end-length(FullIndices.NINO34)+1:end),FullIndices.NINO34,'HadCRUT-Temp','NINO34');
Regress1DTimeSeriesMap(THadCRUT.T(:,:,end-length(FullIndices.NAM)+1:end),FullIndices.NINO34(end-length(FullIndices.NAM)+1:end),'HadCRUT-Temp','NINO34');

%regular longitude:
ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lon')%from 0degE to360degE
THadCRUT.Long = ncread('HadCRUT.4.2.0.0.median.nc','longitude'); %from -180degE to180degE. shift to where it starts increasing

NAMFluxes = load('RegressIndexFluxAllTimeNew.mat','NAM');

contourf(NAMFluxes.NAM.Net(:,:,1))
load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(NAMFluxes.NAM.Net(:,:,1),geoidrefvec,'DisplayType','texturemap');colorbar

load('Indices')

FluxNames = fieldnames(NAMFluxes.NAM);
for i=1:length(FluxNames)
[NAMHemFlux.(FluxNames{i}).NH NAMHemFlux.(FluxNames{i}).SH NAMHemFlux.(FluxNames{i}).HemDif NAMHemFlux.(FluxNames{i}).HemSum] = HemisphericFluxesClimSubtracted(NAMFluxes.NAM.(FluxNames{i}),0,90,['NAMon', FluxNames{i}]);
end
IndexHemFluxes.NAM = NAMHemFlux;
clear NAMFluxes
SAMFluxes = load('RegressIndexFluxAllTimeNew.mat','SAM');
for i=1:length(FluxNames)
[SAMHemFlux.(FluxNames{i}).NH SAMHemFlux.(FluxNames{i}).SH SAMHemFlux.(FluxNames{i}).HemDif SAMHemFlux.(FluxNames{i}).HemSum] = HemisphericFluxesClimSubtracted(SAMFluxes.SAM.(FluxNames{i}),0,90,['SAMon', FluxNames{i}]);
end
IndexHemFluxes.SAM= SAMHemFlux;
clear SAMFluxes

NINOFluxes = load('RegressIndexFluxAllTimeNew.mat','NINO34');
FluxNames = fieldnames(NINOFluxes.NINO34);

for i=1:length(FluxNames)
[NINOHemFlux.(FluxNames{i}).NH NINOHemFlux.(FluxNames{i}).SH NINOHemFlux.(FluxNames{i}).HemDif NINOHemFlux.(FluxNames{i}).HemSum] = HemisphericFluxesClimSubtracted(NINOFluxes.NINO34.(FluxNames{i}),0,90,['NINO34on', FluxNames{i}]);
end
clear NINOFluxes

IndexHemFluxes.NINO34 = NINOHemFlux;

for i=1:length(FluxNames)
IndexHemFluxes.CombinedIndex.(FluxNames{i}).NH = NAMHemFlux.(FluxNames{i}).NH+SAMHemFlux.(FluxNames{i}).NH+NINOHemFlux.(FluxNames{i}).NH;
IndexHemFluxes.CombinedIndex.(FluxNames{i}).SH = NAMHemFlux.(FluxNames{i}).SH+SAMHemFlux.(FluxNames{i}).SH+NINOHemFlux.(FluxNames{i}).SH;
IndexHemFluxes.CombinedIndex.(FluxNames{i}).HemDif = NAMHemFlux.(FluxNames{i}).HemDif+SAMHemFlux.(FluxNames{i}).HemDif+NINOHemFlux.(FluxNames{i}).HemDif;
IndexHemFluxes.CombinedIndex.(FluxNames{i}).HemSum = NAMHemFlux.(FluxNames{i}).HemSum+SAMHemFlux.(FluxNames{i}).HemSum+NINOHemFlux.(FluxNames{i}).HemSum;
end

save('IndexHemFluxes.mat','IndexHemFluxes')

mean(IndexHemFluxes.CombinedIndex.Net.HemDif)
mean(IndexHemFluxes.CombinedIndex.SW.HemDif)
mean(IndexHemFluxes.CombinedIndex.LW.HemDif)

NH = IndexHemFluxes.CombinedIndex.Net.NH;
SH = IndexHemFluxes.CombinedIndex.Net.SH;
HemDif= IndexHemFluxes.CombinedIndex.Net.HemDif;
VariableName = 'CombinedIndexOnNet';
LowerLat = 0;
HigherLat = 90;
time=156;

[AX,H1,H2] = plotyy(1:length(NH),[NH' ...
    SH'],1:length(NH),HemDif);
grid on;
set(gca,'xtick',11:12:time)
set(AX(2),'XTickLabel',[])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year');
set(get(AX(1),'Ylabel'),'FontSize',20,'String','Heat Flux in Latitude Band for NH, SH (Watts/m^2)') 
set(get(AX(2),'Ylabel'),'FontSize',20,'String','NH-SH Difference in Total Heat Flux (Watts/m^2)','FontSize',20) 
difSD = sqrt(var(NH) + var(SH) - 2*getfield(cov(NH,SH), {1,2}));
corrcof = getfield(corrcoef(NH,SH),{1,2});
set(legend(['NH SD = ', num2str(std(NH))],['SH SD = ', num2str(std(SH))] ,sprintf(['NH-SH SD = ', num2str(difSD), '\n Corr Cof = ', num2str(corrcof)])),'Location','BestOutside')
set(H1,'linewidth',3)
set(H2,'linewidth',3)% to change the first line
set(AX,'FontSize',20)
set(gcf, 'Units','inches', 'Position',[0 0 20 10])
set(gca, 'Units','inches', 'Position',[1 1 16 8])
title([num2str(LowerLat),'-',num2str(HigherLat),'deg Anomaly ', VariableName])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[VariableName,'_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.png']);
hold off;

%%%%
for i=1:length(FluxNames)
NAMFluxes.NAM.(FluxNames{i}) = NAMFluxes.NAM.(FluxNames{i}) + SAMFluxes.SAM.(FluxNames{i});
end

NINOFluxes = load('RegressIndexFluxAllTimeNew.mat','NINO34');
for i=1:length(FluxNames)
NAMFluxes.NAM.(FluxNames{i}) = NAMFluxes.NAM.(FluxNames{i}) + NINOFluxes.NINO34.(FluxNames{i});
end
clear NINOFluxes

load LatWeights
for j=1:length(FluxNames)
PlotAsContour = bsxfun(@times,squeeze(mean(CombinedIndexFluxes.CombinedIndexFluxes.(FluxNames{j}),2)),LatWeights(:,2));
contourf(PlotAsContour);colorbar
caxis([-max(max(abs(PlotAsContour))) max(max(abs(PlotAsContour)))])
grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
set(gca,'xtick',11:12:size(PlotAsContour,2))
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year');
set(gca,'ytick',0.5:10*(size(PlotAsContour,1)/180):180.5)
set(gca,'yticklabel',num2cell(-90:10:90))
ylabel('latitude')
title(['WeightedTime-Latitudinal Contour Plot of ', FluxNames{j},' over CombinedIndexFluxes'])
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['WeightedTimeLatRegressContours',FluxNames{j},'_','CombinedIndexFluxes','.png']);
hold off;
end

load('ClimatologySubtractedVars.mat','NetClimatologySubtracted');
load('ClimatologySubtractedVars.mat','LWClimatologySubtracted');
load('ClimatologySubtractedVars.mat','SWClimatologySubtracted');

Net.NetClimatologySubtracted;
load('CombinedIndexFluxAnomaly.mat');

GlobalCorrMapVarNameInput(Net.NetClimatologySubtracted,CombinedIndexFluxes.CombinedIndexFluxes.Net,'Net','CombinedIndexRegress')
GlobalRegressionMap(Net.NetClimatologySubtracted,CombinedIndexFluxes.CombinedIndexFluxes.Net,'Net','CombinedIndexRegress');
GlobalCorrMapVarNameInput(LWClimatologySubtracted,CombinedIndexFluxes.CombinedIndexFluxes.LW,'LW','CombinedIndexRegress');
GlobalCorrMapVarNameInput(SWClimatologySubtracted,CombinedIndexFluxes.CombinedIndexFluxes.SW,'SW','CombinedIndexRegress');



IndicesRemoved.Net = Net.NetClimatologySubtracted-CombinedIndexFluxes.CombinedIndexFluxes.Net;
IndicesRemoved.LW = LWClimatologySubtracted-CombinedIndexFluxes.CombinedIndexFluxes.LW;
IndicesRemoved.SW= SWClimatologySubtracted-CombinedIndexFluxes.CombinedIndexFluxes.SW;

GlobalValuesMinusClimatology(IndicesRemoved.Net,'Net after Index Subtraction')
GlobalValuesMinusClimatology(IndicesRemoved.LW,'LW after Index Subtraction');
GlobalValuesMinusClimatology(IndicesRemoved.SW,'SW after Index Subtraction');


GlobalValuesMinusClimatology(CombinedIndexFluxes.CombinedIndexFluxes.Net,'Explained by the Indices')
GlobalValuesMinusClimatology(CombinedIndexFluxes.CombinedIndexFluxes.LW,'LWFromIndices')
GlobalValuesMinusClimatology(CombinedIndexFluxes.CombinedIndexFluxes.SW,'SWFromIndices')

HemisphericFluxesClimSubtracted(IndicesRemoved.Net,0,90,'Net after Indices Removed')
GlobalCorrMapVarNameInput(Net.NetClimatologySubtracted,IndicesRemoved.Net,'Net','Net Indices Removed')
%they're definitely not in the same units...one has to do more to remove
%theinfluence of the regression coefficients...

load('Indices.mat')
GlobalCorrMapVarNameInput()
NAM = dlmread('NAM.ascii');
NAO = dlmread('NAO.ascii');
SAM = dlmread('SAM.ascii');

TimeLaggedIndices3.NAM = NAM(end-156-3+1:end-3);

Corr1DTimeSeriesMap(LWClimatologySubtracted,TimeLaggedIndices3.NAM,'LW','NAM 3 months ago')
%not really independent.. but you could apply an ARMA model to make
%adjacent time series values more independent of eachother? b/c you obv.
%can't just add regressions from time t to regressions from t-1 because
%there is some correlation..

for n=-3:12
%Corr1DTimeSeriesMap(LWClimatologySubtracted,NAM(end-156-n+1:end-n),'LW',['NAM ', num2str(n) ,' months ago'])
Regress1DTimeSeriesMapWithCorrContours(LWClimatologySubtracted,NAM(min(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-n:max(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-n,3),'LW',['NAM ', num2str(n) ,' months ago'])
Regress1DTimeSeriesMapWithCorrContours(SWClimatologySubtracted,NAM(min(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-n:max(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-n,3),'SW',['NAM ', num2str(n) ,' months ago'])
Regress1DTimeSeriesMapWithCorrContours(NetClimatologySubtracted,NAM(min(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-n:max(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-n,3),'Net',['NAM ', num2str(n) ,' months ago'])
end

ensoTime = ncread('enso.cdf','T');
NINO34 = ncread('enso.cdf','NINO34');
for n=-5:5
%Corr1DTimeSeriesMap(LWClimatologySubtracted,NAM(end-156-n+1:end-n),'LW',['NAM ', num2str(n) ,' months ago'])
Regress1DTimeSeriesMapWithCorrContours(LWClimatologySubtracted,NINO34(min(find(ensoTime > 482 & ensoTime < 638))-n:max(find(ensoTime > 482 & ensoTime < 638))-n),'LW',['NINO34At', num2str(n) ,' months ago'])
Regress1DTimeSeriesMapWithCorrContours(SWClimatologySubtracted,NINO34(min(find(ensoTime > 482 & ensoTime < 638))-n:max(find(ensoTime > 482 & ensoTime < 638))-n),'SW',['NINO34At', num2str(n) ,' months ago'])
Regress1DTimeSeriesMapWithCorrContours(NetClimatologySubtracted,NINO34(min(find(ensoTime > 482 & ensoTime < 638))-n:max(find(ensoTime > 482 & ensoTime < 638))-n),'Net',['NINO34At', num2str(n) ,' months ago'])
end
SAM = dlmread('SAM.ascii');
%SAM = SAM(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2),3);
for n=-5:5
%Corr1DTimeSeriesMap(LWClimatologySubtracted,SAM(end-156-n+1:end-n),'LW',['SAM ', num2str(n) ,' months ago'])
Regress1DTimeSeriesMapWithCorrContours(LWClimatologySubtracted,SAM(min(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-n:max(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-n,3),'LW',['SAM ', num2str(n) ,' months ago'])
Regress1DTimeSeriesMapWithCorrContours(SWClimatologySubtracted,SAM(min(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-n:max(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-n,3),'SW',['SAM ', num2str(n) ,' months ago'])
Regress1DTimeSeriesMapWithCorrContours(NetClimatologySubtracted,SAM(min(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-n:max(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-n,3),'Net',['SAM ', num2str(n) ,' months ago'])
end

MonthsAgo = 40;

for n=-7:MonthsAgo
   ENSOAutoCorr(n+8) = corr(NINO34(min(find(ensoTime > 482 & ensoTime < 638)):max(find(ensoTime > 482 & ensoTime < 638))),NINO34(min(find(ensoTime > 482 & ensoTime < 638))-n:max(find(ensoTime > 482 & ensoTime < 638))-n));
   SAMAutoCorr(n+8)=corr(SAM(min(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-0:max(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-0,3),SAM(min(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-n:max(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2013 & SAM(:,2)==2))-n,3));
   NAMAutoCorr(n+8)=corr(NAM(min(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-0:max(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-0,3),NAM(min(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-n:max(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2))-n,3));
end
plot([ENSOAutoCorr' SAMAutoCorr' NAMAutoCorr'],'LineWidth',3)
title('Autocorrelation at X months ago')
set(gca,'xtick',1:1:7+MonthsAgo)
    set(gca,'xticklabel',num2cell(-7:1:7+MonthsAgo))
legend('ENSO','SAM','NAM')
grid on;