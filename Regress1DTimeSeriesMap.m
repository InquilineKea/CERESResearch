function [Beta,NHMinusSHIndividualLatitudes,FluxSumsAllLatitudeLines] = Regress1DTimeSeriesMap(Flux1,TimeSeries,Name1,Name2)

%time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));
lats = size(Flux1,1);
longs = size(Flux1,2);
time=size(Flux1,3);
Factor = 180/lats;

%TimeSeries = NetDifMonthTotal;
%Flux1 = NetClimatologySubtracted;
TimeSeries3D = ones(1,1,time);
TimeSeries3D(1,1,:) = TimeSeries;
TimeSeries3D = repmat(TimeSeries3D,[180/Factor 360/Factor 1]);

%Beta = bsxfun(@rdivide,time*(sum(Flux1.*TimeSeries3D,3))-sum(Flux1,3).*sum(TimeSeries3D,3),time*sum(TimeSeries3D.*TimeSeries3D,3)-sum(TimeSeries3D,3).*sum(TimeSeries3D,3))*std(TimeSeries);
%multiply by std(TimeSeries) in end so that you come out with units of
Beta = bsxfun(@rdivide,time*(sum(Flux1.*TimeSeries3D,3))-sum(Flux1,3).*sum(TimeSeries3D,3),time*sum(TimeSeries3D.*TimeSeries3D,3)-sum(TimeSeries3D,3).*sum(TimeSeries3D,3));

%Watts/m^2 rather than Watts/m^2/Index

%GlobalCorrs = bsxfun(@rdivide,sum(Flux1.*TimeSeries3D,3),(sqrt(sum(Flux1.^2,3)).*sqrt(sum(TimeSeries3D.^2,3))));
load LatWeights.mat
RegressMean = mean(mean(resizem(Beta,Factor),2).*LatWeights(:,2)); 

FluxSumsAllLatitudeLines = sum(resizem(Beta,Factor),2).*LatWeights(:,2);
NHPlusSHTotalFluxContribution = sum(FluxSumsAllLatitudeLines);
NHTotalFluxContribution = sum(FluxSumsAllLatitudeLines(end/2+1:end));
SHTotalFluxContribution = sum(FluxSumsAllLatitudeLines(1:end/2));
NHMinusSHTotalFluxContribution = NHTotalFluxContribution - SHTotalFluxContribution;
NHMinusSHIndividualLatitudes = FluxSumsAllLatitudeLines(end/2+1:end)-flipud(FluxSumsAllLatitudeLines(1:end/2));
%FLIPUD ONLY WORKS ON ENTRIES WITHIN THE SAME COLUMN. SO MAY NEED TO TRANSPOSE THIS.

HemisphericFluxes = [NHTotalFluxContribution, SHTotalFluxContribution, NHMinusSHTotalFluxContribution,NHPlusSHTotalFluxContribution];

load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(resizem(Beta,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
MakeLizMap;
colormap(lizmap)
caxis([-max(max(Beta)) max(max(Beta))])
title(['Regression of ',Name1,' over Time Series of ',Name2, '. Mean = ', num2str(RegressMean)])
grid on;
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[Name1,'-RegressedOn-',Name2,'.png']);
hold off;

end