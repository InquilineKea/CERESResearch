function [Beta,NHMinusSHIndividualLatitudes,FluxSumsAllLatitudeLines,MaxBeta] = Regress1DTimeSeriesMapWithCorrContours(Flux1,TimeSeries,Name1,Name2)

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
GlobalCorrs = bsxfun(@rdivide,sum(Flux1.*TimeSeries3D,3),(sqrt(sum(Flux1.^2,3)).*sqrt(sum(TimeSeries3D.^2,3))));

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
%geoshow(geoid,geoidrefvec,'DisplayType','surface');colorbar
%geoshow(GlobalCorrs,geoidrefvec,'DisplayType','surface');colorbar
geoshow(resizem(Beta,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
hold on;
[c1,h1] = contourm(GlobalCorrs,geoidrefvec, -1:.1:-.2, '--g','LineWidth',2);
[c2,h2] = contourm(GlobalCorrs,geoidrefvec, .2:.1:1,'g','LineWidth',2);
%[c2,h2] = contourm(GlobalCorrs,geoidrefvec, 1:-.1:.2,'k','LineWidth',1);
% ht1 = clabelm(c1,h1,'LabelSpacing',144);
% ht2 = clabelm(c2,h2,'LabelSpacing',144);

% set(ht1,'Color','r','BackgroundColor','white','FontWeight','bold')
% set(ht2,'Color','r','BackgroundColor','white','FontWeight','bold')
% uistack(ht1,'top')
% 
% uistack(ht2,'top')
% geoshow(contourm(geoid,geoidrefvec, 'LineStyle', 'none'),geoidrefvec,'DisplayType','texturemap');colorbar
MaxBeta = max(max(abs(Beta)));
MakeLizMap;
colormap(lizmap)
caxis([-MaxBeta MaxBeta])
caxis([-15 15])
title(['Regression of ',Name1,' over Time Series of ',Name2, ' w/Corr Contours. Mean = ', num2str(RegressMean)])
grid on;
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[Name1,'-RegressedOn-',Name2,'CorrContour.png']);
hold off;

end