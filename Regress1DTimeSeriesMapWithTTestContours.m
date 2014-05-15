function [Beta,NHMinusSHIndividualLatitudes,FluxSumsAllLatitudeLines,MaxBeta,Alpha] = Regress1DTimeSeriesMapWithTTestContours(Flux1,TimeSeries,Name1,Name2)

%time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));
lats = size(Flux1,1);
longs = size(Flux1,2);
time=size(Flux1,3);
time=length(TimeSeries);
Factor = 180/lats;

%TimeSeries = NetDifMonthTotal;
%Flux1 = NetClimatologySubtracted;
TimeSeries3D = ones(1,1,time);
TimeSeries3D(1,1,:) = TimeSeries;
TimeSeries3D = repmat(TimeSeries3D,[180/Factor 360/Factor 1]);

size(Flux1)
size(TimeSeries3D)

%Beta = bsxfun(@rdivide,time*(sum(Flux1.*TimeSeries3D,3))-sum(Flux1,3).*sum(TimeSeries3D,3),time*sum(TimeSeries3D.*TimeSeries3D,3)-sum(TimeSeries3D,3).*sum(TimeSeries3D,3))*std(TimeSeries);
%multiply by std(TimeSeries) in end so that you come out with units of
Sxy = sum(Flux1.*TimeSeries3D,3);
Sxx = sum(TimeSeries3D.*TimeSeries3D,3);
Syy = sum(Flux1.*Flux1,3);
SSE = Syy-Sxy.*Sxy./Sxx;
Sy = sum(Flux1,3);
Sx = sum(TimeSeries);
Beta = bsxfun(@rdivide,time*(Sxy)-Sx*Sy,time*Sxx-Sx.*Sx);
[x,y] = autocorr(TimeSeries);
AutoCorrTime = length(x(x > 1/exp(1))) + 0.5;
DOF = time/(2*AutoCorrTime);
tScore = Beta*sqrt(DOF-2)./sqrt(SSE/sum((TimeSeries-mean(TimeSeries)).^2));
tTestAlpha = 0.05;
pValue = 1-tcdf(tScore,DOF-1);

Alpha = 1/time*(Sy-Beta*Sx);

%Watts/m^2 rather than Watts/m^2/Index

%tScore = bsxfun(@rdivide,sum(Flux1.*TimeSeries3D,3),(sqrt(sum(Flux1.^2,3)).*sqrt(sum(TimeSeries3D.^2,3))));
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

MakeLizMap;
colormap(lizmap)
load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(resizem(Beta,Factor),geoidrefvec,'DisplayType','texturemap');colorbar

hold on;
[c1,h1] = contourm(pValue,geoidrefvec, (1-tTestAlpha/2):1, '--g','LineWidth',2);
[c2,h2] = contourm(pValue,geoidrefvec, tTestAlpha/2:tTestAlpha, 'g','LineWidth',2);
MaxBeta = max(max(abs(Beta)));
caxis([-MaxBeta MaxBeta])
% caxis([-15 15])
title(['Regression of ',Name1,' over Time Series of ',Name2, ' w/T-Test Contours at 0.05 Sig. Mean = ', num2str(RegressMean)])
grid on;
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[Name1,'_',Name2,'RegressTTest.png']);
hold off;

end

%geoshow(geoid,geoidrefvec,'DisplayType','surface');colorbar
%geoshow(tScore,geoidrefvec,'DisplayType','surface');colorbar
% subplot(2,1,1)
% geoshow(resizem(Alpha,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% subplot(2,1,2)
% geoshow(resizem(Sy,Factor),geoidrefvec,'DisplayType','texturemap');colorbar

% geoshow(resizem(sqrt(SSE/sum((TimeSeries-mean(TimeSeries)).^2)),Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% geoshow(resizem(Beta./sqrt(SSE/sum((TimeSeries-mean(TimeSeries)).^2)),Factor),geoidrefvec,'DisplayType','texturemap');colorbar

%  geoshow(resizem(tScore,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% geoshow(resizem(SSE,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% geoshow(resizem(Sxy,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% [c1,h1] = contourm(tScore,geoidrefvec, 3:7, '--g','LineWidth',2);
% [c2,h2] = contourm(tScore,geoidrefvec, -3:-7, 'g','LineWidth',2);
%[c2,h2] = contourm(tScore,geoidrefvec, 1:-.1:.2,'k','LineWidth',1);
% ht1 = clabelm(c1,h1,'LabelSpacing',144);
% ht2 = clabelm(c2,h2,'LabelSpacing',144);

% set(ht1,'Color','r','BackgroundColor','white','FontWeight','bold')
% set(ht2,'Color','r','BackgroundColor','white','FontWeight','bold')
% uistack(ht1,'top')
% 
% uistack(ht2,'top')
% geoshow(contourm(geoid,geoidrefvec, 'LineStyle', 'none'),geoidrefvec,'DisplayType','texturemap');colorbar
