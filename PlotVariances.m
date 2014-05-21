
load('2014-04-30.mat', 'RawFlux');

% load('2014-04-30.mat', 'MovingAvg.Window12Months');
% RawFlux = MovingAvg.Window12Months;

TempFluxNames = fieldnames(RawFlux);


for i=1:length(TempFluxNames)
f = figure('Visible','off')

MakeLizMap;colormap(lizmap);set(gca,'FontSize',20)
load geoid; ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);
geoshow(var(RawFlux.(TempFluxNames{i}),0,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
load coast;plotm(lat,long,'black')


caxis([0 max(max(var(RawFlux.(TempFluxNames{i}))))])
title(['Variance of ', TempFluxNames{i}])
grid on;
 set(gcf,'Renderer','painters')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
saveas(gcf,[TempFluxNames{i},'12MnthMAVariance.fig'],'fig')
print(gcf,'-dpng','-r300',[TempFluxNames{i},'12MnthMAVariance.png']);

hold off;
end
