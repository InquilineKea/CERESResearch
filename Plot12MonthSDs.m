
load('2014-04-30.mat', 'MovingAvg.Window12Months');

% load('2014-04-30.mat', 'MovingAvg.Window12Months');
% MovingAvg.Window12Months = MovingAvg.Window12Months;

TempFluxNames = fieldnames(MovingAvg.Window12Months);


for i=1:length(TempFluxNames)
f = figure('Visible','off')

MakeLizMap;colormap(lizmap);set(gca,'FontSize',20)
load geoid; ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);
geoshow(std(MovingAvg.Window12Months.(TempFluxNames{i}),0,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
load coast;plotm(lat,long,'black')


caxis([0 max(max(std(MovingAvg.Window12Months.(TempFluxNames{i}))))])
title(['SD of ', TempFluxNames{i}])
grid on;
 set(gcf,'Renderer','painters')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
saveas(gcf,[TempFluxNames{i},'12MnthMA_SD.fig'],'fig')
print(gcf,'-dpng','-r300',[TempFluxNames{i},'12MnthMA_SD.png']);

hold off;
end
