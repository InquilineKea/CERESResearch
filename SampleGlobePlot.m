load('MovingAvg.mat', 'MovingAvg');

f = figure('Visible','off')
% f = figure('Visible','on')
load geoid;ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land,'FaceColor', [0 0 0])


    curPlot = MovingAvg.Window12Months.Net(:,:,1);
    geoshow(curPlot,geoidrefvec,'DisplayType','texturemap');
    
    colorbar
load coast
plotm(lat,long,'black')


set(gcf,'Renderer','painters')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['GEOSHOW_TEST.png']);
hold off;
