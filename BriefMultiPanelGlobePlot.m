MakeLizMap
colormap(lizmap)

VariableName = MovingAvg.Window12Months.(FluxNames{1});

load geoid;
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(mean(VariableName,3),geoidrefvec,'DisplayType','texturemap');h=colorbar
set(gca,'FontSize',20);set(h,'FontSize',20)
view(2)
% title(['Temporal ',num2str(MonthFilterSize),' Month Running Mean of ',VariableName])
title([VariableName, ' Mean (W/m^2)'])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])

% subplot(1, 4, 1 : 3)
% % generate first plot
% subplot(1, 4, 4)
% % generate second plot

set(gca,'FontSize',20)
plot(sind(LatWeights(:,1)),std(squeeze(mean(VariableName,2)),0,2),'LineWidth',3)
%std of latitudinal averaged flux over time. this is prolly preferable
grid on;
% set(gca,'xtick',(0.5-90:10:179.5-90))
set(gca,'xtick',sind((0.5:10*180/180:179.5)-90))
set(gca,'xticklabel',num2cell(-90:10:80))
% title(['Temporal SD of ',num2str(MonthFilterSize),' Month MA in ',VariableName])
title([VariableName, ' SD'])
xlabel('Latitude')
ylabel('SD (W/m^2)')
view(90,-90)
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
