function [FluxDeparture] = MovingAverageFilterAllTiles(Flux,VariableName, MonthFilterSize)

%VariableName = 'Net';
%Flux=Net; %comment out

permuteFlux = permute(Flux,[3 1 2]);

lat = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');
load LatWeights.mat

time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));

MovingAverage = zeros(180,360,time);

days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
allMonthWeights = repmat(days_per_month,1,12);
allMonthWeights = [31 30 31 30 31 31 30 31 30 31 allMonthWeights 31 28];
allMonthWeights(48) = 29; %leap year, http://www.wolframalpha.com/input/?i=months+between+march+2000+and+february+2004
allMonthWeights(96) = 29;

allMonthWeights3D = ones(1,1,time);
allMonthWeights3D(1,1,:) = allMonthWeights;
allMonthWeights3D = repmat(allMonthWeights3D,[180 360 1]);
allMonthWeights3D = permute(allMonthWeights3D,[3 1 2]);

moving_sum = @(n, x) filter(ones(1,n), 1, x);
MovingAverage= bsxfun(@rdivide,moving_sum(MonthFilterSize, permuteFlux .* allMonthWeights3D) ...
    ,moving_sum(MonthFilterSize,allMonthWeights3D)); %or Mean Lat Flux here
MovingAverage = permute(MovingAverage,[2 3 1]);

%MWA means 1st 11 values are junk

MovingAverage = MovingAverage(:,:,MonthFilterSize:end);

NumMonthsAfterMAFiltering = time-(MonthFilterSize-1);

FluxDeparture = zeros(180,360,NumMonthsAfterMAFiltering);
%Flux = LWclear;

for i=1:NumMonthsAfterMAFiltering
    Indices = mod(i,12):12:NumMonthsAfterMAFiltering;
    if mod(i,12)==0
        Indices = 12:12:NumMonthsAfterMAFiltering;
    end
    FluxDeparture(:,:,i) = MovingAverage(:,:,i) - mean(MovingAverage(:,:,Indices),3);
end

MakeLizMap
colormap(lizmap)

load geoid;
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(mean(MovingAverage,3),geoidrefvec,'DisplayType','texturemap');colorbar
title(['Temporal ',num2str(MonthFilterSize),' Month Running Mean of ',VariableName])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['GlobalMean_', num2str(MonthFilterSize) ,'Month_MA_', VariableName, '.png']);

hold off;

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(MovingAverage,0,3),geoidrefvec,'DisplayType','texturemap');colorbar %flux departure; doesnt matter here
title(['Temporal ',num2str(MonthFilterSize),' Month Moving Average SD of ',VariableName])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['GlobalSD', num2str(MonthFilterSize) ,'Month_MA_', VariableName, '.png']);

hold off;

%%CONVERT TO WATTS/M^2!!

[X,Y]=meshgrid(1:size(MovingAverage,3),sind(LatWeights(:,1)));
contourf(X,Y,squeeze(mean(FluxDeparture,2)),30);colorbar
shading flat;
caxis([-max(max(abs(squeeze(mean(FluxDeparture,2))))) max(max(abs(squeeze(mean(FluxDeparture,2)))))])
grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
set(gca,'FontSize',20)
xlabel('Year End')
ylabel('Latitude')
set(gca,'xtick',[12-2-(MonthFilterSize-1):12:time]) %what if Filter Size is 12 months?
[12-2-(MonthFilterSize-1):12:NumMonthsAfterMAFiltering]
set(gca,'XTickLabel',2000:2012)
set(gca,'ytick',sind((0.5:10*180/180:179.5)-90))
set(gca,'yticklabel',num2cell(-90:10:90))
% set(gca,'ytick',10:10:180)
% set(gca,'YTickLabel',-80:10:90)
title(['Mean Deviation of Latitudinal ',num2str(MonthFilterSize),' Month MA of ',VariableName,' (W/m^2)'])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[VariableName, '_', num2str(MonthFilterSize) ,'Month_MA_TimeLat_Contour_','.png']);

hold off;

set(gca,'FontSize',20)
plot(sind(LatWeights(:,1)),std(squeeze(mean(FluxDeparture,2)),0,2),'LineWidth',3)
%std of latitudinal averaged flux over time. this is prolly preferable
grid on;
% set(gca,'xtick',(0.5-90:10:179.5-90))
set(gca,'xtick',sind((0.5:10*180/180:179.5)-90))
set(gca,'xticklabel',num2cell(-90:10:80))
title(['Temporal SD of ',num2str(MonthFilterSize),' Month MA in ',VariableName])
xlabel('Latitude')
ylabel('Temporal SD')
view(90,-90)
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[VariableName, num2str(MonthFilterSize),'Month_MA_Latitudinal_SD','.png']);

hold off;

set(gca,'FontSize',20)
plot(mean(squeeze(std(FluxDeparture,0,3)),2))
% latitudinal mean of ALL STDs along constant latitude circle

LatMWA = squeeze(mean(FluxDeparture,2));
bb = corr(flipud(LatMWA(1:end/2,:))',LatMWA(end/2+1:end,:)')
set(gca,'FontSize',20)
plot(diag(bb),'LineWidth',3)
grid on;
xlabel('Latitude (degrees away from equator)')
ylabel('Correlation Coefficient')
title(['Correlation between between X Latitude N and X Latitude S in MWA of ', VariableName])
view(90,-90)
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[VariableName, num2str(MonthFilterSize) ,'Month_MA_Latitudinal_North_vs_South_Correlation','.png']);

end