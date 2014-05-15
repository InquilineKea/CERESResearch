function [FluxDeparture] = AllTileGlobalValuesMinusClimatology(Flux,VariableName)
%climatology subtracted global values


time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));
load geoid;

FluxDeparture = zeros(180,360,time);
%Flux = LWclear;

for i=1:time
    Indices = mod(i,12):12:time;
    if mod(i,12)==0
        Indices = 12:12:time;
    end
    FluxDeparture(:,:,i) = Flux(:,:,i) - mean(Flux(:,:,Indices),3);
end

%geoshow(mean(Flux(:,:,11:12:end),3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(FluxDeparture,0,3),geoidrefvec,'DisplayType','texturemap');colorbar
title(['Temporal SD of Climatology-Subtracted-',VariableName])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['GlobalSDClimatologySubtracted', VariableName, '.png']);

hold off;

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(mean(Flux,3),geoidrefvec,'DisplayType','texturemap');colorbar
title(['Mean of ',VariableName])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Mean', VariableName, '.png']);

hold off;
% geoshow(FluxDeparture(:,:,1),geoidrefvec,'DisplayType','texturemap');colorbar%March 2000 deviations

%contourf(squeeze(mean(Flux,2)));colorbar

contourf(squeeze(mean(FluxDeparture,2)));colorbar
grid on;
set(gca,'FontSize',20)
xlabel('Time')
ylabel('Latitude')
set(gca,'xtick',11:12:time)
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'ytick',10:10:180)
set(gca,'YTickLabel',-80:10:90)
title(['Mean Deviation of Latitudinal Climatology-Subtracted-',VariableName,' (Watts)'])
caxis([-max(abs(squeeze(mean(FluxDeparture,2)))) max(abs(squeeze(mean(FluxDeparture,2))))])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['ContourLatitudinalClimatologySubtractedTimeVariation',VariableName,'.png']);

hold off;

%they might look small maybe because anomalously warm regions in one area
%may be different from anomalously warm regions in another. but also,
%radiation is going to vary a lot less than temperature.. check LWClear for
%a clear signal of global warming
plot(std(squeeze(mean(FluxDeparture,2)),0,2))%prolly wrong
set(gca,'FontSize',20)
%plot(mean(squeeze(std(FluxDeparture,0,3)),2))
grid on;
set(gca,'xtick',10:10:180)
set(gca,'XTickLabel',-80:10:90)
title(['Temporal SD of Climatology-Subtracted-',VariableName])
xlabel('Latitude')
ylabel('Temporal SD')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['LatitudinalClimatologySubtractedSD',VariableName,'.png']);

hold off;

%plot(polyfit(repmat(1:time,180,1),squeeze(mean(FluxDeparture,2)),1)*[1;0])

% %these should all be around 0
% plot(mean(squeeze(mean(FluxDeparture,2)),2))
% geoshow(mean(NewNet,3),geoidrefvec,'DisplayType','texturemap');colorbar
end
