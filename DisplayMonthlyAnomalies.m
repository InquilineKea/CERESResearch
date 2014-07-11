function blah = DisplayMonthlyAnomalies(Var1,Var2,Var3)

load LatWeights.mat
%zonal weights from http://ceres.larc.nasa.gov/data/zone_weights_lou.txt
ncdisp('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc');
rlutcs = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','rlutcs');
rsutcs = ncread('rsutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc', 'rsutcs');
rsdt = ncread('rsdt_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','rsdt');
rsut = ncread('rsut_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','rsut');
rlut = ncread('rlut_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','rlut');
SW = rsdt-rsut;
LW = -rlut;
SWCF = rsutcs-rsut; %what about shortwave down forcing? since this is what clouds should reflect...
%rsutcs < rsut. less shortwave up in clearsky? because less shortwave is
%reflected in clearsky..
LWCF = rlutcs-rlut; 
%rlutcs > rlut. more longwave up in clearsky
netcs = rsdt - rsutcs - rlutcs; %isnt this just clear sky??
net = rsdt - rsut - rlut;
lat = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');
lon = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lon');
time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));
RawFlux.Net = permute(net,[2 1 3]);
RawFlux.SW = permute(SW,[2 1 3]);
RawFlux.LW = permute(LW,[2 1 3]);
RawFlux.SWCF = permute(SWCF,[2 1 3]);
RawFlux.LWCF = permute(LWCF,[2 1 3]);
RawFlux.SWclear = permute(rsdt-rsutcs,[2 1 3]);
RawFlux.LWclear = permute(-rlutcs,[2 1 3]);
RawFlux.Clear = permute(netcs,[2 1 3]);
RawFlux.Cloud = RawFlux.SWCF + RawFlux.LWCF;
clearvars net rsut rlut rsutcs rlutcs SWCF LWCF temp LW SW netcs rsdt

temp = ncread('t2m.nc','t2m');
temptime = ncread('t2m.nc','time');
tempLats= ncread('t2m.nc','latitude');
tempLongs = ncread('t2m.nc','longitude');
temp = temp(:,:,find(temptime==878016):find(temptime==991296));
temp = permute(temp,[2 1 3]);
E=zeros(180,360,156);
for depth=1:size(temp,3)
  temp(:,:,depth) = flipud(temp(:,:,depth));
  E(:,:,depth)=imresize(temp(:,:,depth),[180 360]);
end
RawFlux.Temp = E;
clear E;
precip = ncread('precip.mon.mean.nc','precip');
lat = ncread('precip.mon.mean.nc','lat');
long = ncread('precip.mon.mean.nc','lon');
preciptime = ncread('precip.mon.mean.nc','time');
precip = permute(precip,[2 1 3]);
precip = precip(:,:,find(preciptime== 73108):find(preciptime== 77828));
for depth=1:size(precip,3)
  precip(:,:,depth) = flipud(precip(:,:,depth));
end
RawFlux.Precip = precip;
clear precip;
clearvars net rsut rlut rsutcs rlutcs SWCF LWCF temp LW SW netcs rsdt

E=zeros(180,360,156);
for depth=1:size(RawFlux.Precip,3)
  E(:,:,depth) = imresize(RawFlux.Precip(:,:,depth),[180 360]);
end
RawFlux.Precip = E;
clear E;

FluxNames = fieldnames(RawFlux);

% FluxNames = {'Net','SW','LW'};

% MonthFilterSize = 3; %or can be 6, 12
% MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).Net = MovingAverageFilterAllTiles(RawFlux.Net,'Net',3);

MonthFilterSize=1;
for i=1:length(FluxNames)
MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}) = MovingAverageFilterAllTiles(RawFlux.(FluxNames{i}),(FluxNames{i}),MonthFilterSize);
    end

clear RawFlux;

Month1=3;Year1=2000;
Month2=2;Year2=2001;
for k=1:size(MovingAvg.Window1Months.Net,3)
f = figure('Visible','off');
% f = figure('Visible','on')
    MakeLizMap;
colormap(lizmap)
subplot(2,2,1)

set(gca,'FontSize',20)
    load geoid;ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);
    curPlot = MovingAvg.Window1Months.Net(:,:,k);
    geoshow(curPlot,geoidrefvec,'DisplayType','texturemap');colorbar
    load coast;plotm(lat,long,'black')
    MaxBeta = max(max(abs(curPlot)));
%     caxis([-MaxBeta MaxBeta])
caxis([-20 20])
    title(['1-Month Net Anomaly on ' num2str(Year1), '-', num2str(Month1)]);
subplot(2,2,2)
    load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);
    curPlot = MovingAvg.Window1Months.(Var1)(:,:,k);
    %curPlot = MovingAvg.Window1Months.SW(:,:,k);
    geoshow(curPlot,geoidrefvec,'DisplayType','texturemap');colorbar
    load coast;plotm(lat,long,'black')
%     MaxBeta = max(max(abs(curPlot)));
%     caxis([-MaxBeta MaxBeta])
     caxis([-20 20])

title(Var1)

subplot(2,2,3)
    load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);
    curPlot = MovingAvg.Window1Months.(Var2)(:,:,k);
    geoshow(curPlot,geoidrefvec,'DisplayType','texturemap');colorbar
    load coast;plotm(lat,long,'black')
%     MaxBeta = max(max(abs(curPlot)));
%     caxis([-MaxBeta MaxBeta])
     caxis([-20 20])

title(Var2)

subplot(2,2,4)
    load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);
    curPlot = MovingAvg.Window1Months.(Var3)(:,:,k);
    geoshow(curPlot,geoidrefvec,'DisplayType','texturemap');colorbar
    load coast;plotm(lat,long,'black')
    MaxBeta = max(max(abs(curPlot)));
     caxis([-4 4])
title(Var3)
 subplots = get(gcf,'Children');
 AllPositions= get(subplots,'Position');

set(subplots(2),'Position',get(subplots(2),'OuterPosition'))
set(subplots(4),'Position',get(subplots(4),'OuterPosition'))
set(subplots(6),'Position',get(subplots(6),'OuterPosition'))
set(subplots(8),'Position',get(subplots(8),'OuterPosition'))

grid on;
% tightfig
set(gcf,'Renderer','painters')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Net',Var1,Var2,Var3,'Anomaly',num2str(k),'.png']);

% 
%     currFrame = getframe;
%     writeVideo(Net12MnthVid,currFrame);
    hold off;
    Month2 = Month2+1;Month1 = Month1+1;
    if eq(Month2,13)
        Month2 = 1;
    Year2 = Year2+1;
    end
    if eq(Month1,13)
        Month1=1;
        Year1=Year1+1;
    end
end
end