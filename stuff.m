load data1.mat

%zonal weights from http://ceres.larc.nasa.gov/data/zone_weights_lou.txt

ncdisp('rlutcs_CERES-EBAF_L4_Ed2-6_200003-201012.nc')

rlutcs = ncread('rlutcs_CERES-EBAF_L4_Ed2-6_200003-201012.nc','rlutcs');
rsutcs = ncread('rsutcs_CERES-EBAF_L4_Ed2-6_200003-201012.nc', 'rsutcs');
rsdt = ncread('rsdt_CERES-EBAF_L4_Ed2-6_200003-201012.nc','rsdt');

net = rsdt - rsutcs - rlutcs;

lat = ncread('rlutcs_CERES-EBAF_L4_Ed2-6_200003-201012.nc','lat');
lon = ncread('rlutcs_CERES-EBAF_L4_Ed2-6_200003-201012.nc','lon');
time = ncread('rlutcs_CERES-EBAF_L4_Ed2-6_200003-201012.nc','time');
% 180*sind(lat)

%Monthly averaged, EBAF product, 3/2000-12/2010, all at TOA: downward shortwave (rsdt), upward shortwave (rsut), upward shortwave in clear-sky (rsutcs), upward longwave (rlut), upward longwave in clear-sky (rlutcs)

%CERES Energy Balanced and Filled (EBAF) Ed2.6r data product.

[LON,LAT]=meshgrid(lon,lat);

dataSum = sum(net,3)/130;
size(dataSum)


%%for processing a single month
data2Plot = net(:,:,1)';
mean(data2Plot,2)
LatWeights(:,2)
size(mean(data2Plot,2).*LatWeights(:,2))
plot(mean(data2Plot,2).*LatWeights(:,2))
plot(sum(data2Plot,2).*LatWeights(:,2),'r')

%%%%%%%
permuteNet = permute(net,[2 1 3]);
colorMarker = {'gx','bx','ro','yd','c^','m-','rs','k<','ch','ko','kx','k-'};
cmap = jet;

%%for processing multiple plots

Year=2001;
MonthIndices = 11+12*(Year-2001):22+12*(Year-2001);



set(gca,'FontSize',20)
for j=MonthIndices
   %plot(LAT(:,2),sum(permuteNet(:,:,j),2).*LatWeights(:,2),colorMarker{mod(j+1:j+1,12)+1})
   plot(LAT(:,2),sum(permuteNet(:,:,j),2).*LatWeights(:,2),'color',cmap((mod(j+1:j+1,12)+1)*5,:),'LineWidth',3)
   latFluxes = sum(permuteNet(:,:,j),2).*LatWeights(:,2);
   difNHSH = [difNHSH sum(latFluxes(91:end)-latFluxes(1:90))];
   %plot(sum(permuteNet(91:end,:,j),2) - sum(flipud(permuteNet(1:90,:,j)),2))
   hold on;   
end
title(num2str(Year));
xlabel('latitude');
ylabel('Total Heat Flux (rsdt - rsutcs - rlutcs)');
set(legend('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'),'Location','BestOutside') %NorthWestOutside

plot(difNHSH)
heatImbalance = sum(sum(permuteNet(:,:,j),2).*LatWeights(:,2))


%%%NH-SH for each 12-month increment

difNHSH = [];

for j=1:130
    latFluxes = sum(permuteNet(:,:,j),2).*LatWeights(:,2);
    difNHSH = [difNHSH sum(latFluxes(91:end)-latFluxes(1:90))];
end


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontSize',20)
windowSize = 12;
AllRunningMeans = filter([31/365 28/365 31/365 30/365 31/365 30/365 31/365 31/365 30/365 31/365 30/365 31/365],1,difNHSH) %but wouldn't this make jan count as april for some starts?
plot(AllRunningMeans(12:end))
grid on;
set(gca,'xtick',[10:12:130])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010})
xlabel('Year where mean ends on');
ylabel('NH-SH Difference in Total Heat Flux (rsdt - rsutcs - rlutcs)');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AllRunningMeans = runmean(difNHSH,12);
% plot(AllRunningMeans)
% grid on;
% set(gca,'xtick',[10:12:130])
% set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010})
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
allMonthWeights = repmat(days_per_month,1,10);
allMonthWeights = [31 30 31 30 31 31 30 31 30 31 allMonthWeights];


moving_sum = @(n, x) filter(ones(1,n), 1, x);
moving_weighted_avg = moving_sum(12, difNHSH .* allMonthWeights) ...
    ./ moving_sum(12, allMonthWeights);

plot(moving_weighted_avg(12:end))

%%%%%%%%%

set(gca,'FontSize',20)
meanNHSH=[];
for jj=1:119
    meanNHSH = [meanNHSH mean(difNHSH(jj:jj+11))];
end
plot(meanNHSH)
grid on;
set(gca,'xtick',[6:12:125])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010})
xlabel('Year where mean is centered on');
ylabel('NH-SH Difference in Total Heat Flux (rsdt - rsutcs - rlutcs)');


%%%%%%%%%%%%%%%%%%%%
Year=2001;
MonthIndices = 11+12*(Year-2001):22+12*(Year-2001);
difNHSH = [];
set(gca,'FontSize',20)
MonthIndicesIncrement = 0; %could change


latFluxes = sum(permuteNet(:,:,j),2).*LatWeights(:,2);

netSum=sum(net(:,:,MonthIndices+MonthIndicesIncrement),3);

%%%%%%%%%%%%%%%%%%

data2Plot = dataSum';
%topHalf = data2Plot(1:180

NorthHalf = data2Plot(91:end,:);
NorthHalf = flipud(NorthHalf);
SouthHalf = data2Plot(1:90,:);

NorthLat = LAT(91:end,:);
SouthLat = LAT(1:90,:);
NorthLat = flipud(NorthLat);
NorthLong = LON(91:end,:);
difHemispheres = data2Plot(91:end,:) - flipud(data2Plot(1:90,:));
NorthMinusSouth = NorthLat - SouthLat;


contourf(NorthLong,NorthLat,difHemispheres);colorbar
%%%
%geoshow(difHemispheres, geoidrefvec, 'DisplayType', 'texturemap');

%%%%HEMISPHERIC DIF MAPS

set(gca,'FontSize',20)
ax = worldmap('World');
setm(ax, 'Origin', [0 180 0])
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
%geoshow(difHemispheres, geoidrefvec, 'DisplayType', 'texturemap');colorbar
%geoshow(LAT(91:end,:) - flipud(LAT(1:90,:)), geoidrefvec, 'DisplayType', 'texturemap');colorbar
geoshow(difHemispheres, geoidrefvec, 'DisplayType', 'texturemap');colorbar
%geoshow(LAT, geoidrefvec, 'DisplayType', 'texturemap');colorbar
%geoshow(LAT(91:end,:), geoidrefvec, 'DisplayType', 'texturemap');colorbar
title('NH-SH Hemispheric Difference of CERES-EBAF TOA (rsdt - rsutcs - rlutcs), averaged over all months 2000-2010')

%%%%%%%%%%%%%



ax = worldmap('World');
setm(ax, 'Origin', [0 180 0])
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(data2Plot, geoidrefvec, 'DisplayType', 'texturemap');

data2Plot

data2Plot(90+45,180+25) %180 by 360 array
data2Plot(90-45,180+25)
contourf(LON,LAT,data2Plot);colorbar

worldmap([-50 50],[160 -30])
load geoid
geoshow(data2Plot, geoidrefvec, 'DisplayType', 'texturemap');
load coast
geoshow(lat, long)

ax = worldmap('World');
setm(ax, 'Origin', [0 180 0])
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(data2Plot, geoidrefvec, 'DisplayType', 'texturemap');colorbar
%lakes = shaperead('worldlakes', 'UseGeoCoords', true);
%geoshow(lakes, 'FaceColor', 'blue')
%rivers = shaperead('worldrivers', 'UseGeoCoords', true);
%geoshow(rivers, 'Color', 'blue')
%cities = shaperead('worldcities', 'UseGeoCoords', true);
%geoshow(cities, 'Marker', '.', 'Color', 'red')
title('CERES-EBAF TOA (rsdt - rsutcs - rlutcs), averaged over all months')

% 
% [7:37:32 PM] Dargan M. W. Frierson: area weight (F) = int f cos(lat) d lat
% [7:38:05 PM] Dargan M. W. Frierson: area weight f = int f cos(lat) d lat
% [7:38:13 PM] Dargan M. W. Frierson: = int f d sin lat
% [7:39:14 PM] Dargan M. W. Frierson: start over
% [7:39:38 PM] Dargan M. W. Frierson: area weight f = int f cos (lat) d lat / (int cos(lat) d lat)
% [7:40:04 PM] Dargan M. W. Frierson: = int f d sin lat / (int d sin lat)
% [7:41:24 PM] Dargan M. W. Frierson: net = SW down - SW up - LW up
% [7:46:17 PM] Dargan M. W. Frierson: 04/2000-03/2001
% [7:46:32 PM] Dargan M. W. Frierson: 05/2000-04/2001
% [7:47:16 PM] Dargan M. W. Frierson: all the way to 06/2012-05/2013
% [7:52:31 PM] *** Call ended, duration 23:0 ***
