load LatWeights.mat

%zonal weights from http://ceres.larc.nasa.gov/data/zone_weights_lou.txt

ncdisp('rlutcs_CERES-EBAF_L4_Ed2-6_200003-201012.nc');

rlutcs = ncread('rlutcs_CERES-EBAF_L4_Ed2-6_200003-201012.nc','rlutcs');
rsutcs = ncread('rsutcs_CERES-EBAF_L4_Ed2-6_200003-201012.nc', 'rsutcs');
rsdt = ncread('rsdt_CERES-EBAF_L4_Ed2-6_200003-201012.nc','rsdt');
rsut = ncread('rsut_CERES-EBAF_L4_Ed2-6_200003-201012.nc','rsut');
rlut = ncread('rlut_CERES-EBAF_L4_Ed2-6_200003-201012.nc','rlut');

SWCF = rsutcs-rsut; %
LWCF = rlut-rlutcs; %
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
permuteSWCF = permute(SWCF,[2 1 3]);
permuteLWCF = permute(LWCF,[2 1 3]);

colorMarker = {'gx','bx','ro','yd','c^','m-','rs','k<','ch','ko','kx','k-'};
cmap = jet;

%%for processing multiple plots
Year=2001;
latFluxesNet = LatWeightedAverage(permuteNet,Year);
latSWCF = LatWeightedAverage(permuteSWCF,Year);
lawLWCF = LatWeightedAverage(permuteLWCF, Year);

%%%%%%%%%%%%

plot(difNHSH)
heatImbalance = sum(sum(permuteNet(:,:,j),2).*LatWeights(:,2))

%%%NH-SH for each 12-month increment

NHMinusSH(permuteNet)
NHMinusSH(permuteSWCF)
NHMinusSH(permuteLWCF)

%meanNHSH is 1x119

norm(valid_moving_weighted_avg-meanNHSH)
norm(valid_moving_weighted_avg-meanNHSH)/norm(valid_moving_weighted_avg)
%6 percent error

%%%%%%%%%%%%%%%

GoodardSeptSeaIceExtent = [6.32 6.75 5.96 6.15 6.05 5.57 5.92 4.3 4.73 5.39 4.93];
Years = [2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010];

%subplot(1,2,1)
set(gca,'FontSize',20)
%1 is march 2000. leap years in february 2004, 2008
plotyy(Years,difNHSH(7:12:end),Years,GoodardSeptSeaIceExtent)
grid on;
%set(gca,'XTickLabel',{2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010})
xlabel('Year');
title('NH-SH Difference in September Heat Flux');
ylabel('Flux Difference');
set(legend('CERES Flux Difference','Goodard Sea Ice Extent'),'Location','BestOutside') 

%axes(ax(1)); ylabel('First y-label');
%axes(ax(2)); ylabel('Goodard Sea Ice Extent (value)');

%corr(difNHSH(7:12:end),GoodardSeptSeaIceExtent)
%subplot(1,2,2)
%plot([6.32 6.75 5.96 6.15 6.05 5.57 5.92 4.3 4.73 5.39 4.93],'bx')



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
