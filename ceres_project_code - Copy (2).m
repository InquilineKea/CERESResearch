load LatWeights.mat

%zonal weights from http://ceres.larc.nasa.gov/data/zone_weights_lou.txt

ncdisp('CERES_EBAF-TOA_Ed2.7_Subset_200003-201302.nc');

rlutcs = ncread('CERES_EBAF-TOA_Ed2.7_Subset_200003-201302.nc','toa_lw_clr_mon');
rsutcs = ncread('CERES_EBAF-TOA_Ed2.7_Subset_200003-201302.nc', 'toa_sw_clr_mon');
%rsdt = ncread('rsdt_CERES-EBAF_L4_Ed2-6_200003-201012.nc','rsdt');
rsut = ncread('CERES_EBAF-TOA_Ed2.7_Subset_200003-201302.nc','toa_sw_all_mon');
%this is toa_outgoing_shortwave_flux... (which the satellite measures..)
%what about the incoming SW flux?
rlut = ncread('CERES_EBAF-TOA_Ed2.7_Subset_200003-201302.nc','toa_lw_all_mon');

SWCF = rsutcs-rsut; %
LWCF = rlut-rlutcs; %
net = ncread('CERES_EBAF-TOA_Ed2.7_Subset_200003-201302.nc','toa_net_all_mon');

lat = ncread('CERES_EBAF-TOA_Ed2.7_Subset_200003-201302.nc','lat');
lon = ncread('CERES_EBAF-TOA_Ed2.7_Subset_200003-201302.nc','lon');
time = length(ncread('CERES_EBAF-TOA_Ed2.7_Subset_200003-201302.nc','time'));
% 180*sind(lat)

%Monthly averaged, EBAF product, 3/2000-12/2010, all at TOA: downward shortwave (rsdt), upward shortwave (rsut), upward shortwave in clear-sky (rsutcs), upward longwave (rlut), upward longwave in clear-sky (rlutcs)

%CERES Energy Balanced and Filled (EBAF) Ed2.6r data product.

[LON,LAT]=meshgrid(lon,lat);
dataSum = sum(net,3)/time;
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
SW = permute(rsut,[2 1 3]);
LW = permute(rlut,[2 1 3]);
SWclear = permute(rsutcs,[2 1 3]);
LWclear = permute(rlutcs,[2 1 3]);
clearvars net SWCF LWCF rsut rlut rsutcs rlutcs

colorMarker = {'gx','bx','ro','yd','c^','m-','rs','k<','ch','ko','kx','k-'};
cmap = jet;

%%for processing multiple plots
Year=2001;
latFluxesNet = LatWeightedAverage(permuteNet,Year,'Net');
latSWCF = LatWeightedAverage(permuteSWCF,Year,'SWCF (rsutcs-rsut)');
latLWCF = LatWeightedAverage(permuteLWCF, Year,'LWCF (rlut-rlutcs)');
latSW = LatWeightedAverage(SW, Year,'SW');
latLW = LatWeightedAverage(LW, Year,'LW');
latSWClear = LatWeightedAverage(LW, Year,'SWClear');
latLWClear = LatWeightedAverage(LW, Year,'LWClear');


%%%%%%%%%%%%

%plot(IDR)
heatImbalance = sum(sum(permuteNet(:,:,j),2).*LatWeights(:,2))

%%%NH-SH for each 12-month increment

[IDR,netNH,netSH,net_MWA,MonthlySubtractedIDR] = NHMinusSH(permuteNet,time);
title('Running Mean in Net, after mean-subtraction')
set(gcf,'paperposition',[1 1 24 12])
print(gcf,'-dpng','-r300',['RunningMean_','net','_','.png']);
[difSWCF,SWCF_NH,SWCF_SH,SWCF_MWA,MonthlySubtractedSWCF] = NHMinusSH(permuteSWCF,time);
title('Running Mean in SWCF, after mean-subtraction') %does this really matter?
set(gcf,'paperposition',[1 1 24 12])
print(gcf,'-dpng','-r300',['RunningMean_','SWCF','_','.png']);
[difLWCF,LWCF_NH,LWCF_SH,LWCF_MWA,MonthlySubtractedLWCF] = NHMinusSH(permuteLWCF,time);
title('Running Mean in LWCF, after mean-subtraction')%does this really matter?
set(gcf,'paperposition',[1 1 24 12])
print(gcf,'-dpng','-r300',['RunningMean_','LWCR','_','.png']);
[IDS,SW_NH,SW_SH,SW_MWA,MonthlySubtractedIDS] = NHMinusSH(SW,time);
title('Running Mean in SW, after mean-subtraction')
set(gcf,'paperposition',[1 1 24 12])
print(gcf,'-dpng','-r300',['RunningMean_','SW','_','.png']);
[IDL,LW_NH,LW_SH,LW_MWA,MonthlySubtractedIDL] = NHMinusSH(LW,time);
title('Running Mean in LW, after mean-subtraction')
set(gcf,'paperposition',[1 1 24 12])
print(gcf,'-dpng','-r300',['RunningMean_','LW','_','.png']);

hold on;
grid on;
set(gca,'FontSize',20)
plot(MonthlySubtractedSWCF(3,:),'LineWidth',3)
plot(MonthlySubtractedLWCF(3,:),'g','LineWidth',3)
plot(MonthlySubtractedIDR(1,:),'r','LineWidth',3)
set(gca,'xtick',11:12:156)
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year');
ylabel('Flux Difference');
set(legend('SH MonthlySubtractedSWCF','SH MonthlySubtractedLWCF','MonthlySubtractedIDR'),'Location','BestOutside') %NorthWestOutside

hold on;
grid on;
set(gca,'FontSize',20)
plot(MonthlySubtractedIDR(2,:),'LineWidth',3)
plot(MonthlySubtractedIDR(3,:),'g','LineWidth',3)
plot(MonthlySubtractedIDR(1,:),'r','LineWidth',3)
set(gca,'xtick',11:12:156)
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year');
ylabel('Flux Difference');
set(legend('NH MonthlySubtracted Net','SH MonthlySubtracted Net','MonthlySubtractedIDR'),'Location','BestOutside') %NorthWestOutside


corrcoef(MonthlySubtractedIDR(1,:),MonthlySubtractedSWCF(2,:))
corrcoef(MonthlySubtractedIDR(1,:),MonthlySubtractedSWCF(3,:))
corrcoef(MonthlySubtractedIDR(1,:),MonthlySubtractedLWCF(2,:))
corrcoef(MonthlySubtractedIDR(1,:),MonthlySubtractedLWCF(3,:))


%IDR_MWA: R1 = IDR. R2 = NH, R3 = SH

hold on;
grid on;
set(gca,'FontSize',20)
plot(LWCF_MWA(2,:)-mean(LWCF_MWA(2,:)),'LineWidth',3)
plot(net_MWA(3,:)-mean(net_MWA(3,:)),'r','LineWidth',3)
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year end');

corrcoef(LWCF_MWA(2,:),net_MWA(3,:))

%meanNHSH is 1x119

norm(valid_moving_weighted_avg-meanNHSH)
norm(valid_moving_weighted_avg-meanNHSH)/norm(valid_moving_weighted_avg)
%6 percent error

%%%%%%%%%%%%%%%%%

netVar = zeros(3,12);
hold on;
for i=1:12
    netVar(1,i) = MonthTrend(IDR,i);
    netVar(2,i) = CalcMonthVar(netNH,i);
    netVar(3,i) = CalcMonthVar(netSH,i);
end
set(legend('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'),'Location','BestOutside')
title('NH-SH Difference in Net Flux of a Given Month, Subtracted by Monthly Mean Over All Years. Jan, Feb offset by 1');

%can also do this for NH,SH individually

SWCFVar = zeros(3,12);
hold on;
for i=1:12
    SWCFVar(1,i) = MonthTrend(difSWCF,i);
    SWCFVar(2,i) = CalcMonthVar(SWCF_NH,i);
    SWCFVar(3,i) = CalcMonthVar(SWCF_SH,i);
end
set(legend('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'),'Location','BestOutside')
title('NH-SH Difference in Net SWCF of a Given Month, Subtracted by Monthly Mean Over All Years. Jan, Feb offset by 1');

LWCFVar = zeros(3,12);
hold on;
for i=1:12
    LWCFVar(1,i) = MonthTrend(difLWCF,i);
    LWCFVar(2,i) = CalcMonthVar(LWCF_NH,i);
    LWCFVar(3,i) = CalcMonthVar(LWCF_SH,i);
end
set(legend('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'),'Location','BestOutside')
title('NH-SH Difference in Net LWCF of a Given Month, Subtracted by Monthly Mean Over All Years. Jan, Feb offset by 1');

LWVar = zeros(3,12);
hold on;
for i=1:12
    LWVar(1,i) = MonthTrend(IDL,i);
    LWVar(2,i) = CalcMonthVar(LW_NH,i);
    LWVar(3,i) = CalcMonthVar(LW_SH,i);
end
set(legend('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'),'Location','BestOutside')
title('NH-SH Difference in LW of a Given Month, Subtracted by Monthly Mean Over All Years. Jan, Feb offset by 1');

SWVar = zeros(3,12);
hold on;
for i=1:12
    SWVar(1,i) = MonthTrend(IDS,i);
    SWVar(2,i) = CalcMonthVar(SW_NH,i);
    SWVar(3,i) = CalcMonthVar(SW_SH,i);
end
set(legend('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'),'Location','BestOutside')
title('NH-SH Difference in SW of a Given Month, Subtracted by Monthly Mean Over All Years. Jan, Feb offset by 1');

hold on;
grid on;
set(gca,'FontSize',20)
plot(netVar(1,:),'LineWidth',3)
plot(SWCFVar(1,:),'g','LineWidth',3)
plot(LWCFVar(1,:),'r','LineWidth',3)
plot(LWVar(1,:),'k','LineWidth',3)
plot(SWVar(1,:),'m','LineWidth',3)
set(legend('net','SWCF','LWCF','LW','SW'),'Location','BestOutside')
title('Variances of each NH-SH component by month')
set(gca,'xtick',1:12)

PlotDecomposedVar(netVar)
PlotDecomposedVar(SWCFVar)
PlotDecomposedVar(LWCFVar)
PlotDecomposedVar(LWVar)
PlotDecomposedVar(SWVar)

style = hgexport('factorystyle');
style.Bounds = 'tight';
hgexport(gcf,'-clipboard',style,'applystyle', true);
%hgexport(gcf, 'figure1.jpg', hgexport('factorystyle'), 'Format', 'jpeg');



%%%%%%%%%%%%%%%

GoodardSeptSeaIceExtent = [6.32 6.75 5.96 6.15 6.05 5.57 5.92 4.3 4.73 5.39 4.93 4.63 3.63];
Years = [2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012];

%subplot(1,2,1)
set(gca,'FontSize',20)
%1 is march 2000. leap years in february 2004, 2008
plotyy(Years,IDR(7:12:end),Years,GoodardSeptSeaIceExtent)
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
IDR = [];
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
