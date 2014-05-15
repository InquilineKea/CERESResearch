load LatWeights.mat

%zonal weights from http://ceres.larc.nasa.gov/data/zone_weights_lou.txt

ncdisp('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc');

%cycle through netcdf fields, initialize structure, set structure name to
%fluxstruct.name = ncread('','FIELDNAME)
%fluxstruct.flux = ncread('',FLUX)
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

%also why no rsdtcs? isnt that more relevant than rsutcs?

lat = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');
lon = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lon');
time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));

% 180*sind(lat)

%Monthly averaged, EBAF product, 3/2000-12/2010, all at TOA: downward shortwave (rsdt), upward shortwave (rsut), upward shortwave in clear-sky (rsutcs), upward longwave (rlut), upward longwave in clear-sky (rlutcs)

%CERES Energy Balanced and Filled (EBAF) Ed2.6r data product.
% 
% [LON,LAT]=meshgrid(lon,lat);
% dataSum = sum(net,3)/time;
% size(dataSum);
% 
% %%for processing a single month
% data2Plot = net(:,:,1)';
% mean(data2Plot,2);
% LatWeights(:,2);
% size(mean(data2Plot,2).*LatWeights(:,2))
% plot(mean(data2Plot,2).*LatWeights(:,2))
% plot(sum(data2Plot,2).*LatWeights(:,2),'r')

%%%%%%%
Net = permute(net,[2 1 3]);
SWCF = permute(SWCF,[2 1 3]);
LWCF = permute(LWCF,[2 1 3]);
SW = permute(SW,[2 1 3]);
LW = permute(LW,[2 1 3]);
SWclear = permute(rsdt-rsutcs,[2 1 3]);
LWclear = permute(-rlutcs,[2 1 3]);
netclear = permute(netcs,[2 1 3]);
TotalCloudForcing = SWCF + LWCF;
clearvars net rsut rlut rsutcs rlutcs

colorMarker = {'gx','bx','ro','yd','c^','m-','rs','k<','ch','ko','kx','k-'};
cmap = jet;

%%for processing multiple plots
% Year=2001;
% %LatWeightedAverage(mean(permuteNet,3),Year,'Net (rsdt - rsut - rlut)');
% LatWeightedAverage(bsxfun(@times, Net,topo<0),Year,'Net (rsdt - rsut - rlut)');%ocean only
% LatWeightedAverage(SWCF,Year,'SWCF (rsutcs-rsut)');
% LatWeightedAverage(LWCF, Year,'LWCF (rlutcs-rlut)');
% LatWeightedAverage(SW, Year,'SW');
% LatWeightedAverage(LW, Year,'LW');
% LatWeightedAverage(SWclear, Year,'SWClear (rsdt-rsutcs)');
% LatWeightedAverage(LWclear, Year,'LWClear (-rlutcs)');
% LatWeightedAverage(netclear, Year,'NetClear (rsdt - rsutcs - rlutcs)');

%%%Land
LatWeightedAverage(bsxfun(@times, SW,topo>=0),Year,'SWLand');%land only
hold off;
LatWeightedAverage(bsxfun(@times, SWCF,topo>=0),Year,'SWCFLand');%land only
hold off;
LatWeightedAverage(bsxfun(@times, LWCF,topo>=0),Year,'LWCFLand');%land only
hold off;
LatWeightedAverage(bsxfun(@times, LW,topo>=0),Year,'LWLand');%land only
hold off;
LatWeightedAverage(bsxfun(@times, Net,topo>=0),Year,'NetLand');%land only
hold off;
LatWeightedAverage(bsxfun(@times, SWclear,topo>=0),Year,'SWClearLand');%land only
hold off;
LatWeightedAverage(bsxfun(@times, LWclear,topo>=0),Year,'LWClearLand');%land only
hold off;
LatWeightedAverage(bsxfun(@times, netclear,topo>=0),Year,'NetClearLand');%land only

%%%Ocean
LatWeightedAverage(bsxfun(@times, SW,topo<0),Year,'SWOcean');%Ocean only
hold off;
LatWeightedAverage(bsxfun(@times, SWCF,topo<0),Year,'SWCFOcean');%Ocean only
hold off;
LatWeightedAverage(bsxfun(@times, LWCF,topo<0),Year,'LWCFOcean');%Ocean only
hold off;
LatWeightedAverage(bsxfun(@times, LW,topo<0),Year,'LWOcean');%Ocean only
hold off;
LatWeightedAverage(bsxfun(@times, Net,topo<0),Year,'NetOcean');%Ocean only
hold off;
LatWeightedAverage(bsxfun(@times, SWclear,topo<0),Year,'SWClearOcean');%land only
hold off;
LatWeightedAverage(bsxfun(@times, LWclear,topo<0),Year,'LWClearOcean');%land only
hold off;
LatWeightedAverage(bsxfun(@times, netclear,topo<0),Year,'NetClearOcean');%land only


%%%%%%%%%%%%

%plot(IDR)
heatImbalance = sum(sum(Net(:,:,j),2).*LatWeights(:,2))

%%%NH-SH for each 12-month increment


plot(xcorr(net_MWA(2,:),net_MWA(3,:)))
plot(xcorr(net_MWA(1,:),net_MWA(3,:)))

[IDR,netNH,netSH,net_MWA,MonthlySubtractedIDR] = NHMinusSH(Net,time,0,90);
title('Running Mean in Net')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','net.png']);

NHMinusSH(Net,time,30,90)
title('30deg-90deg Running Mean in net')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','net30_90.png']);
NHMinusSH(Net,time,0,30)
title('0-30deg Running Mean in net')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','net0_30.png']);

% varfun = @(func, varargin)cellfun(...
%     @(x)evalin('base', sprintf('%s(%s, y)', func, x)), ...
%     arrayfun(@inputname, 2:nargin, 'UniformOutput', false), ...
%     'UniformOutput', false);
% 
% varfun(NHMinusSH(netclear,time,0,90))

%%DO THIS NEXT FOR DIVIDING UP FLUXES INTO 4:
%%0 to 14.5, 14.5 to 30, 30 to 48.59, 48.59 to 90

[IDRcs,netNHcs,netSHcs,net_MWAcs,MonthlySubtractedIDRcs] = NHMinusSH(netclear,time,0,90);
title('Running Mean in Net Clear-Sky')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','netCS.png']);

[difSWCF,SWCF_NH,SWCF_SH,SWCF_MWA,MonthlySubtractedSWCF] = NHMinusSH(SWCF,time,0,90);
title('Running Mean in SWCF') %does this really matter?
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SWCF.png']);

[difLWCF,LWCF_NH,LWCF_SH,LWCF_MWA,MonthlySubtractedLWCF] = NHMinusSH(LWCF,time,0,90);
title('Running Mean in LWCF')%does this really matter?
grid on;
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LWCF.png']);

[IDS,SW_NH,SW_SH,SW_MWA,MonthlySubtractedIDS] = NHMinusSH(SW,time,0,90);
title('Running Mean in SW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SW.png']);

NHMinusSH(SW,time,30,90)
title('30deg-90deg Running Mean in SW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SW30_90.png']);
NHMinusSH(SW,time,0,30)
title('0-30deg Running Mean in SW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SW0_30.png']);

[IDL,LW_NH,LW_SH,LW_MWA,MonthlySubtractedIDL] = NHMinusSH(LW,time,0,90);
title('Running Mean in LW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LW.png']);

NHMinusSH(LW,time,30,90)
title('30deg-90deg Running Mean in LW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LW30_90.png']);
NHMinusSH(LW,time,0,30)
title('0-30deg Running Mean in LW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LW0_30.png']);

[IDLcs,LW_NHcs,LW_SHcs,LW_MWAcs,MonthlySubtractedIDLcs] = NHMinusSH(LWclear,time,0,90);
title('Running Mean in LW Clear')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LWclear.png']);

[IDScs,SW_NHcs,SW_SHcs,SW_MWAcs,MonthlySubtractedIDScs] = NHMinusSH(SWclear,time,0,90);
title('Running Mean in SW Clear')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SWclear.png']);

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
%%%

load('topo.mat','topo','topomap1');
geoshow(topo, geoidrefvec, 'DisplayType', 'texturemap');colorbar
geoshow(topo.*(topo < 0), geoidrefvec, 'DisplayType', 'texturemap');colorbar
geoshow(topo > 0, geoidrefvec, 'DisplayType', 'texturemap');colorbar



topo(topo > 0)

load geoid
subplot(2,1,1)
geoshow(std(Net,0,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
subplot(2,1,2)
%geoshow(std(permuteNet,0,3)-repmat(mean(std(permuteNet,0,3)),180,1), geoidrefvec, 'DisplayType', 'texturemap');colorbar

%cellfun(@(x) geoshow(),)_)

geoshow(mean(Net,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(Net,0,3)-repmat(mean(std(Net,0,3),2),1,360), geoidrefvec, 'DisplayType', 'texturemap');colorbar
title('Net Flux Standard Deviation, after subtracting out zonally-averaged standard deviation (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Globe_','net_Var.png']);

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(netclear,0,3)-repmat(mean(std(netclear,0,3),2),1,360), geoidrefvec, 'DisplayType', 'texturemap');colorbar
title('Net Clear-Sky Standard Deviation, after subtracting out zonally-averaged standard deviation (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Globe_','netCS_Var.png']);

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(SWCF,0,3)-repmat(mean(std(SWCF,0,3),2),1,360), geoidrefvec, 'DisplayType', 'texturemap');colorbar
title('SWCF Standard Deviation, after subtracting out zonally-averaged standard deviation (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Globe_','SWCF_Var.png']);

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(LWCF,0,3)-repmat(mean(std(LWCF,0,3),2),1,360), geoidrefvec, 'DisplayType', 'texturemap');colorbar
title('LWCF Standard Deviation, after subtracting out zonally-averaged standard deviation (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Globe_','LWCF_Var.png']);

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(SW,0,3)-repmat(mean(std(SW,0,3),2),1,360), geoidrefvec, 'DisplayType', 'texturemap');colorbar
title('SW Standard Deviation, after subtracting out zonally-averaged standard deviation (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Globe_','SW_Var.png']);

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(SWclear,0,3)-repmat(mean(std(SWclear,0,3),2),1,360), geoidrefvec, 'DisplayType', 'texturemap');colorbar
title('SW Clear Standard Deviation, after subtracting out zonally-averaged standard deviation (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Globe_','SWClear_Var.png']);

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(LW,0,3)-repmat(mean(std(LW,0,3),2),1,360), geoidrefvec, 'DisplayType', 'texturemap');colorbar
title('LW Standard Deviation, after subtracting out zonally-averaged standard deviation (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Globe_','LW_Var.png']);

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(std(LWclear,0,3)-repmat(mean(std(LWclear,0,3),2),1,360), geoidrefvec, 'DisplayType', 'texturemap');colorbar
title('LW Clear Standard Deviation, after subtracting out zonally-averaged standard deviation (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Globe_','LWClear_Var.png']);

geoshow(std(Net,0,3)-repmat(mean(std(Net,0,3)),180,1), geoidrefvec, 'DisplayType', 'texturemap');colorbar

%sep 2007
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(Net(:,:,91)-mean(Net(:,:,7:12:end),3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(SWclear(:,:,91)-mean(SWclear(:,:,7:12:end),3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

GlobePlot(Net,0)

geoshow(mean(Net(:,:,7:12:end),3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

geoshow(mean(Net(:,:,7:12:end),3)-mean(Net,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

%January 2009 minus other January anomalies
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true); geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(mean(Net(:,:,107),3)-mean(Net(:,:,11:12:156),3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

%maybe just abnormally low SH dif from jan 08 to jan 09?
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true); geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(mean(Net(:,:,107-12),3)-mean(Net(:,:,107),3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true); geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(mean(Net(:,:,108),3)-mean(Net(:,:,12:12:156),3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true); geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(mean(Net(:,:,107-11:107),3)-mean(Net,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar

geoshow(mean(SWCF,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
geoshow(mean(LWCF,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
geoshow(mean(LW,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
geoshow(mean(SW,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
geoshow(mean(LWclear,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
geoshow(mean(SWclear,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
geoshow(mean(netclear,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar



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
%hgexport(gcf, 'figure1.png', hgexport('factorystyle'), 'Format', 'jpeg');



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


latFluxes = sum(Net(:,:,j),2).*LatWeights(:,2);

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

load geoid
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
