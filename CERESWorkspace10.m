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
RawFlux.NetClear = permute(netcs,[2 1 3]);
RawFlux.NetCloud = RawFlux.SWCF + RawFlux.LWCF;
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

for MonthFilterSize =[3,6,12]
    for i=1:length(FluxNames)
MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}) = MovingAverageFilterAllTiles(RawFlux.(FluxNames{i}),(FluxNames{i}),MonthFilterSize);
    end
end
clear RawFlux;

% IndicesMWA.NAM= bsxfun(@rdivide,moving_sum(MonthFilterSize, Indices.NAM .* allMonthWeights) ...
%     ,moving_sum(MonthFilterSize,allMonthWeights)); %or Mean Lat Flux here

for MonthFilterSize=[3,6,12]
IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']) = IndexMovingAverage(MonthFilterSize);
end
SampleWindow = fieldnames(IndicesMvgAvg);
IndexNames = fieldnames(IndicesMvgAvg.(SampleWindow{1}));

%%%for MonthFilterSize = [3,6]
%replace with for MonthFilterSize = 12
% IndexMovingAverage(12) %even for 12 months this works!!

% IndicesMvgAvg.Window3Months.(IndexNames{1})
%for MonthFilterSize=[3,6]
% for MonthFilterSize=[6]
%     for i=1:3
%        for j=1:3
%         [BetaMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{j}).(FluxNames{i}), blah, blah, blah,AlphaMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{j}).(FluxNames{i})]...
%             = Regress1DTimeSeriesMapWithTTestContours(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}) ,...
%         IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{j}),[FluxNames{i}, ' ', num2str(MonthFilterSize),'-Month Mvg Avg'],[IndexNames{j}, ' ', num2str(MonthFilterSize),'-Month Mvg Avg']);
%         close all;
%         MWATimeLatContourRegression(AlphaMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{j}).(FluxNames{i}),BetaMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{j}).(FluxNames{i}),...
%         MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}) ,IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{i}),...
%             [FluxNames{i}],[IndexNames{j}],MonthFilterSize);
%         close all;
%        end
%     end
% end

%%THE SAME

close all;
for MonthFilterSize = [3,6]
    for i=1:3
       for j=1:size(FluxNames)
        j=5
        [Beta.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{i}).(FluxNames{j}), blah, FluxSumAllLats.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{i}).(FluxNames{j}), blah,...
            Alpha.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{i}).(FluxNames{j})] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{j}) ,...
        IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{i}),[FluxNames{j}, ' ', num2str(MonthFilterSize),'-Month Mvg Avg'],[IndexNames{i}, ' ', num2str(MonthFilterSize),'-Month Mvg Avg']);
        close all;
%         RegressedFlux.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{i}).(FluxNames{j})=...
            MWATimeLatContourRegression(Alpha.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{i}).(FluxNames{j})...
            ,Beta.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{i}).(FluxNames{j}),MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{j}) ,...
            IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{i}),...
            [FluxNames{j}],[IndexNames{i}],MonthFilterSize);
        close all;
        end
    end
end

%Time Series of each latitudinal band

% LowerLat = [0,15,30,45];
% UpperLat = [15,30,45,90];
PrecipAsymmetry = MovingAvgTimeSeries.Window12Months.Precip.HemDif0to14./MovingAvgTimeSeries.Window12Months.Precip.Global0to14;
PrecipAsymmetry = MovingAvgTimeSeries.Window12Months.Precip.HemDif0to90./MovingAvgTimeSeries.Window12Months.Precip.Global0to90;

plot(MovingAvgTimeSeries.Window12Months.Precip.Global0to90)
hold on;
plot(MovingAvgTimeSeries.Window12Months.Precip.HemDif0to90,'r')

for MonthFilterSize=[3,6]
    for i=1:size(FluxNames)
    [MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).NH0to90,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).SH0to90,...
    MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).HemDif0to90,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).Global0to90,...
    MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).AsymIndex0to90] ...
      = FluxInLatitudinalBand(RawFlux.(FluxNames{i}),0,90,FluxNames{i},MonthFilterSize);
    [MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).NH0to20,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).SH0to20,...
    MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).HemDif0to20,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).Global0to20,...
    MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).AsymIndex0to20] ...
      = FluxInLatitudinalBand(RawFlux.(FluxNames{i}),0,20,FluxNames{i},MonthFilterSize);
      [MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).NH20to90,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).SH20to90,...
    MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).HemDif20to90,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).Global20to90,...
    MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).AsymIndex20to90] ...
      = FluxInLatitudinalBand(RawFlux.(FluxNames{i}),20,90,FluxNames{i},MonthFilterSize);
%= FluxInLatitudinalBand(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}),0,90,FluxNames{i},MonthFilterSize);
%  hold on; plot(MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).AsymIndex0to90)
    for k=0:3
    [MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).(['NH',num2str(round(rad2deg(asin(k*0.25)))),'to',num2str(round(rad2deg(asin((k+1)*0.25))))]),...
        MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).(['SH',num2str(round(rad2deg(asin(k*0.25)))),'to',num2str(round(rad2deg(asin((k+1)*0.25))))]),...
        MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).(['HemDif',num2str(round(rad2deg(asin(k*0.25)))),'to',num2str(round(rad2deg(asin((k+1)*0.25))))]),...
        MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).(['Global',num2str(round(rad2deg(asin(k*0.25)))),'to',num2str(round(rad2deg(asin((k+1)*0.25))))]),...
        MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}).(['AsymIndex',num2str(round(rad2deg(asin(k*0.25)))),'to',num2str(round(rad2deg(asin((k+1)*0.25))))])] ...
        = FluxInLatitudinalBand(RawFlux.(FluxNames{i}),rad2deg(asin(k*0.25)),rad2deg(asin((k+1)*0.25)),FluxNames{i},MonthFilterSize);
    %= FluxInLatitudinalBand(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}),rad2deg(asin(k*0.25)),rad2deg(asin((k+1)*0.25)),FluxNames{i},MonthFilterSize);
    end
    end
end

% save('2014-04-30.mat','MovingAvgTimeSeries')

%%%this is where getting data for ALL the fluxes can be helpful. 
% 
% PlotLatRegressCoefs(MovingAvg.Window3Months,MovingAvgTimeSeries.Window3Months.Net,FluxNames,3,'Net')
% PlotLatRegressCoefs(MovingAvg.Window6Months,MovingAvgTimeSeries.Window6Months.Net,FluxNames,6,'Net')
% 
% PlotLatRegressCoefs(MovingAvg.Window3Months,MovingAvgTimeSeries.Window3Months.SW,FluxNames,3,'SW')
% PlotLatRegressCoefs(MovingAvg.Window6Months,MovingAvgTimeSeries.Window6Months.SW,FluxNames,6,'SW')
% 
% PlotLatRegressCoefs(MovingAvg.Window3Months,MovingAvgTimeSeries.Window3Months.LW,FluxNames,3,'LW')
% PlotLatRegressCoefs(MovingAvg.Window6Months,MovingAvgTimeSeries.Window6Months.LW,FluxNames,6,'LW')

PlotLatRegressCoefs(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']),IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']),FluxNames,MonthFilterSize,[''])
FluxNameToCompareWith = fieldnames(MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']));

PlotLatRSquared(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']),IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']),FluxNames,MonthFilterSize,[''])

FluxTupleSelections{1} = {'Net';'SWCF';'LWCF'};
FluxTupleSelections{2} = {'LW';'LWclear';'LWCF'};

FluxTupleSelections{1} = {'SWclear';'LWclear';'SWCF';'LWCF';'Net'};

FluxTupleSelections{1} = {'Temp';'Precip';'NetClear';'NetCloud';'Net'};
for jj = 1:length(FluxTupleSelections)
SubsetFluxIndices = [];
for i = 1:length(FluxTupleSelections{jj})
    SubsetFluxIndices = [SubsetFluxIndices find(strcmp(FluxTupleSelections{jj}(i),FluxNames))];
end
SubsetFluxNames = FluxNames(SubsetFluxIndices);
for h = 1:size(SubsetFluxNames)
    i = SubsetFluxIndices(h);
    SubsetFlux.(FluxNames{i}) = MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i});
end
% PlotLatRSquared(SubsetFlux,IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']),SubsetFluxNames,MonthFilterSize,[SubsetFluxNames{:}])
CompareFluxLatTimeSeries(MovingAvgTimeSeries,SubsetFluxNames,SubsetFluxIndices,FluxNames,MonthFilterSize,0,90,'HemDif')
CompareFluxLatTimeSeries(MovingAvgTimeSeries,SubsetFluxNames,SubsetFluxIndices,FluxNames,MonthFilterSize,0,14,'HemDif')
CompareFluxLatTimeSeries(MovingAvgTimeSeries,SubsetFluxNames,SubsetFluxIndices,FluxNames,MonthFilterSize,14,30,'HemDif')
CompareFluxLatTimeSeries(MovingAvgTimeSeries,SubsetFluxNames,SubsetFluxIndices,FluxNames,MonthFilterSize,30,49,'HemDif')
CompareFluxLatTimeSeries(MovingAvgTimeSeries,SubsetFluxNames,SubsetFluxIndices,FluxNames,MonthFilterSize,49,90,'HemDif')

clearvars SubsetFlux SubsetFluxNames SubsetFluxIndices
end


PlotLatRSquared(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']),MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).NetClear,FluxNames,MonthFilterSize,'NetClear')
PlotLatRSquared(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']),MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Temp,FluxNames,MonthFilterSize,'Temperature')



%%%%%%%%%
%2014-05-05
MovingAvgTimeSeries.Window12Months.Net.GlobalAvg=(MovingAvgTimeSeries.Window12Months.Net.Global0to90)/2;
TempFluxNames = fieldnames(MovingAvg.Window12Months);

for i=1:length(TempFluxNames)
 Regress1DTimeSeriesMapWithTTestContours(MovingAvg.Window12Months.(TempFluxNames{i}), MovingAvgTimeSeries.Window12Months.Net.HemDif,...
     TempFluxNames{i}, 'Net NH-SH');
Regress1DTimeSeriesMapWithTTestContours(MovingAvg.Window12Months.(TempFluxNames{i}), MovingAvgTimeSeries.Window12Months.Net.GlobalAvg,...
    TempFluxNames{i}, 'Global Mean Net Flux');
Regress1DTimeSeriesMapWithTTestContours(MovingAvg.Window12Months.(TempFluxNames{i}), MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to90,...
    TempFluxNames{i}, 'Global Precip Asymmetry');
Regress1DTimeSeriesMapWithTTestContours(MovingAvg.Window12Months.(TempFluxNames{i}), MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14,...
    TempFluxNames{i}, '0to14deg Precip Asymmetry');


end

PlotLatRSquared(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']),MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net,FluxNames,MonthFilterSize,'Net')
PlotLatRSquared(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']),MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Precip,FluxNames,MonthFilterSize,'Precip')

[AX,H1,H2] = plotyy(1:length(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14'),MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14',...
    1:length(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14'),MovingAvgTimeSeries.Window12Months.Net.HemDif0to14')
legend(['Precip Asym 0-14deg'],['Net Hemdif 0-14deg'])
grid on;
set(gca,'xtick',12-2-(MonthFilterSize-1):12:time)
set(AX(2),'XTickLabel',[])
set(gca,'XTickLabel',2000:2012)
xlabel('Year End');
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['NetPrecipAsym0to14','.png']);


[NewCorr,NewRSq]= PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to20,MovingAvgTimeSeries.Window12Months.Net.HemDif20to90,'PrecipAsym0to20','NetHemDif20to90',12)
PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to20,MovingAvgTimeSeries.Window12Months.Net.HemDif0to20,'PrecipAsym0to20','NetHemDif0to20',12)

%negative correlation of -0.34 when precip leads Net by 6 months
% 1:11 vs 2:12 (1 month lag)
% 1:10 vs 3:12 (2 months lag)
%Timelags just decay for longwave..


for i=1:length(TempFluxNames)

MakeLizMap;colormap(lizmap);set(gca,'FontSize',20)
load geoid; ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true); geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5]);
geoshow(var(MovingAvg.Window12Months.(TempFluxNames{i}),0,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
caxis([0 max(max(var(MovingAvg.Window12Months.(TempFluxNames{i}))))])
title(['Variance of ', TempFluxNames{i}])
grid on;
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[TempFluxNames{i},'12MnthMAVariance.png']);
hold off;
end




%end 2014-05-05


%%%%%%%%%

% %%%%
% NetClearCloudNames = {'Net','NetClear','NetCloud'};
% for i = 1:3
%    NetClearCloud.(NetClearCloudNames{i}) = MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(NetClearCloudNames{i});
% end
% 
% PlotLatRSquared(NetClearCloud,IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']),fieldnames(NetClearCloud),MonthFilterSize,[''])

    
for MonthFilterSize=[3,6]
for i=1:3
    PlotLatRegressCoefs(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']),MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}),FluxNames,MonthFilterSize,FluxNameToCompareWith{i})
end
end

%strcat('net',fieldnames(MovingAvgTimeSeries.Window3Months.Net))
% 
% corr(IndicesMvgAvg.Window3Months.(IndexNames{1})',   MovingAvgTimeSeries.Window3Months.(FluxNames{1}).NH')
% corr(IndicesMvgAvg.Window3Months.(IndexNames{2})',   MovingAvgTimeSeries.Window3Months.(FluxNames{1}).SH')
% corr(IndicesMvgAvg.Window3Months.(IndexNames{3})',   MovingAvgTimeSeries.Window3Months.(FluxNames{1}).HemDif')
% corr(IndicesMvgAvg.Window3Months.(IndexNames{3})',   MovingAvgTimeSeries.Window3Months.(FluxNames{1}).HemSum_EquatorTo20')
% corr(IndicesMvgAvg.Window3Months.(IndexNames{3})',   MovingAvgTimeSeries.Window3Months.(FluxNames{1}).SH_EquatorTo20')
% corr(IndicesMvgAvg.Window3Months.(IndexNames{3})',   MovingAvgTimeSeries.Window3Months.(FluxNames{1}).NH_EquatorTo20')


%%%

% 
% MonthFilterSize = 6; %or can be 6, 12
% IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']) = IndexMovingAverage(MonthFilterSize);
% IndexMovingAverage(12) %even for 12 months this works!!
% SampleWindow = fieldnames(IndicesMvgAvg);
% IndexNames = fieldnames(IndicesMvgAvg.(SampleWindow{1}));
% close all;
% for i=1:3
% MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}) = MovingAverageFilterAllTiles(RawFlux.(FluxNames{i}),FluxNames{i},MonthFilterSize);
% close all;
% end


% [MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.NH,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.SH,...
%     MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemDif,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemSum] ...
%     = FluxInLatitudinalBand(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).Net,0,90,'Net',MonthFilterSize);
% 
% [MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.NH_EquatorTo20,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.SH_EquatorTo20,...
%     MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemDif_EquatorTo20,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemSum_EquatorTo20] ...
%     = FluxInLatitudinalBand(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).Net,0,20,'Net',MonthFilterSize);

% IndicesMvgAvg.Window3Months.(IndexNames{1})

% corr(IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{1})',   MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{1}).NH')
% corr(IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{2})',   MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{1}).SH')
% corr(IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{3})',   MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{1}).HemDif')
% corr(IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{3})',   MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{1}).HemSum_EquatorTo20')
% corr(IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{3})',   MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{1}).SH_EquatorTo20')
% corr(IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']).(IndexNames{3})',   MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{1}).NH_EquatorTo20')


% MovingAvg6Month.Net = MovingAverageFilterAllTiles(RawFlux.Net,'Net',MonthFilterSize);
% MovingAvg6Month.SW = MovingAverageFilterAllTiles(RawFlux.SW,'SW',MonthFilterSize);
% MovingAvg6Month.LW = MovingAverageFilterAllTiles(RawFlux.Net,'LW',MonthFilterSize);
% Indices_6MonthMA = IndexMovingAverage(MonthFilterSize);
% 

%%%%
MonthFilterSize = 12; %or can be 6, 12
IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']) = IndexMovingAverage(MonthFilterSize);
SampleWindow = fieldnames(IndicesMvgAvg);
IndexNames = fieldnames(IndicesMvgAvg.(SampleWindow{1}));
close all;
for i=1:3
MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i}) = MovingAverageFilterAllTiles(RawFlux.(FluxNames{i}),FluxNames{i},MonthFilterSize);
close all;
end

close all;
for i=1:3
   for j=1:3
    [MA_12Month_Beta.(IndexNames{j}).(FluxNames{i}), blah, blah, blah,MA_12Month_Alpha.(IndexNames{j}).(FluxNames{i})] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg.Window12Months.(FluxNames{i}) ,...
    IndicesMvgAvg.Window12Months.(IndexNames{j}),[FluxNames{i}, ' ', num2str(MonthFilterSize),'-Month Mvg Avg'],[IndexNames{j}, ' ', num2str(MonthFilterSize),'-Month Mvg Avg']);
    close all;
    RegressedFlux = MWATimeLatContourRegression(MA_12Month_Alpha.(IndexNames{j}).(FluxNames{i}),MA_12Month_Beta.(IndexNames{j}).(FluxNames{i}),MovingAvg.Window12Months.(FluxNames{i}) ,IndicesMvgAvg.Window12Months.(IndexNames{i}),...
        [FluxNames{i}, ' ', num2str(MonthFilterSize),'-Month Mvg Avg'],[IndexNames{j}, ' ', num2str(MonthFilterSize),'-Month Mvg Avg'],MonthFilterSize);
    close all;
   end
end

% mean(RegressedFlux,3) %this really is 0 everywhere.
% contourf(mean(RegressedFlux,3));colorbar
% mean(mean(RegressedFlux,3),2)



% [MA_3Month_Beta_Net.NAM, blah, blah, blah,MA_3Month_Alpha_Net.NAM] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.Net,Indices_3MonthMA.NAM,'Net 3-Month MA','NAM 3-Month MA');
% [MA_3Month_Beta_Net.NAM, blah, blah, blah,MA_3Month_Alpha_Net.NAM] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.Net,Indices_3MonthMA.NAM,['Net ',num2str(MonthFilterSize),'-Month MA'],...
%     ['NAM', num2str(MonthFilterSize),'-Month MA'])

% [MA_3Month_Beta_Net.SAM, blah, blah, blah,MA_3Month_Alpha_Net.SAM] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.Net,Indices_3MonthMA.SAM,'Net 3-Month MA','SAM 3-Month MA');
% [MA_3Month_Beta_Net.NINO34, blah, blah, blah,MA_3Month_Alpha_Net.NINO34] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.Net,Indices_3MonthMA.NINO34,'Net 3-Month MA','NINO34 3-Month MA');
% 
% close all;
% MWATimeLatContourRegression(MA_3Month_Alpha_Net.NAM,MA_3Month_Beta_Net.NAM,MovingAvg3Month,Indices_3MonthMA.NAM,'Net 3-Month MA','NAM 3-Month MA',MonthFilterSize);
% close all;
% MWATimeLatContourRegression(MA_3Month_Alpha_Net.SAM,MA_3Month_Beta_Net.SAM,MovingAvg3Month,IndicesMWA.SAM,'Net 3-Month MA','SAM 3-Month MA',MonthFilterSize);
% close all;
% MWATimeLatContourRegression(MA_3Month_Alpha_Net.NINO34,MA_3Month_Beta_Net.NINO34,MovingAvg3Month,IndicesMWA.NINO34,'Net 3-Month MA','NINO34 3-Month MA',MonthFilterSize);
% 
% %%%
% [MA_3Month_Beta_LW.NAM, blah, blah, blah,MA_3Month_Alpha_LW.NAM] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.LW,Indices_3MonthMA.NAM,'LW 3-Month MA','NAM 3-Month MA')
% [MA_3Month_Beta_LW.SAM, blah, blah, blah,MA_3Month_Alpha_LW.SAM] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.LW,Indices_3MonthMA.SAM,'LW 3-Month MA','SAM 3-Month MA')
% [MA_3Month_Beta_LW.NINO34, blah, blah, blah,MA_3Month_Alpha_LW.NINO34] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.LW,Indices_3MonthMA.NINO34,'LW 3-Month MA','NINO34 3-Month MA')
% hold off;
% MWATimeLatContourRegression(MA_3Month_Alpha_LW.NAM,MA_3Month_Beta_LW.NAM,MovingAvg3Month,IndicesMWA.NAM,'LW 3-Month MA','NAM 3-Month MA',MonthFilterSize);
% MWATimeLatContourRegression(MA_3Month_Alpha_LW.SAM,MA_3Month_Beta_LW.SAM,MovingAvg3Month,IndicesMWA.SAM,'LW 3-Month MA','SAM 3-Month MA',MonthFilterSize);
% close all;
% MWATimeLatContourRegression(MA_3Month_Alpha_LW.NINO34,MA_3Month_Beta_LW.NINO34,MovingAvg3Month,IndicesMWA.NINO34,'LW 3-Month MA','NINO34 3-Month MA',MonthFilterSize);
% 
% %%%
% [MA_3Month_Beta_SW.NAM, blah, blah, blah,MA_3Month_Alpha_SW.NAM] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.SW,Indices_3MonthMA.NAM,'SW 3-Month MA','NAM 3-Month MA')
% close all;
% [MA_3Month_Beta_SW.SAM, blah, blah, blah,MA_3Month_Alpha_SW.SAM] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.SW,Indices_3MonthMA.SAM,'SW 3-Month MA','SAM 3-Month MA')
% close all;
% [MA_3Month_Beta_SW.NINO34, blah, blah, blah,MA_3Month_Alpha_SW.NINO34] = Regress1DTimeSeriesMapWithTTestContours(MovingAvg3Month.SW,Indices_3MonthMA.NINO34,'SW 3-Month MA','NINO34 3-Month MA')
% hold off;
% MWATimeLatContourRegression(MA_3Month_Alpha_SW.NAM,MA_3Month_Beta_SW.NAM,MovingAvg3Month,IndicesMWA.NAM,'SW 3-Month MA','NAM 3-Month MA',MonthFilterSize);
% close all;
% MWATimeLatContourRegression(MA_3Month_Alpha_SW.SAM,MA_3Month_Beta_SW.SAM,MovingAvg3Month,IndicesMWA.SAM,'SW 3-Month MA','SAM 3-Month MA',MonthFilterSize);
% close all;
% MWATimeLatContourRegression(MA_3Month_Alpha_SW.NINO34,MA_3Month_Beta_SW.NINO34,MovingAvg3Month,IndicesMWA.NINO34,'SW 3-Month MA','NINO34 3-Month MA',MonthFilterSize);

%%%%DO THIS BELOW FOR 3 MONTHS and 6 MONTHS to get the TEXT FILES!!
% [MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.NH,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.SH,...
%     MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemDif,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemSum] ...
%     = FluxInLatitudinalBand(MovingAvg3Month.Net,0,90,'Net',MonthFilterSize);
% 
% [MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.NH_EquatorTo15,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.SH_EquatorTo15,...
%     MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemDif_EquatorTo15,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemSum_EquatorTo15] ...
%     = FluxInLatitudinalBand(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).Net,0,15,'Net',MonthFilterSize);


% [MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.NH_EquatorTo20,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.SH_EquatorTo20,...
%     MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemDif_EquatorTo20,MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).Net.HemSum_EquatorTo20] ...
%     = FluxInLatitudinalBand(MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).Net,0,20,'Net',MonthFilterSize);


% SubsetFluxIndices = find(strcmp('Net',FluxNames) | strcmp('Temp',FluxNames) | strcmp('Precip',FluxNames));
% 
% SubsetFluxIndices = find(strcmp('Net',FluxNames) | strcmp('SW',FluxNames) | strcmp('LW',FluxNames));
% 
% SubsetFluxIndices = find(strcmp('LW',FluxNames) | strcmp('LWclear',FluxNames) | strcmp('LWCF',FluxNames));
% SubsetFluxIndices = find(strcmp('SW',FluxNames) | strcmp('SWclear',FluxNames) | strcmp('SWCF',FluxNames));
% 
% SubsetFluxIndices = find(strcmp('NetCloud',FluxNames) | strcmp('SWCF',FluxNames) | strcmp('LWCF',FluxNames));
% SubsetFluxIndices = find(strcmp('NetClear',FluxNames) | strcmp('SWclear',FluxNames) | strcmp('LWclear',FluxNames));
% 
% SubsetFluxIndices = find(strcmp('Net',FluxNames) | strcmp('SW',FluxNames) | strcmp('LW',FluxNames)| strcmp('LWCF',FluxNames));
% SubsetFluxIndices = find(strcmp('Net',FluxNames) | strcmp('SWCF',FluxNames) | strcmp('Precip',FluxNames)| strcmp('LWCF',FluxNames));
% SubsetFluxIndices = find(strcmp('Net',FluxNames) | strcmp('SWCF',FluxNames) | strcmp('LWCF',FluxNames)| strcmp('SWclear',FluxNames)| strcmp('LWclear',FluxNames));
% SubsetFluxIndices = find(strcmp('Net',FluxNames) | strcmp('SWCF',FluxNames) | strcmp('LWCF',FluxNames)| strcmp('SWclear',FluxNames)| strcmp('LWclear',FluxNames));
% 
% SubsetFluxNames = FluxNames(SubsetFluxIndices);
% for h = 1:size(SubsetFluxNames)
%     i = SubsetFluxIndices(h);
%     SubsetFlux.(FluxNames{i}) = MovingAvg.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{i});
% end
% PlotLatRSquared(SubsetFlux,IndicesMvgAvg.(['Window',num2str(MonthFilterSize),'Months']),SubsetFluxNames,MonthFilterSize,[SubsetFluxNames{:}])
% clearvars SubsetFlux SubsetFluxNames SubsetFluxIndices
% 

% plot(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14'-mean(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14'))
% hold on;
% plot(MovingAvgTimeSeries.Window12Months.Net.HemDif0to14'-mean(MovingAvgTimeSeries.Window12Months.Net.HemDif0to14'),'r')
% 
% corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14',MovingAvgTimeSeries.Window12Months.Net.HemDif0to14')
% % corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to90',MovingAvgTimeSeries.Window12Months.Net.HemDif0to90')
% 
% corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14',MovingAvgTimeSeries.Window12Months.Net.HemDif0to14')
% corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14',MovingAvgTimeSeries.Window12Months.Net.HemDif0to90')
% 
% ExtratropicalComposite = MovingAvgTimeSeries.Window12Months.Net.HemDif14to30+MovingAvgTimeSeries.Window12Months.Net.HemDif30to49+MovingAvgTimeSeries.Window12Months.Net.HemDif49to90;
% ExtratropicalComposite = MovingAvgTimeSeries.Window12Months.Net.HemDif20to90;
% 
% corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14',ExtratropicalComposite')
% corr(MovingAvgTimeSeries.Window12Months.Net.HemDif0to90',ExtratropicalComposite')

% Flux1 = MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to20;
% Flux2 = MovingAvgTimeSeries.Window12Months.Net.HemDif20to90;

% MaxLag = 10;
% TimeLeadCorr = zeros(1,MaxLag);
% TimeLagCorr = zeros(1,MaxLag);
% 
% for i=1:MaxLag
%     TimeLeadCorr(i) = corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14(1:end-i)',MovingAvgTimeSeries.Window12Months.Net.HemDif0to14(i+1:end)');
%     TimeLagCorr(i)=corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14(i+1:end)',MovingAvgTimeSeries.Window12Months.Net.HemDif0to14(1:end-i)');
% end
% TimeLeadCorr = fliplr(TimeLeadCorr);
% tempCorr = corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14',MovingAvgTimeSeries.Window12Months.Net.HemDif0to14');
% 
% TimeLaggedCorr = [TimeLeadCorr tempCorr TimeLagCorr];
% 
% TimeLaggedCorr(MaxLag+1)
% MonthsBeforeAfter = [-MaxLag:MaxLag]
% plot(MonthsBeforeAfter,TimeLaggedCorr)
% MonthOfMaxCorr =  find(max(abs(TimeLaggedCorr))==abs(TimeLaggedCorr));
% TimeLaggedCorr(MonthOfMaxCorr)

% a = PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14,MovingAvgTimeSeries.Window12Months.Net.HemDif0to14,'PrecipAsym0to14','NetHemDif0to14',12)

% 
% for i=1:13
%     i
%     'precip leads Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14(1:end-i)',ExtratropicalComposite(i+1:end)')
%     'precip lags Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14(i+1:end)',ExtratropicalComposite(1:end-i)') %Crrln of 0.50 when precip lags extra-tropical net
%     %by 9-10 months
%     corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14(i+1:end)',MovingAvgTimeSeries.Window12Months.Net.HemDif0to90(1:end-i)')
% %crrln decreases when hemdif from 0 to 90 is included
% end
% 
% for i=1:13
%     i
%     'precip leads Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Net.HemDif0to14(1:end-i)',ExtratropicalComposite(i+1:end)')
%     'precip lags Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Net.HemDif0to14(i+1:end)',ExtratropicalComposite(1:end-i)') %Crrln of 0.50 when precip lags extra-tropical net
% %time lag between tropical/extratropical maxe out at 5-6 months both ways..
% end
% 
% %%%%%
% for i=1:10
%     i
%     'precip leads Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to20(1:end-i)',MovingAvgTimeSeries.Window12Months.Net.HemDif0to20(i+1:end)')
%     'precip lags Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to20(i+1:end)',MovingAvgTimeSeries.Window12Months.Net.HemDif0to20(1:end-i)')
% end
% %negative correlation of -0.34 when precip leads Net by 6 months
% % 1:11 vs 2:12 (1 month lag)
% % 1:10 vs 3:12 (2 months lag)
% %Timelags just decay for longwave..
% 
% for i=1:13
%     i
%     'precip leads Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to20(1:end-i)',ExtratropicalComposite(i+1:end)')
%     'precip lags Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to20(i+1:end)',ExtratropicalComposite(1:end-i)') %Crrln of 0.50 when precip lags extra-tropical net
%     %by 9-10 months
%     corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to20(i+1:end)',MovingAvgTimeSeries.Window12Months.Net.HemDif0to90(1:end-i)')
% %crrln decreases when hemdif from 0 to 90 is included
% end
% 
% for i=1:13
%     i
%     'precip leads Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Net.HemDif0to20(1:end-i)',ExtratropicalComposite(i+1:end)')
%     'precip lags Net by i months'
%     corr(MovingAvgTimeSeries.Window12Months.Net.HemDif0to20(i+1:end)',ExtratropicalComposite(1:end-i)') %Crrln of 0.50 when precip lags extra-tropical net
% %time lag between tropical/extratropical maxe out at 5-6 months both ways..
% end
% 
% 
% %%%%%%%
% 
% i=6;
% TimePts =1:length(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14(1:end-i));
% plotyy(TimePts,MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14(1:end-i)',TimePts,MovingAvgTimeSeries.Window12Months.Net.HemDif0to14(i+1:end)')
% i=10;
% TimePts =1:length(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14(1:end-i));
% plotyy(TimePts,MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14(1:end-i)',TimePts,ExtratropicalComposite(i+1:end)')
% 
% corr(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14',MovingAvgTimeSeries.Window12Months.Precip.HemDif0to14')
% 
% TimePts = 1:length(MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14);
% 
% plotyy(TimePts',MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14',TimePts',MovingAvgTimeSeries.Window12Months.Net.HemDif0to14')
% plotyy(TimePts',MovingAvgTimeSeries.Window12Months.Precip.AsymIndex0to14',TimePts',MovingAvgTimeSeries.Window12Months.Net.HemDif0to90')
