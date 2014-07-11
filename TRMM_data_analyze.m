

cd TRMM

%1440 longs, 400 lats (from -50 to 50)
lats = ncread('3B43.20120101.7.nc','latitude');
longs = ncread('3B43.20120101.7.nc','longitude');

% lats(length(lats)/2)
% lats(length(lats)+1/2)

% TRMMPrecip = cat(3,TRMMPrecip1,TRMMPrecip2);

dirls = dir('*.nc');

RawTRMMPrecip =zeros(length(lats),length(longs),length(dirls));

for filenum=1:length(dirls)
    RawTRMMPrecip(:,:,filenum) = ncread(dirls(filenum).name,'pcp')';
end

time = length(dirls);

MonthFilterSize = 12;
cd ..

% [TRMMPrecip.NH0to20,TRMMPrecip.SH0to20, TRMMPrecip.HemDif0to20,TRMMPrecip.HemSum0to20,TRMMPrecip.AsymIndex0to20] = GeneralizedFluxInLatitudinalBand(RawTRMMPrecip,0,20,'TRMMPrecip',12,lats);
[TRMMPrecip.NH0to15,TRMMPrecip.SH0to15, TRMMPrecip.HemDif0to15,TRMMPrecip.HemSum0to15,TRMMPrecip.AsymIndex0to15] = GeneralizedFluxInLatitudinalBand(RawTRMMPrecip,0,15,'TRMMPrecip',12,lats);
% [TRMMPrecip.NH0to30,TRMMPrecip.SH0to30, TRMMPrecip.HemDif0to30,TRMMPrecip.HemSum0to30,TRMMPrecip.AsymIndex0to30] = GeneralizedFluxInLatitudinalBand(RawTRMMPrecip,0,30,'TRMMPrecip',12,lats);
% 
% %MW this first, and then take weighted average sums?
% MovingAverage = zeros(length(lats),length(longs),time);
% days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
% allMonthWeights = repmat(days_per_month,1,12);
% allMonthWeights = [31 30 31 30 31 31 30 31 30 31 allMonthWeights 31 28];
% allMonthWeights(48) = 29; %leap year, http://www.wolframalpha.com/input/?i=months+between+march+2000+and+february+2004
% allMonthWeights(96) = 29;
% allMonthWeights3D = ones(1,1,time);
% allMonthWeights3D(1,1,:) = allMonthWeights;
% allMonthWeights3D = repmat(allMonthWeights3D,[length(lats) length(longs) 1]);
% % allMonthWeights3D = permute(allMonthWeights3D,[3 1 2]);
% moving_sum = @(n, x) filter(ones(1,n), 1, x);
% MovingAverage= bsxfun(@rdivide,moving_sum(MonthFilterSize, TRMMPrecip .* allMonthWeights3D) ...
%     ,moving_sum(MonthFilterSize,allMonthWeights3D)); %or Mean Lat Flux here
% % MovingAverage = permute(MovingAverage,[2 3 1]);
% %MWA means 1st 11 values are junk
% clear allMonthWeights3D;
% MovingAverage = MovingAverage(:,:,MonthFilterSize:end);
% NumMonthsAfterMAFiltering = time-(MonthFilterSize-1);
% TRMM_MA_Precip = zeros(length(lats),length(longs),NumMonthsAfterMAFiltering);
% %Flux = LWclear;
% for i=1:NumMonthsAfterMAFiltering
%     Indices = mod(i,12):12:NumMonthsAfterMAFiltering;
%     if mod(i,12)==0
%         Indices = 12:12:NumMonthsAfterMAFiltering;
%     end
%     TRMM_MA_Precip(:,:,i) = MovingAverage(:,:,i) - mean(MovingAverage(:,:,Indices),3);
% end
% clear MovingAverage;
% 
% 
% %%%%%%
% 
% TRMMPrecipLatMeans = squeeze(mean(TRMM_MA_Precip,2));
% LatWeightsTRMM = cosd(lats)./sum(cosd(lats));
% TRMMPrecipModified = zeros(length(lats),NumMonthsAfterMAFiltering);
% for i =1:NumMonthsAfterMAFiltering
% TRMMPrecipModified(:,i) = TRMMPrecipLatMeans(:,i).*LatWeightsTRMM;
% end
% 
% TRMMPrecip.SH0to20 = sum(TRMMPrecipModified(lats > -20 & lats < 0,:));
% TRMMPrecip.NH0to20 = sum(TRMMPrecipModified(lats > 0 & lats < 20,:));
% TRMMPrecip.Global0to20 = sum(TRMMPrecipModified(lats > -20 & lats < 20,:))/2;
% TRMMPrecip.HemDif0to20 = sum(TRMMPrecipModified(lats > 0 & lats < 20,:)) - sum(TRMMPrecipModified(lats > -20 & lats < 0,:));
% % TRMMPrecip.Asym0to20 = 2*TRMMPrecip.HemDif0to20./TRMMPrecip.Global0to20;
% 
% TRMMPrecip.Global0to15 = sum(TRMMPrecipModified(lats > -15 & lats < 15,:))/2;
% TRMMPrecip.HemDif0to15 = sum(TRMMPrecipModified(lats > 0 & lats < 15,:)) - sum(TRMMPrecipModified(lats > -15 & lats < 0,:));
% % TRMMPrecip.Asym0to15 = 2*TRMMPrecip.HemDif0to15./TRMMPrecip.Global0to15;
% 
% TRMMPrecip.Global0to30 = sum(TRMMPrecipModified(lats > -30 & lats < 30,:))/2;
% TRMMPrecip.HemDif0to30 = sum(TRMMPrecipModified(lats > 0 & lats < 30,:)) - sum(TRMMPrecipModified(lats > -30 & lats < 0,:));
% % TRMMPrecip.Asym0to30 = 2*TRMMPrecip.HemDif0to30./TRMMPrecip.Global0to30;

cd ..
save('TRMMPrecip.mat','TRMMPrecip')

load('TRMMPrecip.mat')
MonthFilterSize = 12;
% load('2014-04-30.mat','MovingAvgTimeSeries')
AllNetFields = fieldnames(MovingAvgTimeSeries.Window12Months.Net);
AllPrecipFields = fieldnames(MovingAvgTimeSeries.Window12Months.Precip);
AllTRMMPrecipFields = fieldnames(TRMMPrecip);

% clear fid;clear fid2;
FluxFile = ['NetVsTRMMPrecipMaxTimeLaggedCorr', num2str(MonthFilterSize), '_MonthMA','.txt'];
fid = fopen(FluxFile,'w');
% TimeLaggedCorrFile = ['AllTimeLaggedCorr_NetVsTRMMPrecip', num2str(MonthFilterSize), '_MonthMA','.txt'];
TimeLaggedCorrCSV =['AllTimeLaggedCorr_NetVsTRMMPrecip', num2str(MonthFilterSize), '_MonthMA','.csv'];
TimeLaggedpValCSV =['AllTimeLaggedpVal_NetVsTRMMPrecip', num2str(MonthFilterSize), '_MonthMA','.csv'];
fid2 = fopen(TimeLaggedCorrCSV,'w');fclose(fid2);
fid3 = fopen(TimeLaggedpValCSV,'w');fclose(fid3);
%         dlmwrite(TimeLaggedCorrCSV,[0,-25:25],'-append');
for i=1:length(AllNetFields)
    for j=1:length(AllTRMMPrecipFields)
        [MaxCorr,RealMonthOfMaxCorr,AllTimeLaggedCorr,AllTimeLaggedPVal] = PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i}),TRMMPrecip.(AllTRMMPrecipFields{j}),['Net',AllNetFields{i}],['TRMMPrecip',AllTRMMPrecipFields{j}],12,1);
        fprintf(fid, '%s\n',['Net',AllNetFields{i},' vs ', 'TRMMPrecip',AllTRMMPrecipFields{j}, ' Corr = ', num2str(MaxCorr), ' at ',AllTRMMPrecipFields{j},' Leading Net by ', num2str(RealMonthOfMaxCorr), ' Months Time Lag']);
%         fprintf(fid2, '%s\n',['Net',AllNetFields{i},' vs ', 'TRMMPrecip',AllTRMMPrecipFields{j}, ' Corr = ', num2str(AllTimeLaggedCorr)]);
fid2 = fopen(TimeLaggedCorrCSV,'a');
fprintf(fid2,['Net',AllNetFields{i},' vs. TRMMPrecip',AllTRMMPrecipFields{j},',']);
fclose(fid2);
        dlmwrite(TimeLaggedCorrCSV,AllTimeLaggedCorr,'-append');
fid3 = fopen('AllTimeLaggedPVal.csv','a');
fprintf(fid3,['Net',AllNetFields{i},' vs. TRMMPrecip',AllTRMMPrecipFields{j},',']);
fclose(fid3);
        dlmwrite('AllTimeLaggedPVal.csv',AllTimeLaggedPVal,'-append');
%         corr(MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})',MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})')
    end
end
fprintf(fid, '**************************\n');
fclose(fid);

% fclose(fid2);

MonthFilterSize=12;
csvID = fopen(['AllTimeLaggedCorr_NetVsTRMMPrecip', num2str(MonthFilterSize), '_MonthMA','.csv'],'r');
C = textscan(csvID,'%s %*[^\n]','delimiter', ',')
CompareNames=cellstr(C{1});
% index = regexp(CompareNames,'PrecipAsymIndex0to[1234](.*)');
index = regexp(CompareNames,'PrecipAsymIndex0to[12](.*)');
% index2 = strfind(CompareNames,'PrecipAsym');
% index3 = regexp(CompareNames,'NetHemDif[1234](.*)to[349]');
index3 = regexp(CompareNames,'SH[1234](.*)to90');
Index = find(not(cellfun('isempty', index)) &  not(cellfun('isempty', index3)) );
M = csvread(['AllTimeLaggedCorr_NetVsTRMMPrecip', num2str(MonthFilterSize), '_MonthMA','.csv'],0,1)
cmap = jet;
hold on;
% LineStyles = {'-','--',':'};
LineStyles = {'-','--'};
for i = 1:length(Index)
plot([-25:25],M(Index(i),:),'linewidth',3,'linestyle',LineStyles{mod(i,2)+1},'color',cmap(i*floor(64/ceil(length(Index))),:))
end
grid on;
set(gca,'FontSize',20)
xlabel('Time Lag, with TRMM Precip leading Net (Months)')
xlabel('<-- Net Leading TRMM Precip (Months)            ||           Precip Leading Net (Months) -->')
ylabel('Correlation between Net and TRMM Precip')
legend(CompareNames(Index),'FontSize',8,'Location','NorthEast')
	set(gca,'GridLineStyle','--')
	set(gcf,'paperposition',[0 0 20 10])
	print(gcf,'-dpng','-r300',['NetVsTRMMPrecipTimeLaggedCorr.png']);
fclose(csvID);

% clear fid;clear fid2;
FluxFile = ['NetVsPrecipMaxTimeLaggedCorr', num2str(MonthFilterSize), '_MonthMA','.txt'];
fid = fopen(FluxFile,'w');
% TimeLaggedCorrFile = ['AllTimeLaggedCorr_NetVsPrecip', num2str(MonthFilterSize), '_MonthMA','.txt'];
TimeLaggedCorrCSV =['AllTimeLaggedCorr_NetVsPrecip', num2str(MonthFilterSize), '_MonthMA','.csv'];
TimeLaggedpValCSV =['AllTimeLagged_pVal_NetVsPrecip', num2str(MonthFilterSize), '_MonthMA','.csv'];

fid2 = fopen(TimeLaggedCorrCSV,'w');fclose(fid2);
fid3 = fopen(TimeLaggedpValCSV,'w');fclose(fid3);
%         dlmwrite(TimeLaggedCorrCSV,[0,-25:25],'-append');
for i=1:length(AllNetFields)
    for j=1:length(AllPrecipFields)
        [MaxCorr,RealMonthOfMaxCorr,AllTimeLaggedCorr,TimeLaggedpVal] = PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i}),MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{j}),'Net',AllPrecipFields{j},12,1);
        fprintf(fid, '%s\n',['Net',AllNetFields{i},' vs ', 'Precip',AllPrecipFields{j}, ' Corr = ', num2str(MaxCorr), ' at ',AllPrecipFields{j},' Leading Net by ', num2str(RealMonthOfMaxCorr), ' Months Time Lag']);
%         fprintf(fid2, '%s\n',['Net',AllNetFields{i},' vs ', 'Precip',AllPrecipFields{j}, ' Corr = ', num2str(AllTimeLaggedCorr)]);
fid2 = fopen(TimeLaggedCorrCSV,'a');
fprintf(fid2,['Net',AllNetFields{i},' vs. Precip',AllPrecipFields{j},',']);
fclose(fid2);
dlmwrite(TimeLaggedCorrCSV,AllTimeLaggedCorr,'-append');

fid3 = fopen(TimeLaggedpValCSV,'a');
fprintf(fid3,['Net',AllNetFields{i},' vs. Precip',AllPrecipFields{j},',']);
fclose(fid3);
dlmwrite(TimeLaggedpValCSV,AllTimeLaggedCorr,'-append');

%         corr(MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})',MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})')
    end
end
fprintf(fid, '**************************\n');
fclose(fid);

MonthFilterSize=12;
csvID = fopen('AllTimeLaggedCorr_NetVsPrecip12_MonthMA.csv','r');
C = textscan(csvID,'%s %*[^\n]','delimiter', ',');
CompareNames=cellstr(C{1});
% index = regexp(CompareNames,'PrecipAsymIndex0to[1234](.*)');
% index = regexp(CompareNames,'PrecipAsymIndex0to[12](.*)');
index = regexp(CompareNames,'PrecipAsymIndex0to14');
% index2 = strfind(CompareNames,'PrecipAsym');
% index3 = regexp(CompareNames,'NetHemDif[1234](.*)to[349]');
 index3 = regexp(CompareNames,'SH[1234](.*)to[349]');
Index = find(not(cellfun('isempty', index)) &  not(cellfun('isempty', index3)) );
M = csvread('AllTimeLaggedCorr_NetVsPrecip12_MonthMA.csv',0,1)
cmap = jet;
hold on;
% LineStyles = {'-','--',':'};
LineStyles = {'-','--'};
for i = 1:length(Index)
plot([-25:25],M(Index(i),:),'linewidth',3,'linestyle',LineStyles{mod(i,2)+1},'color',cmap(i*floor(64/ceil(length(Index))),:))
end
grid on;
set(gca,'FontSize',20)
xlabel('<-- Net Leading GPCP Precip (Months)            ||           Precip Leading Net (Months) -->')
ylabel('Correlation between Net and GPCP Precip')
legend(CompareNames(Index),'FontSize',8,'Location','NorthEast')
	set(gca,'GridLineStyle','--')
	set(gcf,'paperposition',[0 0 20 10])
	print(gcf,'-dpng','-r300',['NetVsGPCPPrecipTimeLaggedCorr.png']);
fclose(csvID);

%%%%%%%%%%
%2014-07-08

[xc,lags] = xcorr(TRMMPrecip.AsymIndex0to15,TRMMPrecip.AsymIndex0to15)

[xc,lags] = xcorr(TRMMPrecip.AsymIndex0to15,MovingAvgTimeSeries.Window12Months.Net.NH15to90)
[xc,lags] = xcorr(detrend(TRMMPrecip.AsymIndex0to15),detrend(MovingAvgTimeSeries.Window12Months.Net.NH15to90))
[xc,lags] = xcorr(detrend(TRMMPrecip.AsymIndex0to15),detrend(MovingAvgTimeSeries.Window12Months.Net.SH15to90))
%[xc,lags] = xcorr(detrend(TRMMPrecip.AsymIndex0to15),detrend(MovingAvgTimeSeries.Window12Months.Net.HemDif15to90))
% [xc,lags] = xcorr(TRMMPrecip.AsymIndex0to15,MovingAvgTimeSeries.Window12Months.Net.NH15to90)
[xc,lags] = xcorr(MovingAvgTimeSeries.Window12Months.Net.HemDif15to90,TRMMPrecip.AsymIndex0to15)
X = MovingAvgTimeSeries.Window12Months.Net.HemDif15to90;Y = TRMMPrecip.AsymIndex0to15;
X=detrend(MovingAvgTimeSeries.Window12Months.Net.HemDif15to90);Y = detrend(TRMMPrecip.AsymIndex0to15);
 X = detrend(MovingAvgTimeSeries.Window12Months.Net.HemDif15to90(1:end-5));Y=detrend(MovingAvgTimeSeries.Window12Months.Net.HemDif15to90(6:end))
%  X = MovingAvgTimeSeries.Window12Months.Net.HemDif15to90(1:end-5);Y=MovingAvgTimeSeries.Window12Months.Net.HemDif15to90(6:end);
[xc,lags] = xcorr(X,Y,50,'coeff');[m,i] = max(abs(xc));tau = lags(i);plot(lags,xc);grid on;tau
xlabel('<-- Net Leading TRMM Precip (Months)            ||           Precip Leading Net (Months) -->')
X = detrend(MovingAvgTimeSeries.Window12Months.Net.NH15to90);Y=detrend(MovingAvgTimeSeries.Window12Months.Net.NH15to90);
X=detrend(MovingAvgTimeSeries.Window12Months.Net.SH15to90);Y = detrend(TRMMPrecip.AsymIndex0to15);
X=detrend(IndicesMWA.NINO34);Y = detrend(TRMMPrecip.AsymIndex0to15);
X=detrend(IndicesMWA.NINO34);Y = detrend(MovingAvgTimeSeries.Window12Months.Net.SH15to90);


%%%%%%%%

csvID = fopen(['AllTimeLaggedCorr_NetVsTRMMPrecip', num2str(MonthFilterSize), '_MonthMA','.csv'],'r');
C = textscan(csvID,'%s %*[^\n]','delimiter', ',')
CompareNames=cellstr(C{1});
% index = regexp(CompareNames,'PrecipAsymIndex0to[1234](.*)');
indexTRMMPrecip1 = regexp(CompareNames,'PrecipAsymIndex0to15');

csvID2 = fopen('AllTimeLaggedCorr_NetVsPrecip12_MonthMA.csv','r');
C2 = textscan(csvID2,'%s %*[^\n]','delimiter', ',');
CompareNames2=cellstr(C2{1});
indexGPCPPrecip1 = regexp(CompareNames2,'PrecipAsymIndex0to14');
% index2 = strfind(CompareNames,'PrecipAsym');
indexTRMMPrecip2 = regexp(CompareNames,'NetHemDif15to90');
indexGPCPPrecip2 = regexp(CompareNames2,'NetHemDif15to90');
IndexTRMM = find(not(cellfun('isempty', indexTRMMPrecip2)) &  not(cellfun('isempty', indexTRMMPrecip1)) );
IndexGPCP = find(not(cellfun('isempty', indexGPCPPrecip2)) &  not(cellfun('isempty', indexGPCPPrecip1)) );
M = csvread(['AllTimeLaggedCorr_NetVsTRMMPrecip', num2str(MonthFilterSize), '_MonthMA','.csv'],0,1)
M2 = csvread('AllTimeLaggedCorr_NetVsPrecip12_MonthMA.csv',0,1)
cmap = jet;
hold on;
% LineStyles = {'-','--',':'};
LineStyles = {'-','--'};
for i = 1:length(IndexTRMM)
plot([-25:25],M(IndexTRMM(i),:),'linewidth',3,'linestyle',LineStyles{mod(i,2)+1},'color',cmap(i*floor(64/ceil(length(Index))),:))
end
hold on;
for j=1:length(IndexGPCP)
plot([-25:25],M2(IndexGPCP(j),:),'linewidth',3,'linestyle',LineStyles{mod(i+j,2)+1},'color',cmap((i+j)*floor(64/ceil(length(Index))),:))    
end
grid on;
set(gca,'FontSize',20)
xlabel('<-- Net Leading TRMM Precip (Months)            ||           Precip Leading Net (Months) -->')
ylabel('Correlation between Net and GPCP Precip')
legend({'TRMM Precip Asym (0-15deg) vs extratropical Net (15-90deg) NH-SH Correlation','GPCP Precip Asym (0-15deg) vs extratropical Net (15-90deg) NH-SH Correlation'},'FontSize',8,'Location','NorthEast')
	set(gca,'GridLineStyle','--')
	set(gcf,'paperposition',[0 0 20 10])
	print(gcf,'-dpng','-r300',['TropicalPrecipVsExtratropicalNetTimeLag.png']);
fclose(csvID);


%%%%%%%%%%


FluxFile = ['PrecipVsTRMMPrecipMaxTimeLaggedCorr', num2str(MonthFilterSize), '_MonthMA','.txt'];
fid = fopen(FluxFile,'w');
for i=1:length(AllPrecipFields)
    for j=1:length(AllTRMMPrecipFields)
        [MaxCorr,RealMonthOfMaxCorr,AllTimeLaggedCorr] = PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{i}),TRMMPrecip.(AllTRMMPrecipFields{j}),'Precip','AllTRMMPrecipFields{j}',12,0);
        fprintf(fid, '%s\n',['Precip',AllPrecipFields{i},' vs TRMMPrecip',AllTRMMPrecipFields{j}, ' Corr = ', num2str(MaxCorr), ' at ',AllTRMMPrecipFields{j},' Leading Precip by ', num2str(RealMonthOfMaxCorr), ' Months Time Lag']);
%         corr(MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})',MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})')
    end
end
fprintf(fid, '**************************\n');
fclose(fid);



