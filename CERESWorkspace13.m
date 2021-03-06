AllNetFields = fieldnames(MovingAvgTimeSeries.Window12Months.Net);
AllPrecipFields = fieldnames(MovingAvgTimeSeries.Window12Months.Precip);

FluxFile = ['PrecipVsNetCorr', num2str(MonthFilterSize), '_MonthMA_','.txt'];
fid = fopen(FluxFile,'w');
for i=1:length(AllPrecipFields)
    for j=1:length(AllNetFields)
        fprintf(fid, '%s\n',['Precip',AllPrecipFields{i},'vsNet',AllNetFields{j}, ' Corr = ', num2str( ...
            corr(MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{i})',MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{j})') )]);
%         corr(MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{i})',MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})')
        
    end
end
fprintf(fid, '**************************\n');
fclose(fid);

FluxFile = ['PrecipVsNetMaxTimeLaggedCorr', num2str(MonthFilterSize), '_MonthMA','.txt'];
fid = fopen(FluxFile,'w');
for i=1:length(AllPrecipFields)
    for j=1:length(AllNetFields)
        [MaxCorr,RealMonthOfMaxCorr] = PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{i}),MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{j}),'Precip','Net',12,0);
        
        fprintf(fid, '%s\n',['Precip',AllPrecipFields{i},' vs Net',AllNetFields{j}, ' Corr = ', num2str(MaxCorr), ' at Net Leading Precip by ', num2str(RealMonthOfMaxCorr), ' Months Time Lag']);
%         corr(MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{i})',MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})')
    end
end
fprintf(fid, '**************************\n');
fclose(fid);

AllIndexFields = IndexNames;
FluxFile = ['PrecipVsIndicesMaxTimeLaggedCorr', num2str(MonthFilterSize), '_MonthMA','.txt'];
fid = fopen(FluxFile,'w');
for i=1:length(AllPrecipFields)
    for j=1:length(AllIndexFields)
        [MaxCorr,RealMonthOfMaxCorr] = PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{i}),IndicesMvgAvg.Window12Months.(AllIndexFields{j}),'Precip','AllIndexFields{j}',12,0);
        fprintf(fid, '%s\n',['Precip',AllPrecipFields{i},' vs ',AllIndexFields{j}, ' Corr = ', num2str(MaxCorr), ' at ',AllIndexFields{j},' Leading Precip by ', num2str(RealMonthOfMaxCorr), ' Months Time Lag']);
%         corr(MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{i})',MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})')
    end
end
fprintf(fid, '**************************\n');
fclose(fid);

FluxFile = ['NetVsIndicesMaxTimeLaggedCorr', num2str(MonthFilterSize), '_MonthMA','.txt'];
fid = fopen(FluxFile,'w');
for i=1:length(AllNetFields)
    for j=1:length(AllIndexFields)
        [MaxCorr,RealMonthOfMaxCorr] = PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i}),IndicesMvgAvg.Window12Months.(AllIndexFields{j}),'Net','AllIndexFields{j}',12,0);
        fprintf(fid, '%s\n',['Net',AllNetFields{i},' vs ',AllIndexFields{j}, ' Corr = ', num2str(MaxCorr), ' at ',AllIndexFields{j},' Leading Net by ', num2str(RealMonthOfMaxCorr), ' Months Time Lag']);
%         corr(MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})',MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})')
    end
end
fprintf(fid, '**************************\n');
fclose(fid);



for TimeLag=8:10
FluxFile = ['PrecipVsNetCorrsAt',num2str(TimeLag),'Months', num2str(MonthFilterSize), '_MonthMA','.txt'];
fid = fopen(FluxFile,'w');
for i=1:length(AllPrecipFields)
    for j=1:length(AllNetFields)
        [MaxCorr,RealMonthOfMaxCorr,TimeSeriesOfCorrs] = PlotFluxSidebySide(MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{i}),MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{j}),'Precip','Net',12,0);
        CorrAtTimeLag = TimeSeriesOfCorrs((length(TimeSeriesOfCorrs)-1)/2+TimeLag);
        
        fprintf(fid, '%s\n',['Precip',AllPrecipFields{i},' vs Net',AllNetFields{j}, ' Corr = ', num2str(CorrAtTimeLag), ' at Net Leading Precip by ', num2str(TimeLag), ' Months Time Lag']);
%         corr(MovingAvgTimeSeries.Window12Months.Precip.(AllPrecipFields{i})',MovingAvgTimeSeries.Window12Months.Net.(AllNetFields{i})')
        
    end
end
fprintf(fid, '**************************\n');
fclose(fid);
end