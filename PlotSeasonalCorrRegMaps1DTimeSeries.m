function mapss = PlotSeasonalCorrRegMaps1DTimeSeries(Flux1,TimeSeries)

NumYears = size(TimeSeries)/12;

%SPRING
SpringIndices = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(1:3,1,13);
%WINTER
WinterIndices = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(10:12,1,13);
%AUTUMN
AutumnIndices = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(7:9,1,13);
%SUMMER
SummerIndices = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(4:6,1,13);

subplot(2,2,1)
Corr1DTimeSeriesMap(Flux1(:,:,WinterIndices),TimeSeries(WinterIndices));
subplot(2,2,2)
Corr1DTimeSeriesMap(Flux1(:,:,SpringIndices),TimeSeries(SpringIndices));
subplot(2,2,3)
Corr1DTimeSeriesMap(Flux1(:,:,SummerIndices),TimeSeries(SummerIndices));
subplot(2,2,4)
Corr1DTimeSeriesMap(Flux1(:,:,AutumnIndices),TimeSeries(AutumnIndices));
title(subplot(2,2,1),' Winter Crrln')
title(subplot(2,2,2),' Spring Crrln')
title(subplot(2,2,3),' Summer Crrln')
title(subplot(2,2,4),' Autumn Crrln')
set(gcf,'paperposition',[0 0 20 10])
   [ax,h3]=suplabel(['AllSeasonTimeSeriesCorr-',num2str(inputname(1)),'-',num2str(inputname(2))] ,'t'); 
   set(h3,'FontSize',30) 
print(gcf,'-dpng','-r300',['AllSeasonTimeSeriesCorr_',num2str(inputname(1)),'-',num2str(inputname(2)),'.png']);
hold off;

subplot(2,2,1)
Regress1DTimeSeriesMap(Flux1(:,:,WinterIndices),TimeSeries(WinterIndices));
subplot(2,2,2)
Regress1DTimeSeriesMap(Flux1(:,:,SpringIndices),TimeSeries(SpringIndices));
subplot(2,2,3)
Regress1DTimeSeriesMap(Flux1(:,:,SummerIndices),TimeSeries(SummerIndices));
subplot(2,2,4)
Regress1DTimeSeriesMap(Flux1(:,:,AutumnIndices),TimeSeries(AutumnIndices));
title(subplot(2,2,1),' Winter Regre')
title(subplot(2,2,2),' Spring Regre')
title(subplot(2,2,3),' Summer Regre')
title(subplot(2,2,4),' Autumn Regre')
set(gcf,'paperposition',[0 0 20 10])
   [ax,h3]=suplabel(['AllSeasonTimeSeriesRegres-',num2str(inputname(1)),'-',num2str(inputname(2))] ,'t'); 
   set(h3,'FontSize',30) 
print(gcf,'-dpng','-r300',['AllSeasonTimeSeriesRegress_',num2str(inputname(1)),'-',num2str(inputname(2)),'.png']);


end