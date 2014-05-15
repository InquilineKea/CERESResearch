function mapss = PlotSeasonalCorrRegMaps1DTimeSeries(Flux1,TimeSeries)

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
   [ax,h3]=suplabel(['AllSeasonTimeSeriesCorr_',num2str(inputname(1)),'-',num2str(inputname(2))] ,'t'); 
   set(h3,'FontSize',30) 

set(gcf,'paperposition',[0 0 20 10])
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
   [ax,h3]=suplabel(['AllSeasonTimeSeriesRegres_',num2str(inputname(1)),'-',num2str(inputname(2))] ,'t'); 
   set(h3,'FontSize',30) 

set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['AllSeasonTimeSeriesRegress_',num2str(inputname(1)),'-',num2str(inputname(2)),'.png']);


end