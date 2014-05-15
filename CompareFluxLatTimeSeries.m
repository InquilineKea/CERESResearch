function blah = CompareFluxLatTimeSeries(MovingAvgTimeSeries,SubsetFluxNames,SubsetFluxIndices,FluxNames,MonthFilterSize,LowerLat,HigherLat,Hemisphere)
time = length(MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{1}).(['Global', num2str(LowerLat),'to',num2str(HigherLat)]));

cmap= jet;Style = repmat({'--','-'},1,length(SubsetFluxNames));
Markers = repmat({'none','none'},1,length(SubsetFluxNames));
for i = 1:size(SubsetFluxNames)
       k = SubsetFluxIndices(i);
h=        plot(MovingAvgTimeSeries.(['Window',num2str(MonthFilterSize),'Months']).(FluxNames{k}).([Hemisphere, num2str(LowerLat),'to',num2str(HigherLat)]));
set(gca,'xtick',12-2-(MonthFilterSize-1):12:time)
set(gca,'XTickLabel',2000:2012)
xlabel('Year End');
        set(h(1),'color',cmap(i*floor(64/ceil(length(SubsetFluxNames))),:),'LineStyle',Style{i},'Marker',Markers{i},'LineWidth',3)
            if numel(h) > 1
            set(h(2),'Marker','d','LineStyle','None','MarkerSize',13,'MarkerFaceColor',cmap(i*floor(64/ceil(length(SubsetFluxNames))),:),'MarkerEdgeColor','k');
            set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

            end
            %when nothing is significant, then h(2) is completely empty...
            %hold on;
            %plot(sind(LatWeights(Significant,1)),TempFlux(Significant),'d','MarkerSize',13,'MarkerFaceColor',cmap(j*floor(64/ceil(length(FluxNames))),:),'color',cmap(j*floor(64/ceil(length(FluxNames))),:))
        grid on;
        hold on;
        set(gca,'FontSize',18)
        set(gca,'GridLineStyle','--')
hold on;
end
legend(SubsetFluxNames)


title([num2str(LowerLat),'-',num2str(HigherLat),'deg ', num2str(MonthFilterSize),' Month Moving Avg ', Hemisphere])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[SubsetFluxNames{:},num2str(MonthFilterSize),'MonthMA',Hemisphere,num2str(LowerLat),'-',num2str(HigherLat),'.png']);
hold off;

end