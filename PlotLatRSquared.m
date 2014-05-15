function blah = PlotLatRSquared(Flux3D,Indices,FluxNames,MonthFilterSize,TimeSeriesName)
load LatWeights
cmap = jet;
Style = repmat({'--','-'},1,length(FluxNames));
Markers = repmat({'none','none'},1,length(FluxNames));

IndexIterator = fieldnames(Indices);

lats = size(Flux3D,1);
longs = size(Flux3D,2);
time=length(Indices.(IndexIterator{1}));
Factor = 180/lats;
lats = size(Flux3D,1);
Factor = 180/lats;

TotalRegressSum = zeros(length(FluxNames),1);
for j=1:size(IndexIterator)
        Flux2Legend = FluxNames;
        for i=1:size(FluxNames)
        Flux1 = Flux3D.(FluxNames{i});
        LatFlux = squeeze(mean(Flux1,2)); %180x156
        TimeSeries = Indices.(IndexIterator{j});
        TimeSeries3D = ones(1,1,time);
        TimeSeries3D(1,:) = TimeSeries;
        TimeSeries3D = repmat(TimeSeries3D,[180 1]);
        
        [LatBeta,pValues, LatAlpha] = RegressLatvsTimeTTestContours(Flux3D.(FluxNames{i}),Indices.(IndexIterator{j}),FluxNames{i},IndexIterator{j});
        Significant = pValues > 0.975 | pValues < 0.025;
    % size(LatAlpha)
    % size(LatBeta)
    % size(TimeSeries)
    % % size(LatBeta.*TimeSeries)
    % size(repmat(LatAlpha,[1,length(TimeSeries)]) )
        
        RegressedFlux = repmat(LatAlpha,[1,length(TimeSeries)]) + bsxfun(@times,repmat(LatBeta,[1,length(TimeSeries)]),TimeSeries);
        
        TimeMean =  repmat(mean(LatFlux,2),[1,length(TimeSeries)]) ;
        SSTotal = sum((LatFlux-TimeMean ).^2,2);
        SSRes = sum((LatFlux-RegressedFlux).^2,2);%%why is this so high? Why is SS_Res >> SS_Total?
        NumMonthsAfterMAFiltering = time-(MonthFilterSize-1);

         R_Squared = 1 - SSRes./SSTotal;
         
        figure(1)
            h=plot(sind(LatWeights(:,1)),R_Squared,sind(LatWeights(Significant,1)),R_Squared(Significant));
set(h(1),'color',cmap(i*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{i},'Marker',Markers{i},'LineWidth',3)
            if numel(h) > 1
            set(h(2),'Marker','d','LineStyle','None','MarkerSize',13,'MarkerFaceColor',cmap(i*floor(64/ceil(length(FluxNames))),:),'MarkerEdgeColor','k');
            set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

            end
            %when nothing is significant, then h(2) is completely empty...
            %hold on;
            %plot(sind(LatWeights(Significant,1)),TempFlux(Significant),'d','MarkerSize',13,'MarkerFaceColor',cmap(j*floor(64/ceil(length(FluxNames))),:),'color',cmap(j*floor(64/ceil(length(FluxNames))),:))
        grid on;
        hold on;
        set(gca,'FontSize',18)
        set(gca,'GridLineStyle','--')
        set(gca,'xtick',sind(-90:10:90))
%         set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{j}),2),1))/180:180.5)
        set(gca,'xticklabel',num2cell(-90:10:90))
            xlabel('latitude')

        Flux2Legend{i} = [Flux2Legend{i},' = ',num2str(round(TotalRegressSum(i)))];
        
        figure(2)
            g=plot(sind(LatWeights(:,1)),LatBeta,sind(LatWeights(Significant,1)),LatBeta(Significant));
set(g(1),'color',cmap(i*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{i},'Marker',Markers{i},'LineWidth',3)
            if numel(g) > 1
            set(g(2),'Marker','d','LineStyle','None','MarkerSize',13,'MarkerFaceColor',cmap(i*floor(64/ceil(length(FluxNames))),:),'MarkerEdgeColor','k');
            set(get(get(g(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

            end
            %when nothing is significant, then h(2) is completely empty...
            %hold on;
            %plot(sind(LatWeights(Significant,1)),TempFlux(Significant),'d','MarkerSize',13,'MarkerFaceColor',cmap(j*floor(64/ceil(length(FluxNames))),:),'color',cmap(j*floor(64/ceil(length(FluxNames))),:))
        grid on;
        hold on;
        set(gca,'FontSize',18)
        set(gca,'GridLineStyle','--')
        set(gca,'xtick',sind(-90:10:90))
%         set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{j}),2),1))/180:180.5)
        set(gca,'xticklabel',num2cell(-90:10:90))
            xlabel('latitude')

%         Flux2Legend{i} = [Flux2Legend{i},' = ',num2str(round(TotalRegressSum(i)))];
        end
    figure(1)
    legend(Flux2Legend,'Location','EastOutside')

    title([num2str(MonthFilterSize),'-Month Mvg Avg. R-Squared. of each flux over ',TimeSeriesName, IndexIterator{j}, ...
        ' (W/m). 95% T-Test Sig Marked. Weighted sum of all lat. regressions in legend'])
        set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',[num2str(MonthFilterSize),'Month_MA_LatRSquared',TimeSeriesName,'_',IndexIterator{j},'.png']);
    hold off;
    
    figure(2)
    legend(Flux2Legend,'Location','EastOutside')
    title([num2str(MonthFilterSize),'-Month Mvg Avg. Regression. of each flux over ', TimeSeriesName,IndexIterator{j}, ...
        ' (W/m). 95% T-Test Sig Marked. Weighted sum of all lat. regressions in legend']);
        set(gcf,'paperposition',[0 0 20 10])

    print(gcf,'-dpng','-r300',[num2str(MonthFilterSize),'Month_MA_LatRegress',TimeSeriesName,'_',IndexIterator{j},'.png']);
    hold off;

end


end