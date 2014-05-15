function blah = PlotLatRegressCoefs(Flux3D,Indices,FluxNames,MonthFilterSize,TimeSeriesName)
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

        for i=1:3%for each of the fluxes: Net, LW, SW

        [LatBeta,pValues, LatAlpha] = RegressLatvsTimeTTestContours(Flux3D.(FluxNames{i}),Indices.(IndexIterator{j}),FluxNames{i},IndexIterator{j});
        Significant = pValues > 0.975 | pValues < 0.025;
        %IS THIS A SUM OR A MEAN?? IF SUM, IT WOULD BE W/M. IF MEAN, IT WOULD BE
        %W/M^2, RIGHT??
        %         FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2).*LatWeights(:,2);
                %FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2);
                %beta=Regress1DTimeSeriesMapWithTTestContours(Flux3D.(FluxNames{i}),Indices.(IndexIterator{j}),FluxNames{i},IndexIterator{j});
        Flux1 = Flux3D.(FluxNames{i});
        TimeSeries = Indices.(IndexIterator{j});
        %TimeSeries = NetDifMonthTotal;
        %Flux1 = NetClimatologySubtracted;
        TimeSeries3D = ones(1,1,time);
        TimeSeries3D(1,1,:) = TimeSeries;
        % size(TimeSeries3D)
        TimeSeries3D = repmat(TimeSeries3D,[180 360 1]);
        % size(TimeSeries3D)
        %Beta = bsxfun(@rdivide,time*(sum(Flux1.*TimeSeries3D,3))-sum(Flux1,3).*sum(TimeSeries3D,3),time*sum(TimeSeries3D.*TimeSeries3D,3)-sum(TimeSeries3D,3).*sum(TimeSeries3D,3))*std(TimeSeries);
        %multiply by std(TimeSeries) in end so that you come out with units of

%         Sxy = sum(Flux1.*TimeSeries3D,3);
%         Sxx = sum(TimeSeries3D.*TimeSeries3D,3);
%         Syy = sum(Flux1.*Flux1,3);
%         SSE = Syy-Sxy.*Sxy./Sxx;
%         Sy = sum(Flux1,3);
%         Sx = sum(TimeSeries);
%         Beta = bsxfun(@rdivide,time*(Sxy)-Sx*Sy,time*Sxx-Sx.*Sx);
% 
%         FluxSumsAllLats.(IndexIterator{j}).(FluxNames{i}) = sum(Beta,2); %will this result in the correct units on the left axis??
%         TempFlux = FluxSumsAllLats.(IndexIterator{j}).(FluxNames{i});
        TotalRegressSum(i) = sum(LatBeta.*LatWeights(:,2));

            h=plot(sind(LatWeights(:,1)),LatBeta,sind(LatWeights(Significant,1)),LatBeta(Significant));
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
        Flux2Legend{i} = [Flux2Legend{i},' = ',num2str(round(TotalRegressSum(i)))];
        end
        legend(Flux2Legend,'Location','EastOutside')
    title([num2str(MonthFilterSize),'-Month Mvg Avg. Regre. of each index over ',TimeSeriesName,' ', IndexIterator{j}, ' (W/m). 95% T-Test Sig Marked. Weighted sum of all lat. regressions in legend'])
    xlabel('latitude')
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',[num2str(MonthFilterSize),'Month_MA_LatRegress',TimeSeriesName,'_',IndexIterator{j},'.png']);
    hold off;

end


end