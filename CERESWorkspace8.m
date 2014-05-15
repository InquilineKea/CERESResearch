load MonthlyFluxDeparturesStruct.mat
FluxNames = fieldnames(MonthlyFluxDepartures);
load LatitudinalDifFluxTimeSeriesAsField.mat
load Indices.mat
Indices = rmfield(Indices,'PNA');
Indices = rmfield(Indices,'SAO');
Indices = rmfield(Indices,'T');
IndexNames = fieldnames(Indices);
MonthlyFluxDepartures = rmfield(MonthlyFluxDepartures,'Clouds');
MonthlyFluxDepartures = rmfield(MonthlyFluxDepartures,'Precip');
FluxNames = fieldnames(MonthlyFluxDepartures);


SnowCoverIndexAllYears = csvread('SnowCoverNH.csv');
IceCoverAllYears = csvread('IceCoverExtent.csv');

%dividing this into seasons might be easier to do here..

SnowCover= SnowCoverIndexAllYears(find(SnowCoverIndexAllYears(:,1) == 2000 & SnowCoverIndexAllYears(:,2)==3):find(SnowCoverIndexAllYears(:,1) == 2013& SnowCoverIndexAllYears(:,2)==2),:);
IceCover= IceCoverAllYears(find(IceCoverAllYears(:,1) == 2000 & IceCoverAllYears(:,2)==3):find(IceCoverAllYears(:,1) == 2013& IceCoverAllYears(:,2)==2),:);
SnowCover = SnowCover(:,3);
IceCover = IceCover(:,3);
SnowCover = SubtractClimatologyFromTimeSeries(SnowCover);
IceCover = SubtractClimatologyFromTimeSeries(IceCover);

Indices.Ice = IceCover';
Indices.Snow = SnowCover';

%WINTER
SeasonIndices.Winter = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(10:12,1,13);
%SPRING
SeasonIndices.Spring = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(1:3,1,13);
%SUMMER
SeasonIndices.Summer = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(4:6,1,13);
%AUTUMN
SeasonIndices.Autumn = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(7:9,1,13);
SeasonNames = fieldnames(SeasonIndices);

for k =1:length(SeasonNames)
   corr(SnowCover(SeasonIndices.(SeasonNames{k})),IceCover(SeasonIndices.(SeasonNames{k})))
end

corr(SnowCover,IceCover)

%shouldn't everything be normalized to month first? i mean.. yes.. the
%seasonal averaging normalizes things a bit but not completely..

load RegressIndexFluxSeasonal.mat

for i=1:8
   [a,a,a,MaxBeta(i)] = Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.Snow,FluxNames{i},'NH Snow Cover')
   [a,a,a,MaxBeta(i)] = Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.Ice,FluxNames{i},'NH Ice Cover')
end

for i =1:9
   [RegressIndexFlux.Snow.(FluxNames{i}), LatRegressFluxDif.Snow.(FluxNames{i}),FluxSumAllLats.Snow.(FluxNames{i})]= Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.Snow,FluxNames{i},'NH Snow Cover')
   [RegressIndexFlux.Ice.(FluxNames{i}), LatRegressFluxDif.Ice.(FluxNames{i}),FluxSumAllLats.Ice.(FluxNames{i})] = Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.Ice,FluxNames{i},'NH Ice Cover')  
end


for j=1:length(FluxNames)
    for k=1:length(SeasonNames)
TempFlux = MonthlyFluxDepartures.(FluxNames{j});
% size(TempFlux(SeasonIndices.(SeasonNames{k}),:,:))
    [RegressIndexFlux.Ice.(FluxNames{j}).(SeasonNames{k}),LatRegressFluxDif.Ice.(FluxNames{j}).(SeasonNames{k}),FluxSumAllLats.Ice.(FluxNames{j}).(SeasonNames{k})] = Regress1DTimeSeriesMapWithTTestContours(TempFlux(:,:,SeasonIndices.(SeasonNames{k})),IceCover(SeasonIndices.(SeasonNames{k}),:),FluxNames{j},[SeasonNames{k},' Ice Extent']);
    [RegressIndexFlux.Snow.(FluxNames{j}).(SeasonNames{k}),LatRegressFluxDif.Snow.(FluxNames{j}).(SeasonNames{k}),FluxSumAllLats.Snow.(FluxNames{j}).(SeasonNames{k})]=Regress1DTimeSeriesMapWithTTestContours(TempFlux(:,:,SeasonIndices.(SeasonNames{k})),SnowCover(SeasonIndices.(SeasonNames{k}),:),FluxNames{j},[SeasonNames{k},' Snow Extent']);
    end
end
IndexNames = fieldnames(RegressIndexFlux);
RegressIndexFlux.Ice= orderfields(RegressIndexFlux.Ice);
RegressIndexFlux.Snow= orderfields(RegressIndexFlux.Snow);

save('RegressIndexFluxSeasonal.mat','RegressIndexFlux')


for i=5:6
    for k =1:length(SeasonNames)
        for j=1:length(FluxNames)
            lats = size(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k}),1);
            Factor = 180/lats;
            TempMonthFlux = MonthlyFluxDepartures.(FluxNames{j});
            TempIndex = Indices.(IndexNames{i});
            [LatBeta,pValues] = RegressLatvsTimeTTestContours(TempMonthFlux(:,:,SeasonIndices.(SeasonNames{k})),TempIndex(SeasonIndices.(SeasonNames{k})),FluxNames{j},IndexNames{i});
            Significant = pValues > 0.975 | pValues < 0.025;
    %IS THIS A SUM OR A MEAN?? IF SUM, IT WOULD BE W/M. IF MEAN, IT WOULD BE
    %W/M^2, RIGHT??
    %         FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2).*LatWeights(:,2);
            FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k}),Factor),2);
            TempFlux = FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k});
            if strcmp(FluxNames{j},'Net')
                %h = plot(FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{j},'Marker',Markers{j},'LineWidth',9);
                netHandle=plot(sind(LatWeights(:,1)),TempFlux,sind(LatWeights(Significant,1)),TempFlux(Significant));
                set(netHandle(1),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{j},'Marker',Markers{j},'LineWidth',3)
                if numel(netHandle) > 1
                set(netHandle(2),'Marker','d','LineStyle','None','MarkerSize',13,'MarkerFaceColor',cmap(j*floor(64/ceil(length(FluxNames))),:),'MarkerEdgeColor','k');
                set(get(get(netHandle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

                end
                %             hold on;
                %h = plot(sind(LatWeights(:,1)),TempFlux(Significant),'d','MarkerSize',13,'MarkerFaceColor',cmap(j*floor(64/ceil(length(FluxNames))),:),'color',cmap(j*floor(64/ceil(length(FluxNames))),:));

            else
                %plot(FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{j},'Marker',Markers{j},'LineWidth',3)    
    %             h= plot(sind(LatWeights(:,1)),FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{j},'Marker',Markers{j},'LineWidth',3,...
    %                 sind(LatWeights(Significant,1)),TempFlux(Significant),'d','MarkerSize',13,'MarkerFaceColor',cmap(j*floor(64/ceil(length(FluxNames))),:),'color',cmap(j*floor(64/ceil(length(FluxNames))),:))
                h=plot(sind(LatWeights(:,1)),TempFlux,sind(LatWeights(Significant,1)),TempFlux(Significant));
                set(h(1),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{j},'Marker',Markers{j},'LineWidth',3)
                if numel(h) > 1
                set(h(2),'Marker','d','LineStyle','None','MarkerSize',13,'MarkerFaceColor',cmap(j*floor(64/ceil(length(FluxNames))),:),'MarkerEdgeColor','k');
                set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

                end
                %when nothing is significant, then h(2) is completely empty...
                %hold on;
                %plot(sind(LatWeights(Significant,1)),TempFlux(Significant),'d','MarkerSize',13,'MarkerFaceColor',cmap(j*floor(64/ceil(length(FluxNames))),:),'color',cmap(j*floor(64/ceil(length(FluxNames))),:))
            end
            grid on;
            hold on;
            set(gca,'FontSize',18)
            set(gca,'GridLineStyle','--')
            set(gca,'xtick',sind(-90:10:90))
    %         set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{j}),2),1))/180:180.5)
            set(gca,'xticklabel',num2cell(-90:10:90))
        end
        legend(FluxNames,'Location','EastOutside')
        uistack(netHandle, 'top') %MUST come AFTER legend
        title(['Total Latitudinal Regre. of ', IndexNames{i}, ' in ',SeasonNames{k},' (W/m). 95% T-Test Sig Marked with diamonds'])
        xlabel('latitude')
        set(gcf,'paperposition',[0 0 20 10])
        print(gcf,'-dpng','-r300',['TotalLatRegress',IndexNames{i},SeasonNames{k},'.png']);
        hold off;
    end
end

%%this is just seasonal though. what about yearly? NH-SH flux dif?

fid = fopen('AllNHMinusSHFluxes.txt','w');
for i=1:6
%     IndexNames{i}
%         SeasonNames{k}
        for j=1:length(FluxNames)
    for k =1:length(SeasonNames)
%             FluxNames{j}
        TempFlux = FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k});
        AvgNHMinusSHFlux.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k}) = sum(TempFlux(end/2+1:end))-sum(TempFlux(1:end/2));
        fprintf(fid, '%s\n',[IndexNames{i}, FluxNames{j}, SeasonNames{k}, ' NH - SH = ', num2str(AvgNHMinusSHFlux.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k}))]);
        end
    end
end
fclose(fid);

AvgNHMinusSHFlux.Ice.Net.Winter
