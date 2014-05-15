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

NormIndices.NAM = Indices.NAM/std(Indices.NAM);
NormIndices.SAM= Indices.SAM/std(Indices.SAM);
NormIndices.NAO = Indices.NAO/std(Indices.NAO);
NormIndices.NINO34= Indices.NINO34/std(Indices.NINO34);

[a,b] = RegressLatvsTimeTTestContours(MonthlyFluxDepartures.Net,Indices.NAM);

%%%%%

load('RegressIndexFlux.mat')
load LatWeights.mat

%test = structfun(@(x)(orderfields(x)),RegressIndexFlux);
RegressIndexFlux.NAM = orderfields(RegressIndexFlux.NAM);
RegressIndexFlux.SAM = orderfields(RegressIndexFlux.SAM);
RegressIndexFlux.NAO= orderfields(RegressIndexFlux.NAO);
RegressIndexFlux.NINO34= orderfields(RegressIndexFlux.NINO34);
FluxNames = fieldnames(RegressIndexFlux.NINO34);
MonthlyFluxDepartures = orderfields(MonthlyFluxDepartures);

Style = repmat({'--','-'},1,length(FluxNames));
Markers = repmat({'none','none'},1,length(FluxNames));
i=1;j=1;

%WINTER
SeasonIndices.Winter = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(10:12,1,13);
%SPRING
SeasonIndices.Spring = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(1:3,1,13);
%SUMMER
SeasonIndices.Summer = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(4:6,1,13);
%AUTUMN
SeasonIndices.Autumn = 12*ceil(find(repmat(1:3,1,13))/3)-12+repmat(7:9,1,13);
SeasonNames = fieldnames(SeasonIndices);


for i=1:length(IndexNames) %regress all fields [and indices] over total NH net flux and total SH net flux. might do it for total NH+SH flux too.
    for j=1:length(FluxNames)
            TempFlux = MonthlyFluxDepartures.(FluxNames{j});
            TempIndices = Indices.(IndexNames{i});
        for k=1:length(SeasonNames)
    [RegressIndexFlux.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k}), LatRegressFluxDif.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k}),FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}).(SeasonNames{k})] = Regress1DTimeSeriesMapWithTTestContours(TempFlux(:,:,SeasonIndices.(SeasonNames{k})),TempIndex(SeasonIndices.(SeasonNames{k})),FluxNames{j},[IndexNames{i}, ' during ', SeasonNames{k}]);
    %NH Total, SH Total, NH-SH Total, NH+SH Total

    %WorldRegress2 = Regress1DTimeSeriesMap(MonthlyFluxDepartures.(FluxNames{i}), FluxTimeSeries.Net.AllLats.N,MonthlyFluxDepartures.(FluxNames{i}),FluxNames{i},[Flux1DName{1}, LatName{5}, HemisphereName{1},'H']);
    % clear WorldRegress1;
    % clear WorldRegress2
        end
    end
end
%this entire loop took 30 minutes

%Struct field assignment overwrites a value with class "double" ??

RegressIndexFlux.NAM = orderfields(RegressIndexFlux.NAM);
RegressIndexFlux.SAM = orderfields(RegressIndexFlux.SAM);
RegressIndexFlux.NAO= orderfields(RegressIndexFlux.NAO);
RegressIndexFlux.NINO34= orderfields(RegressIndexFlux.NINO34);

save('RegressIndexFluxSeasonal.mat','RegressIndexFlux')

for i=1:4
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
