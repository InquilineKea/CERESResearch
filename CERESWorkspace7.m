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

for i =1:9
   [RegressIndexFlux.Snow.(FluxNames{i}), LatRegressFluxDif.Snow.(FluxNames{i}),FluxSumAllLats.Snow.(FluxNames{i})]= Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.Snow,FluxNames{i},'NH Snow Cover')
   [RegressIndexFlux.Ice.(FluxNames{i}), LatRegressFluxDif.Ice.(FluxNames{i}),FluxSumAllLats.Ice.(FluxNames{i})] = Regress1DTimeSeriesMapWithTTestContours(MonthlyFluxDepartures.(FluxNames{i}),Indices.Ice,FluxNames{i},'NH Ice Cover')  
end


%test = structfun(@(x)(orderfields(x)),RegressIndexFlux);
RegressIndexFlux.NAM = orderfields(RegressIndexFlux.NAM);
RegressIndexFlux.SAM = orderfields(RegressIndexFlux.SAM);
RegressIndexFlux.NAO= orderfields(RegressIndexFlux.NAO);
RegressIndexFlux.NINO34= orderfields(RegressIndexFlux.NINO34);
FluxNames = fieldnames(RegressIndexFlux.NINO34);
MonthlyFluxDepartures = orderfields(MonthlyFluxDepartures);

IndexNames = fieldnames(RegressIndexFlux);

RegressIndexFlux.Snow = orderfields(RegressIndexFlux.Snow,RegressIndexFlux.SAM);
RegressIndexFlux.Ice = orderfields(RegressIndexFlux.Ice,RegressIndexFlux.SAM);

save('RegressIndexFlux.mat','RegressIndexFlux')

Style = repmat({'--','-'},1,length(FluxNames));
Markers = repmat({'none','none'},1,length(FluxNames));
i=1;j=1;
%%for everything at once
for i=1:6
    for j=1:length(FluxNames)
        lats = size(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),1);
        Factor = 180/lats;
        [LatBeta,pValues] = RegressLatvsTimeTTestContours(MonthlyFluxDepartures.(FluxNames{j}),Indices.(IndexNames{i}),FluxNames{j},IndexNames{i});
        Significant = pValues > 0.975 | pValues < 0.025;
%IS THIS A SUM OR A MEAN?? IF SUM, IT WOULD BE W/M. IF MEAN, IT WOULD BE
%W/M^2, RIGHT??
%         FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2).*LatWeights(:,2);
        FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2);
        TempFlux = FluxSumsAllLats.(IndexNames{i}).(FluxNames{j});
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
    title(['Total Latitudinal Regre. of ', IndexNames{i}, ' (W/m). 95% T-Test Sig Marked with diamonds'])
    xlabel('latitude')
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['TotalLatRegress',IndexNames{i},'.png']);
    hold off;
end

%%%try more limited set

FluxNamesAll = fieldnames(RegressIndexFlux.NINO34);
FluxNames = FluxNames(1:4);
FluxNames = FluxNamesAll(5:end);
FluxNames = FluxNamesAll;

for i=1:6
    for j=5:length(FluxNames)
        lats = size(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),1);
        Factor = 180/lats;
        [LatBeta,pValues] = RegressLatvsTimeTTestContours(MonthlyFluxDepartures.(FluxNames{j}),Indices.(IndexNames{i}),FluxNames{j},IndexNames{i});
        Significant = pValues > 0.975 | pValues < 0.025;
%IS THIS A SUM OR A MEAN?? IF SUM, IT WOULD BE W/M. IF MEAN, IT WOULD BE
%W/M^2, RIGHT??
%         FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2).*LatWeights(:,2);
        FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2);
        TempFlux = FluxSumsAllLats.(IndexNames{i}).(FluxNames{j});
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
    title(['Total Latitudinal Regre. of ', IndexNames{i}, ' (W/m). 95% T-Test Sig Marked with diamonds'])
    xlabel('latitude')
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['TotalLatRegress',IndexNames{i},'.png']);
    hold off;
end

%%%%%%%%%%%%%%%CUSTOMIZED VERSION

for i=1:6
    for j=5:length(FluxNames)
        lats = size(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),1);
        Factor = 180/lats;
        [LatBeta,pValues] = RegressLatvsTimeTTestContours(MonthlyFluxDepartures.(FluxNames{j}),Indices.(IndexNames{i}),FluxNames{j},IndexNames{i});
        Significant = pValues > 0.975 | pValues < 0.025;
%IS THIS A SUM OR A MEAN?? IF SUM, IT WOULD BE W/M. IF MEAN, IT WOULD BE
%W/M^2, RIGHT??
%         FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2).*LatWeights(:,2);
        FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2);
        TempFlux = FluxSumsAllLats.(IndexNames{i}).(FluxNames{j});
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
    legend(FluxNames(5:end),'Location','EastOutside')
%     uistack(netHandle, 'top') %MUST come AFTER legend
    title(['Total Latitudinal Regre. of ', IndexNames{i}, ' (W/m). 95% T-Test Sig Marked with diamonds'])
    xlabel('latitude')
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['TotalLatRegress',IndexNames{i},'.png']);
    hold off;
end

[beta,pvalue]=RegressLatRegionvsTimeTTestContours(MonthlyFluxDepartures.Net,Indices.NAM,'Net','NAM',30,60)
[beta,pvalue]=RegressLatRegionvsTimeTTestContours(MonthlyFluxDepartures.Net,Indices.NAM,'Net','NAM',50,90)


[beta,pvalue]=RegressLatRegionvsTimeTTestContours(MonthlyFluxDepartures.Net,Indices.Ice,'Net','Ice',70,90)


