
load LatitudinalDifFluxTimeSeriesAsField.mat
load Indices.mat
Indices = rmfield(Indices,'PNA');
Indices = rmfield(Indices,'SAO');
Indices = rmfield(Indices,'T');
IndexNames = fieldnames(Indices);

for i = 1:size(FluxNames)
NHMinusSH(MonthlyFluxDepartures.(FluxNames{i}),156,0,90,FluxNames{i});
end

for i = 1:size(FluxNames)
HemisphericFluxesClimSubtracted(MonthlyFluxDepartures.(FluxNames{i}),0,90,FluxNames{i});
end

for i = 1:size(FluxNames)
ContourLatPlot(MonthlyFluxDepartures.(FluxNames{i}),FluxNames{i});
end
%%%%%%%%%%%


load Indices.mat
Indices = rmfield(Indices,'PNA');
Indices = rmfield(Indices,'SAO');
Indices = rmfield(Indices,'T');
IndexNames = fieldnames(Indices);


days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
allMonthWeights = repmat(days_per_month,1,12);
allMonthWeights = [31 30 31 30 31 31 30 31 30 31 allMonthWeights 31 28];
allMonthWeights(48) = 29; %leap year, http://www.wolframalpha.com/input/?i=months+between+march+2000+and+february+2004
allMonthWeights(96) = 29;


moving_sum = @(n, x) filter(ones(1,n), 1, x);
IndicesMWA.NAM= bsxfun(@rdivide,moving_sum(12, Indices.NAM .* allMonthWeights) ...
    ,moving_sum(12,allMonthWeights)); %or Mean Lat Flux here

cmap = jet;
for i = 1:length(IndexNames)
    IndicesMWA.(IndexNames{i})= bsxfun(@rdivide,moving_sum(12, Indices.(IndexNames{i}) .* allMonthWeights) ...
    ,moving_sum(12,allMonthWeights)); %or Mean Lat Flux here
    Temp = IndicesMWA.(IndexNames{i});
    IndicesMWA.(IndexNames{i}) = Temp(12:end);
    IndexMWA2Plot = IndicesMWA.(IndexNames{i});
    plot(IndexMWA2Plot(12:end),'color',cmap(i*floor(64/ceil(length(IndexNames))),:),'LineWidth',3)
    grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
    set(gca,'FontSize',20)
    xlabel('Year End')
    ylabel('Index MWA Anomaly')
    set(gca,'xtick',[11:12:119+26])
    set(gca,'XTickLabel',{2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
legend(IndexNames,'Location','NorthEastOutside')
    hold on;
    title('MWA in Indices over Time')
end

set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['MWA_Indices','.png']);

load('MWA.mat')

[MWABeta_Net.NAM, blah, blah, blah,Alpha_Net.NAM] = Regress1DTimeSeriesMapWithTTestContours(MWA.Net,IndicesMWA.NAM,'Net MWA','NAM MWA')
[MWABeta_Net.SAM, blah, blah, blah] = Regress1DTimeSeriesMapWithTTestContours(MWA.Net,IndicesMWA.SAM,'Net MWA','SAM MWA')
[MWABeta_Net.NINO34, blah, blah, blah] = Regress1DTimeSeriesMapWithTTestContours(MWA.Net,IndicesMWA.NINO34,'Net MWA','NINO34 MWA')

MWATimeLatContourRegression(Alpha_Net.NAM,MWABeta_Net.NAM,MWA.Net,IndicesMWA.NAM,'Net MWA','NAM MWA',12);
MWATimeLatContourRegression(MWABeta_Net.SAM,IndicesMWA.SAM,'Net MWA','SAM MWA');
MWATimeLatContourRegression(MWABeta_Net.NINO34,IndicesMWA.NINO34,'Net MWA','NINO34 MWA');
%%%%%%%%%%%

load MonthlyFluxDeparturesStruct.mat
MonthlyFluxDepartures = rmfield(MonthlyFluxDepartures,'Clouds');
MonthlyFluxDepartures = rmfield(MonthlyFluxDepartures,'Precip');
FluxNames = fieldnames(MonthlyFluxDepartures);


HemisphericFluxesClimSubtracted(MonthlyFluxDepartures.Net,0,90,'Net')
HemisphericFluxesClimSubtracted(MonthlyFluxDepartures.SW,0,90,'SW')
HemisphericFluxesClimSubtracted(MonthlyFluxDepartures.LW,0,90,'LW')
ContourLatPlot(MonthlyFluxDepartures.Net,'Net')


HemisphericFluxesClimSubtracted(MonthlyFluxDepartures.Net,0,15,'Net')
HemisphericFluxesClimSubtracted(MonthlyFluxDepartures.SW,0,15,'SW')
HemisphericFluxesClimSubtracted(MonthlyFluxDepartures.LW,0,15,'LW')


for i = 1:length(FluxNames)
   
    ContourLatPlot(MonthlyFluxDepartures.(FluxNames{i}),FluxNames{i})

    
end

NHMinusSH(Net,0,90,'Net')

NHMinusSH(MonthlyFluxDepartures.Net,0,15,'Net')

NHMinusSH(MonthlyFluxDepartures.Net,45,90,'Net')

NHMinusSH(MonthlyFluxDepartures.SW,0,15,'SW')
NHMinusSH(MonthlyFluxDepartures.LW,0,15,'LW')



NHMinusSH(SW,0,90,'SW')

NHMinusSH(LW,0,90,'LW')




