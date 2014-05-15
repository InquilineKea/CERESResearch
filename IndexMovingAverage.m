function IndicesMWA = IndexMovingAverage(MonthFilterSize)
load Indices.mat
Indices = rmfield(Indices,'PNA');
Indices = rmfield(Indices,'SAO');
Indices = rmfield(Indices,'T');
Indices = rmfield(Indices,'NAO');

IndexNames = fieldnames(Indices);

days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
allMonthWeights = repmat(days_per_month,1,12);
allMonthWeights = [31 30 31 30 31 31 30 31 30 31 allMonthWeights 31 28];
allMonthWeights(48) = 29; %leap year, http://www.wolframalpha.com/input/?i=months+between+march+2000+and+february+2004
allMonthWeights(96) = 29;
moving_sum = @(n, x) filter(ones(1,n), 1, x);
cmap = jet;
for i = 1:length(IndexNames)
    IndicesMWA.(IndexNames{i})= bsxfun(@rdivide,moving_sum(MonthFilterSize, Indices.(IndexNames{i}) .* allMonthWeights) ...
    ,moving_sum(MonthFilterSize,allMonthWeights)); %or Mean Lat Flux here
    Temp = IndicesMWA.(IndexNames{i});
    IndicesMWA.(IndexNames{i}) = Temp(MonthFilterSize:end);
    IndexMWA2Plot = IndicesMWA.(IndexNames{i});
    plot(IndexMWA2Plot,'color',cmap(i*floor(64/ceil(length(IndexNames))),:),'LineWidth',3)
    set(gca,'FontSize',20)
    xlabel('Year End')
    ylabel('Index MWA Anomaly')
    set(gca,'xtick',[12-2-(MonthFilterSize-1):12:length(IndicesMWA.(IndexNames{i}))]) %what if Filter Size is 12 months?
    set(gca,'XTickLabel',2000:2012)
        grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')

legend(IndexNames,'Location','NorthEastOutside')
    hold on;
    title([num2str(MonthFilterSize),' Month Moving Average in Indices over Time'])
end

% 12-2-(MonthFilterSize-1):12:length(IndicesMWA.(IndexNames{i}))
% length(IndicesMWA.(IndexNames{i}))
%start is 2000-03.
%if MovingFilterSize = 1, then 12-3+1 months. if MovingFilterSize = 6,
%12-3+1-6 months
%if =6, then first end-date is August (5 months from Aug to Dec). 12-2-6+1.

set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[num2str(MonthFilterSize),'Month_Moving_Average_Indices','.png']);
hold off;

end