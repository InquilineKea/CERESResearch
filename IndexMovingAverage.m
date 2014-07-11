function IndicesMWA = IndexMovingAverage(MonthFilterSize)

NAM = dlmread('NAM.ascii');
% NAO = dlmread('NAO.ascii');
SAM = dlmread('SAM.ascii');
% PNA = dlmread('PNA.ascii');
% T = dlmread('monthly.land_ocean.90S.90N.df_1901-2000mean.dat');
%NAM(find(NAM(:,1) >= 2000),3)
%find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2013 & NAM(:,2)==2)
NAM = NAM(find(NAM(:,1) == 2000 & NAM(:,2)==3):find(NAM(:,1) == 2014 & NAM(:,2)==2),3);
% NAO = NAO(find(NAO(:,1) == 2000 & NAO(:,2)==3):find(NAO(:,1) == 2013 & NAO(:,2)==2),3);
SAM = SAM(find(SAM(:,1) == 2000 & SAM(:,2)==3):find(SAM(:,1) == 2014 & SAM(:,2)==2),3);
% PNA = PNA(find(PNA(:,1) == 2000 & PNA(:,2)==3):find(PNA(:,1) == 2013 & PNA(:,2)==2),3);
% T = T(find(T(:,1) == 2000 & T(:,2)==3):find(T(:,1) == 2013 & T(:,2)==2),3);
ensoTime = ncread('enso.cdf','T');
NINO34 = ncread('enso.cdf','NINO34');
NINO34 = NINO34(find(ensoTime > 482 & ensoTime < 638+12));

Indices.NAM = NAM';
Indices.SAM = SAM';
% Indices.NAO = NAO';
Indices.NINO34 = NINO34';
% Indices.PNA = PNA';
% Indices.SAO = double(SAO');
% Indices.T = T';
save('Indices.mat','Indices')

% load Indices.mat
% Indices = rmfield(Indices,'PNA');
% Indices = rmfield(Indices,'SAO');
% Indices = rmfield(Indices,'T');
% Indices = rmfield(Indices,'NAO');

IndexNames = fieldnames(Indices);

days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
allMonthWeights = repmat(days_per_month,1,13);
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
    set(gca,'XTickLabel',2000:2013)
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