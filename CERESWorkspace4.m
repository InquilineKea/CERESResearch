load LatitudinalDifFluxTimeSeriesAsField.mat

% FluxLatitude.SW = [FluxTimeSeries.SW.Lat49_90.S;FluxTimeSeries.SW.Lat30_49.S;...
%     FluxTimeSeries.SW.Lat15_30.S;FluxTimeSeries.SW.Lat0_15.S;FluxTimeSeries.SW.Lat0_15.N;
%     FluxTimeSeries.SW.Lat15_30.N;FluxTimeSeries.SW.Lat30_49.N;FluxTimeSeries.SW.Lat49_90.N;...
%     FluxTimeSeries.SW.AllLats.S;FluxTimeSeries.SW.AllLats.N;FluxTimeSeries.SW.Lat49_90.Dif;...
%     FluxTimeSeries.SW.Lat30_49.Dif;FluxTimeSeries.SW.Lat15_30.Dif;FluxTimeSeries.SW.Lat0_15.Dif];

fields = fieldnames(FluxTimeSeries);
% FluxTimeSeries.(fields{1})
% FluxTimeSeries.(fields{1}).Lat49_90.Dif
% 
for i=1:length(fields)
   FluxLatitude.(fields{i}) = [FluxTimeSeries.(fields{i}).Lat49_90.S;FluxTimeSeries.(fields{i}).Lat30_49.S;...
    FluxTimeSeries.(fields{i}).Lat15_30.S;FluxTimeSeries.(fields{i}).Lat0_15.S;FluxTimeSeries.(fields{i}).Lat0_15.N;
    FluxTimeSeries.(fields{i}).Lat15_30.N;FluxTimeSeries.(fields{i}).Lat30_49.N;FluxTimeSeries.(fields{i}).Lat49_90.N;...
    FluxTimeSeries.(fields{i}).AllLats.S;FluxTimeSeries.(fields{i}).AllLats.N;FluxTimeSeries.(fields{i}).Lat49_90.Dif;...
    FluxTimeSeries.(fields{i}).Lat30_49.Dif;FluxTimeSeries.(fields{i}).Lat15_30.Dif;FluxTimeSeries.(fields{i}).Lat0_15.Dif];
end

% LatBand = fieldnames(FluxTimeSeries.(fields{1}));
% LatBand2 = cellfun(@strcat,LatBand,'S')

load Indices.mat
Indices = rmfield(Indices,'PNA');
Indices = rmfield(Indices,'SAO');
Indices = rmfield(Indices,'T');
IndexNames = fieldnames(Indices);

corr_mat=corr([FluxLatitude.(fields{1}); NAO']');
corr_mat(15,:);

corr_mat=corr([FluxLatitude.(fields{1}); Indices.IndexNames{1}]');
pcolor(corr_mat);colorbar

for i =1:length(fields)
    for j=1:length(IndexNames)
    corr_mat=corr([FluxLatitude.(fields{i}); Indices.(IndexNames{j})]');
    plot(corr_mat(1:end-1,end),'rs','MarkerFaceColor','g','markersize', 20)
    set(gca,'FontSize',18)
    title(['Correlation Coefficients between each Latitudinal Band in ' fields{i},' and ', IndexNames{j} ,' index'])
    xlabel('Latitudinal Band')
    ylabel('Correlation Coefficient')
    grid on;
    set(gca,'xtick',1:1:15)
    set(gca,'XTickLabel',{'90-49S','49-30S','30-15S','15-0S','0-15N','15-30N','30-49N','49-90N','SH','NH','90-49Dif','49-30Dif','30-15Dif','15-0Dif'})
    set(gca,'GridLineStyle','--')
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['LatitudeCorrs_', fields{i}, '_and_', IndexNames{j},'.png']);
    hold off;
    end
end

load('ClimatologySubtractedVars','SWClimatologySubtracted');
load('ClimatologySubtractedVars','TempClimatologySubtracted');
load('ClimatologySubtractedVars','PrecipClimatologySubtracted');
load('ClimatologySubtractedVars','SWCFClimatologySubtracted');

GlobalLatCorrelationPlot(SWClimatologySubtracted);
GlobalLatCorrelationPlot(TempClimatologySubtracted);
GlobalLatCorrelationPlot(PrecipClimatologySubtracted);
GlobalLatCorrelationPlot(SWCFClimatologySubtracted);


load MonthlyFluxDeparturesStruct.mat

cmap = jet;

FluxNames = fieldnames(MonthlyFluxDepartures);
for i=1:length(FluxNames)
    LatMeans = squeeze(mean(MonthlyFluxDepartures.(FluxNames{i}),2));
    for j=1:length(IndexNames)
        plot(corr(LatMeans', Indices.(IndexNames{j})'),'color',cmap(j*floor(64/ceil(length(IndexNames))),:),'LineWidth',3);
        grid on;
        set(gca,'GridLineStyle','--')
        set(gca,'FontSize',18)
        set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{i}),2),1))/180:180.5)
        set(gca,'xticklabel',num2cell(-90:10:90))
        xlabel('latitude')
        ylabel('correlation')
        title(['Correlation between ',FluxNames{i},' and each global index'])
        hold on;
    end
    legend(IndexNames)
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['AllLatIndexCorrs_',FluxNames{i},'.png']);
    hold off;
end
for i=1:length(FluxNames)
    GlobalLatCorrelationPlot(MonthlyFluxDepartures.(FluxNames{i}),FluxNames{i})
end

for i=1:length(FluxNames)
    for j=1:length(IndexNames)
        PlotIndexThreshold(MonthlyFluxDepartures.(FluxNames{i}),Indices.(IndexNames{j}),90,FluxNames{i},IndexNames{j});
    end
end

for i=1:length(FluxNames)
    for j=1:length(IndexNames)
        PlotIndexThresholdLower(MonthlyFluxDepartures.(FluxNames{i}),Indices.(IndexNames{j}),10,FluxNames{i},IndexNames{j});
    end
end

blah=squeeze(mean(MonthlyFluxDepartures.(FluxNames{5}),2));

regress(blah(1,:)',[Indices.NAM' Indices.NAO' Indices.NINO34' Indices.NAM'.*Indices.NINO34'])
regress(blah(1,:)',[Indices.NAM' Indices.NAO' Indices.NINO34' ])

% arrayfun(@regress,blah',[Indices.NAM' Indices.NAO' Indices.NINO34' ])
% % arrayfun(@regress,blah',repmat(Indices.NAM,180,1))
% 
% X = [Indices.NAM' Indices.SAM' Indices.NAO' Indices.NINO34'];
%Beta = inv(X'*X)*X'*blah(1,:)';

%     Beta = inv(X'*X)*X'*blah';
%     plot(Beta')
%     legend(IndexNames(1:4))
%     set(gca,'xtick',1:10:181)
%     set(gca,'xticklabel',num2cell(-90:10:90))
%     plot(abs(Beta'))
%     %X = [Indices.NAM' Indices.SAM' Indices.NAO' Indices.NINO34' Indices.NAM'.*Indices.NINO34' Indices.SAM'.*Indices.NINO34' Indices.NAM'.*Indices.NAO'];
%     X = [Indices.NAM' Indices.SAM' Indices.NAO' Indices.NINO34' Indices.NAM'.*Indices.NINO34'];
%     %Beta = inv(X'*X)*X'*blah(1,:)';
for i=1:length(FluxNames)
    Beta = inv(X'*X)*X'*squeeze(mean(MonthlyFluxDepartures.(FluxNames{i}),2))';   
    plot(Beta','LineWidth',3)
    grid on;
    IndexLegend = IndexNames(1:4);
    IndexLegend{5} = 'NAM*NINO';
    legend(IndexLegend,'Location','NorthWestOutside')
    set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{i}),2),1))/180:180.5)
    set(gca,'xticklabel',num2cell(-90:10:90))
    set(gca,'FontSize',18)
    title(['Latitudinal Regression Coefficient of each Index for ', FluxNames{i}])
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['IndexLatRegressions_',FluxNames{i},'.png']);
    hold off;
end
%Warning: Matrix is close to singular or badly scaled.
%Results may be inaccurate. RCOND =  7.016633e-17. 
%^is this error due to NAM+NINO b/c that is NOT linearly independnet and
%multiple regression DEPENDS on linear independence?

% X = [Indices.NAM' Indices.SAM' Indices.NAO' Indices.NINO34' Indices.PNA'];
% %Beta = inv(X'*X)*X'*blah(1,:)';
% Beta = inv(X'*X)*X'*blah';
% plot(Beta','LineWidth',3)
% legend(IndexNames(1:4),'PNA')
% set(gca,'xtick',1:10:181)
% set(gca,'xticklabel',num2cell(-90:10:90))

blah = reshape(MonthlyFluxDepartures.(FluxNames{5}),180*360,156);
X = [Indices.NAM' Indices.SAM' Indices.NAO' Indices.NINO34' Indices.PNA' Indices.NAM'.*Indices.NINO34' Indices.SAM'.*Indices.NINO34' Indices.NAM'.*Indices.NAO'];
%Beta = inv(X'*X)*X'*blah(1,:)';
Beta = inv(X'*X)*X'*blah';
Beta2=reshape(Beta,size(X,2),180,360);
Beta2=permute(Beta2,[2 3 1]);
load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(Beta2(:,:,1),geoidrefvec,'DisplayType','texturemap');colorbar

% 
% arrayfun(@(x)regress(x,Indices.NAM') ,blah')
% bsxfun(@times,inv(X'*X)*X',blah')
% 
% 
% NewRegress = @(x) regress(x,[Indices.NAM' Indices.NAO' Indices.NINO34' Indices.NAM'.*Indices.NINO34']);
% arrayfun(NewRegress,blah)

subplot(2,1,1)

for i=1:length(FluxNames)
    contourf(squeeze(mean(MonthlyFluxDepartures.(FluxNames{i}),2)));colorbar%colorbar('LOCATION','SouthOutside')
    caxis([-max(max(abs(squeeze(mean(MonthlyFluxDepartures.(FluxNames{i}),2))))) max(max(abs(squeeze(mean(MonthlyFluxDepartures.(FluxNames{i}),2)))))])
    grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
    set(gca,'xtick',11:12:size(MonthlyFluxDepartures.(FluxNames{i}),3))
    set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
    xlabel('Year');
    set(gca,'ytick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{i}),2),1))/180:180.5)
    set(gca,'yticklabel',num2cell(-90:10:90))
    ylabel('latitude')
    hold on;
    title(['Time-Latitudinal Contour Plot of ', FluxNames{i}])
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['TimeLatContours_',FluxNames{i},'.png']);
    hold off;
end

subplot(2,1,2)
cmap = jet;
for j = 1:4
    plot(Indices.(IndexNames{j})','color',cmap(j*floor(64/ceil(4)),:),'LineWidth',3)
    hold on;
    grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
end
set(gca,'xtick',11:12:size(MonthlyFluxDepartures.LW,3))
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
legend(IndexNames)
set(gca,'xtick',11:12:size(MonthlyFluxDepartures.LW,3))
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year');
ylabel('Index');
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['AllIndicesOverTime','.png']);

%%%%%%%%%%
% contour3(CorMatrix)
%[c,h]=contourf(CorMatrix);colorbar
LatMeans = squeeze(mean(SWClimatologySubtracted,2));

% contour3(SWClimatologySubtracted)
CorMatrix = corr([LatMeans']);

contourf(CorMatrix);colorbar
grid on;
 set(gca,'xtick',1:10:181)
 set(gca,'xticklabel',num2cell(-90:10:90))
  set(gca,'ytick',1:10:181)
 set(gca,'yticklabel',num2cell(-90:10:90))
caxis([-1 1])
%clabel(c,h)

% pcolor(CorMatrix);colorbar
%  set(gca,'xtick',1:10:181)
%  set(gca,'xticklabel',num2cell(-90:10:90))
%   set(gca,'ytick',1:10:181)
%  set(gca,'yticklabel',num2cell(-90:10:90))
% caxis([-1 1])
% 
% surf(CorMatrix);colorbar
%  set(gca,'xtick',1:10:181)
%  set(gca,'xticklabel',num2cell(-90:10:90))
%   set(gca,'ytick',1:10:181)
%  set(gca,'yticklabel',num2cell(-90:10:90))
% caxis([-1 1])




% set(gca,'xtick',sind(-90:30:90))
% set(gca,'xticklabel',{'90S','60S','30S','0','30N','60N','90N'})

set(gca,'xtick',1:10:180)
set(gca,'XTickLabel',{'90-49S','49-30S','30-15S','15-0S','0-15N','15-30N','30-49N','49-90N','SH','NH','90-49Dif','49-30Dif','30-15Dif','15-0Dif'})

