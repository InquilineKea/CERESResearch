
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
FluxTimeSeries.Net.AllLats.N

NormIndices.NAM = Indices.NAM/std(Indices.NAM);
NormIndices.SAM= Indices.SAM/std(Indices.SAM);
NormIndices.NAO = Indices.NAO/std(Indices.NAO);
NormIndices.NINO34= Indices.NINO34/std(Indices.NINO34);


Flux1DName = fieldnames(FluxTimeSeries);
LatName = fieldnames(FluxTimeSeries.(Flux1DName{1}));
HemisphereName = fieldnames(FluxTimeSeries.(Flux1DName{1}).(LatName{1}));


for i=1:length(IndexNames) %regress all fields [and indices] over total NH net flux and total SH net flux. might do it for total NH+SH flux too.
    for j=1:length(FluxNames)
    [RegressIndexFlux.(IndexNames{i}).(FluxNames{j}), LatRegressFluxDif.(IndexNames{i}).(FluxNames{j}),FluxSumsAllLats.(IndexNames{i}).(FluxNames{j})] = Regress1DTimeSeriesMap(MonthlyFluxDepartures.(FluxNames{j}), Indices.(IndexNames{i}),FluxNames{j},IndexNames{i});
    %NH Total, SH Total, NH-SH Total, NH+SH Total

    %WorldRegress2 = Regress1DTimeSeriesMap(MonthlyFluxDepartures.(FluxNames{i}), FluxTimeSeries.Net.AllLats.N,MonthlyFluxDepartures.(FluxNames{i}),FluxNames{i},[Flux1DName{1}, LatName{5}, HemisphereName{1},'H']);
    % clear WorldRegress1;
    % clear WorldRegress2
    end
end
save('RegressIndexFlux.mat','RegressIndexFlux')
%reconstruct 180x360x156 matrix of how index contributes to a given flux
%need to convert 180x360 into 180x360x156
load LatWeights
load RegressIndexFlux
for i=1:length(IndexNames) %regress all fields [and indices] over total NH net flux and total SH net flux. might do it for total NH+SH flux too.
    for j=1:length(FluxNames)
    RegressIndexFluxAllTime.(IndexNames{i}).(FluxNames{j})=bsxfun(@times,RegressIndexFlux.(IndexNames{i}).(FluxNames{j}), permute(Indices.(IndexNames{i})',[3 2 1]));
    %then lat-time contour plot it
    [X,Y]=meshgrid(1:156,sind(LatWeights(:,1)));
    ToContour = squeeze(mean(RegressIndexFluxAllTime.(IndexNames{i}).(FluxNames{j}),2));
    contourf(X,Y,ToContour);colorbar
%     figure(2)
%     contourf(squeeze(mean(RegressIndexFluxAllTime.(IndexNames{i}).(FluxNames{j}),2)));
%     blah=squeeze(mean(RegressIndexFluxAllTime.(IndexNames{i}).(FluxNames{j}),2));
    caxis([-max(max(abs(ToContour))) max(max(abs(ToContour)))])
    MakeLizMap
    colormap(lizmap)
    grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
    set(gca,'xtick',11:12:size(MonthlyFluxDepartures.(FluxNames{i}),3))
    set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
    xlabel('Year');
%     plot(sind(LatWeights(:,1)))
%     plot(-1:1)
    set(gca,'ytick',sind((0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{i}),2),1))/180:179.5)-90))
    set(gca,'yticklabel',num2cell(-90:10:90))
    ylabel('latitude')
    title(['Time-Latitudinal Contour Plot of ', FluxNames{j},' over ',IndexNames{i} ' (W/m^2)'])
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['TimeLatRegressContours',FluxNames{j},'_',IndexNames{i},'.png']);
    hold off;
    clear RegressIndexFluxAllTime.(IndexNames{i}).(FluxNames{j});
    end
end

save('RegressIndexFluxAllTimeNew.mat','RegressIndexFluxAllTime',VERSION,'-v7.3')
load RegressIndexFluxAllTimeNew
%Weighted Contour Plot Contributions
RegressIndexFluxAllTime.NAM = NAM;
clear NAM;
RegressIndexFluxAllTime.NAO = NAO;
clear NAO;
RegressIndexFluxAllTime.NINO34=NINO34;
clear NINO34;
RegressIndexFluxAllTime.SAM=SAM;
clear SAM;
for i=1:length(IndexNames) %regress all fields [and indices] over total NH net flux and total SH net flux. might do it for total NH+SH flux too.
    for j=1:length(FluxNames)
    PlotAsContour = bsxfun(@times,squeeze(mean(RegressIndexFluxAllTime.(IndexNames{i}).(FluxNames{j}),2)),LatWeights(:,2));
    contourf(PlotAsContour);colorbar
    caxis([-max(max(abs(PlotAsContour))) max(max(abs(PlotAsContour)))])
    MakeLizMap
    colormap(lizmap)
    grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
    set(gca,'xtick',11:12:size(PlotAsContour,2))
    set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
    xlabel('Year');
    set(gca,'ytick',0.5:10*(size(PlotAsContour,1)/180):180.5)
    set(gca,'yticklabel',num2cell(-90:10:90))
    ylabel('latitude')
    title(['WeightedTime-Latitudinal Contour Plot of ', FluxNames{j},' over ',IndexNames{i}])
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['WeightedTimeLatRegressContours',FluxNames{j},'_',IndexNames{i},'.png']);
    hold off;
    end
end
clear RegressIndexFluxAllTime

%also, should do regressions latitudinally and integrate the regression
%data over the Hemisphere to see if their field helps contribute to the
%total NH net flux (should multiply that by the std of total NH net flux)
%this can't justbe regressed on net flux though. it should also be
%regressed on all other types of fluxes (and flux differences).
%after all we want to know if an index/flux ultimately contributes to NH-SH
%difference and HOW
%[MonthlyFluxDeparture.FluxNames{i}.HemisBudgetsFromRegressOnNet, MonthlyFluxDeparture.FluxNames{i}.NetLatitudeFluxDif] = Regress1DTimeSeriesMap(MonthlyFluxDepartures.(FluxNames{i}), FluxTimeSeries.Net.AllLats.N,MonthlyFluxDepartures.(FluxNames{i}),FluxNames{i},[Flux1DName{1}, LatName{5}, HemisphereName{1},'H']);
%and eventually find a time-dependent way to show the evolution of the
%fields over time. i mean.. you just get regression coefficients for the
%whole world. then multiply these regression coefficients over time to get
%time series of size 180x360x156
% plot(Indices.NetLatitudeFluxDif{2})
% plot(Indices.NetLatitudeFluxDif{3})
% plot(Indices.NetLatitudeFluxDif{4})

for i=1:4
    for j=1:length(FluxNames)
        plot(FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}))
        grid on;
        set(gca,'FontSize',18)
        title(['Total Latitudinal Regre. of ',FluxNames{j},' over ' IndexNames{i}])
        set(gca,'GridLineStyle','--')
        set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{j}),2),1))/180:180.5)
        set(gca,'xticklabel',num2cell(-90:10:90))
        xlabel('latitude')
        set(gcf,'paperposition',[0 0 20 10])
        print(gcf,'-dpng','-r300',['TotalLatRegress',FluxNames{j},'Over',IndexNames{i},'.png']);
        hold off;
    end
end
% for i=1:4
%     for j=1:size(FluxNames)
%         temp = FluxSumsAllLats.(IndexNames{i}).(FluxNames{j});
%         plot(cumsum(flipud(temp(1:end/2))),'LineWidth',5)%SH: -90 to 0
%         set(gca,'FontSize',18)
%         hold on;
%         grid on;
%         set(gca,'GridLineStyle','--')
%         plot(cumsum(temp(end/2+1:end)),'g','LineWidth',5) %NH: 90 to 0
%         legend({'SH','NH'})
%         title(['Cum. Sum (eq. to pole) of Latitudinal Regre. of ',FluxNames{j},' over ' IndexNames{i}])
%         set(gcf,'paperposition',[0 0 20 10])
%         print(gcf,'-dpng','-r300',['CumSumLatRegress',FluxNames{j},'Over',IndexNames{i},'.png']);
%         hold off;
%     end
% end

% for i=1:4
%     for j=1:size(FluxNames)
%         temp = FluxSumsAllLats.(IndexNames{i}).(FluxNames{j});
%         plot(cumsum(flipud(temp(1:end/2))))%SH: -90 to 0
%         set(gca,'FontSize',18)
%         hold on;
%         grid on;
%         set(gca,'GridLineStyle','--')
%         plot(cumsum(temp(end/2+1:end)),'g') %NH: 90 to 0
%         legend({'SH','NH'})
%         title(['Cum. Sum (eq. to pole) of Latitudinal Regre. of ',FluxNames{j},' over ' IndexNames{i}])
%         set(gcf,'paperposition',[0 0 20 10])
%         print(gcf,'-dpng','-r300',['CumSumLatRegress',FluxNames{j},'Over',IndexNames{i},'.png']);
%         hold off;
%     end
% end

load('RegressIndexFlux.mat')
load LatWeights.mat

%test = structfun(@(x)(orderfields(x)),RegressIndexFlux);
RegressIndexFlux.NAM = orderfields(RegressIndexFlux.NAM);
RegressIndexFlux.SAM = orderfields(RegressIndexFlux.SAM);
RegressIndexFlux.NAO= orderfields(RegressIndexFlux.NAO);
RegressIndexFlux.NINO34= orderfields(RegressIndexFlux.NINO34);
FluxNames = fieldnames(RegressIndexFlux.NINO34);

Style = repmat({'--','-'},1,length(FluxNames));
Markers = repmat({'none','none'},1,length(FluxNames));
for i=1:4
    for j=1:length(FluxNames)
        lats = size(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),1);
        Factor = 180/lats;
%IS THIS A SUM OR A MEAN?? IF SUM, IT WOULD BE W/M. IF MEAN, IT WOULD BE
%W/M^2, RIGHT??
%         FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2).*LatWeights(:,2);
        FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}) = sum(resizem(RegressIndexFlux.(IndexNames{i}).(FluxNames{j}),Factor),2);
        if strcmp(FluxNames{j},'Net')
            %h = plot(FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{j},'Marker',Markers{j},'LineWidth',9);
            h = plot(sind(LatWeights(:,1)),FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{j},'Marker',Markers{j},'LineWidth',9);
        else
            %plot(FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{j},'Marker',Markers{j},'LineWidth',3)    
            plot(sind(LatWeights(:,1)),FluxSumsAllLats.(IndexNames{i}).(FluxNames{j}),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineStyle',Style{j},'Marker',Markers{j},'LineWidth',3)
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
        uistack(h, 'top') %MUST come AFTER legend
        title(['Total Latitudinal Regre. of ', IndexNames{i}, ' (W/m)'])
        xlabel('latitude')
        set(gcf,'paperposition',[0 0 20 10])
        print(gcf,'-dpng','-r300',['TotalLatRegress',IndexNames{i},'.png']);
        hold off;
end

%%%try this again with significance markers for latitudes of regression
%%%t-test significance.

regression(squeeze(sum(MonthlyFluxDepartures.Net,2)),Indices.NAM)



for i=1:4
    for j=1:9
        temp = FluxSumsAllLats.(IndexNames{i}).(FluxNames{j});
        %plot(cumsum(flipud(temp(1:end/2))./cumsum(LatWeights(end/2+1:end,2))),'LineWidth',5)%SH: -90 to 0
        plot(cumsum(flipud(temp(1:end/2))),'LineWidth',5)%SH: -90 to 0
        %cumsum in denominator doesn't work like that!
        set(gca,'FontSize',18)
        hold on;
        grid on;
        set(gca,'GridLineStyle','--')
       % plot(cumsum(temp(end/2+1:end)./cumsum(LatWeights(end/2+1:end,2))),'g','LineWidth',5) %NH: 90 to 0
        plot(cumsum(temp(end/2+1:end)),'g','LineWidth',5) %NH: 90 to 0

        legend({'SH','NH'})
        title(['Cum. Sum (eq. to pole) of Latitudinal Regre. of ',FluxNames{j},' over ' IndexNames{i}])
        set(gcf,'paperposition',[0 0 20 10])
        print(gcf,'-dpng','-r300',['CumSumLatRegress',FluxNames{j},'Over',IndexNames{i},'.png']);
        hold off;
    end
end

for i=1:4
    for j=1:size(FluxNames)
        RegressMean = mean(mean(resizem(Beta,Factor),2).*LatWeights(:,2)); 
        FluxSumsAllLatitudeLines = sum(resizem(Beta,Factor),2).*LatWeights(:,2);
        NHPlusSHTotalFluxContribution = sum(FluxSumsAllLatitudeLines);
        NHTotalFluxContribution = sum(FluxSumsAllLatitudeLines(end/2+1:end));
        SHTotalFluxContribution = sum(FluxSumsAllLatitudeLines(1:end/2));
        NHMinusSHTotalFluxContribution = NHTotalFluxContribution - SHTotalFluxContribution;
        LatRegressFluxDif.(IndexNames{i}).(FluxNames{j}) = LatRegressFluxDif.(IndexNames{i}).(FluxNames{j});
        plot(LatRegressFluxDif.(IndexNames{i}).(FluxNames{j}),'color',cmap(j*floor(64/ceil(length(FluxNames))),:),'LineWidth',3)
        grid on;
        set(gca,'FontSize',18)
        set(gca,'GridLineStyle','--')
        set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{i}),2),1))/180:180.5)
        set(gca,'xticklabel',num2cell(0:10:90))
    end
    title(['Latitudinal Regre. Dif of ', IndexNames{i}])
    xlabel('latitude')
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['DifLatRegress',IndexNames{i},'.png']);
    hold off;
end

for i=1:length(FluxNames) %regress all fields [and indices] over total NH net flux and total SH net flux. might do it for total NH+SH flux too.
    for j=i+1:length(FluxNames)
    [FluxFluxRegressBeta.(FluxNames{i}).(FluxNames{j}), FluxFluxLatRegressDif.(FluxNames{i}).(FluxNames{j}),FluxFluxSumsAllLats.(FluxNames{i}).(FluxNames{j})] = GlobalRegressionMap(MonthlyFluxDepartures.(FluxNames{i}), MonthlyFluxDepartures.(FluxNames{j}),FluxNames{i},FluxNames{j});
    end
end

save('FluxFluxRegressBeta.mat','FluxFluxRegressBeta')

for i=1:length(FluxNames)
    for j=i+1:length(FluxNames)
        plot(FluxFluxSumsAllLats.(FluxNames{i}).(FluxNames{j}))
        grid on;
        set(gca,'FontSize',18)
        title(['Total Latitudinal Regre. of ',FluxNames{i},' over ' FluxNames{j}])
        set(gca,'GridLineStyle','--')
        length(FluxFluxSumsAllLats.(FluxNames{i}).(FluxNames{j})) %for xtick when u want to do precip
        set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{j}),2),1))/180:180.5)
        set(gca,'xticklabel',num2cell(-90:10:90))
        xlabel('latitude')
        set(gcf,'paperposition',[0 0 20 10])
        print(gcf,'-dpng','-r300',['TotalLatRegress',FluxNames{i},'Over',FluxNames{j},'.png']);
        hold off;
    end
end
for i=1:length(FluxNames)
    for j=i+1:length(FluxNames)
        temp = FluxFluxSumsAllLats.(FluxNames{i}).(FluxNames{j});
        plot(cumsum(flipud(temp(1:end/2))))%SH: -90 to 0
        set(gca,'FontSize',18)
        hold on;
        grid on;
        set(gca,'GridLineStyle','--')
        plot(cumsum(temp(end/2+1:end)),'g') %NH: 90 to 0
        legend({'SH','NH'})
        title(['Cum. Sum (eq. to pole) of Latitudinal Regre. of ',FluxNames{i},' over ' FluxNames{j}])
        set(gcf,'paperposition',[0 0 20 10])
        print(gcf,'-dpng','-r300',['CumSumLatRegress',FluxNames{i},'Over',FluxNames{j},'.png']);
        hold off;
    end
end
for i=1:length(FluxNames)
    for j=i+1:length(FluxNames)
        plot(FluxFluxLatRegressDif.(FluxNames{i}).(FluxNames{j}))
        grid on;
        set(gca,'FontSize',18)
        set(gca,'GridLineStyle','--')
        set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{i}),2),1))/180:180.5)
        set(gca,'xticklabel',num2cell(0:10:90))
        title(['Latitudinal Regre. Dif of ',FluxNames{i},' over ' FluxNames{j}])
        xlabel('latitude')
        set(gcf,'paperposition',[0 0 20 10])
        print(gcf,'-dpng','-r300',['DifLatRegress',FluxNames{i},'Over',FluxNames{j},'.png']);
        hold off;
    end
end




NormNAM = Indices.NAM/norm(Indices.NAM);
NormSAM = Indices.SAM/norm(Indices.SAM);
NormNINO34 = Indices.NINO34/norm(Indices.NINO34);
NormNAO = Indices.NAO/norm(Indices.NAO);

NAMNINODif = NormNAM-NormNINO34;
SAMNINODif = NormSAM-NormNINO34;
NAONINODif =NormNAO-NormNINO34;

NormalizedDif.NAMMinusNINO34 = NAMNINODif;
NormalizedDif.SAMMinusNINO34 = SAMNINODif;
NormalizedDif.NAOMinusNINO34 = NAONINODif;
NormalizedDif.NAMPlusNINO34= NormNAM+NormNINO34;
NormalizedDif.SAMPlusNINO34 = NormSAM+NormNINO34;
NormalizedIndexNames = fieldnames(NormalizedDif);

for j=1:length(NormalizedIndexNames)
    plot(NormalizedDif.(NormalizedIndexNames{j})','color',cmap(j*floor(64/ceil(length(NormalizedIndexNames))),:),'LineWidth',3);
    hold on;
    grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
end
set(gca,'xtick',11:12:size(MonthlyFluxDepartures.LW,3))
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
legend(NormalizedIndexNames)
set(gca,'xtick',11:12:size(MonthlyFluxDepartures.LW,3))
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year');
ylabel('Index');

% CalculateEOFs(MonthlyFluxDepartures.Temp,'Temp',1);
% CalculateEOFs(MonthlyFluxDepartures.Precip,'Precip',1);
% CalculateEOFs(MonthlyFluxDepartures.LWCF,'LWCF',1);
%%%%%%%%%%%%%%%%EOF CALCULATIONS%%%%%%%%%%%%%%%%%%
for i=1:length(FluxNames)
   [U.(FluxNames{i}),S.(FluxNames{i}),V.(FluxNames{i})] = CalculateEOFs(MonthlyFluxDepartures.(FluxNames{i}),FluxNames{i},1);
end
V.Net=reshape(V.Net,180,360,6);

NetTest = reshape(MonthlyFluxDepartures.Net,64800,156);
NetTestMeanSubtracted = NetTest - repmat(mean(NetTest),64800,1);
NetCovariance = NetTestMeanSubtracted'*NetTestMeanSubtracted;
trace(NetCovariance)
S.Net(1,1)^2/trace(NetCovariance) %fraction of variance held in first mode= 0.025... very very low

% FluxReshaped = structfun(@(x) reshape(x),64800,156, MonthlyFluxDepartures,'UniformOutput',false);
% Inputs to STRUCTFUN must be scalar structures.

for i=1:length(FluxNames)
   V.(FluxNames{i}) =  reshape(V.(FluxNames{i}),180,360,6);
end

save('USV.mat','U','S','V')
load USV.mat

corr(U.Net(:,1),Indices.NINO34')
plot(max(S.Net))
Svals = max(S.Net);
Svals(1)/sum(Svals); %this is much less than 20% of variance
TotalNetEOFRegression = zeros(180,360,156);
for i =1:6
OneEOF = TimeLatContourRegression(V.Net(:,:,i)*S.Net(i,i), U.Net(:,i)','Net',['EOF',num2str(i),' of Net']);
%why multiply V and S??
HemisphericFluxesClimSubtracted(OneEOF,0,90,['Net EOF',num2str(i)]);
TotalNetEOFRegression = TotalNetEOFRegression + OneEOF;
% TimeLatContourRegression(V.Net(:,:,2), U.Net(:,2)','Net','EOF2 of Net')
% TimeLatContourRegression(V.Net(:,:,3), U.Net(:,3)','Net','EOF3 of Net')
% TimeLatContourRegression(V.Net(:,:,4), U.Net(:,4)','Net','EOF4 of Net')
% TimeLatContourRegression(V.Net(:,:,5), U.Net(:,5)','Net','EOF5 of Net')
end

%how can the correlation coefficient be 1???



GlobalRegressionMap(MonthlyFluxDepartures.Net, V.Net,'Net','1st  EOFs');

load LatitudinalDifFluxTimeSeries.mat

fid = fopen('EOFs.txt','w');
for h=1:length(FluxNames)
for i=1:6
V_Temp = V.(FluxNames{h});
U_Temp = U.(FluxNames{h});
TempEOF = sum(V_Temp(:,:,i),2);
    EOF_NH(i) = sum(TempEOF(91:end));
    EOF_SH(i) = sum(TempEOF(1:90));
    EOF_HemDif(i) = EOF_NH(i)-EOF_SH(i);
    EOF_Global(i) = sum(TempEOF);
    EOF_NFlux_Corr(i)=corr(U_Temp(:,i),FluxTimeSeries.Net.AllLats.N');
    EOF_SFlux_Corr(i)=corr(U_Temp(:,i),FluxTimeSeries.Net.AllLats.S');
    EOF_DifFlux_Corr(i)=corr(U_Temp(:,i),FluxTimeSeries.Net.AllLats.Dif');
    EOF_SumFlux_Corr(i)=corr(U_Temp(:,i),FluxTimeSeries.Net.AllLats.N'+FluxTimeSeries.Net.AllLats.S');
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' NH = ', num2str(EOF_NH(i))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' SH = ', num2str(EOF_SH(i))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' NH-SH Dif = ', num2str(EOF_HemDif(i))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' Global = ', num2str(EOF_Global(i))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' Correlation with NH Net Flux = ', num2str(EOF_NFlux_Corr(i))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' Correlation with SH Net Flux = ', num2str(EOF_SFlux_Corr(i))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' Correlation with NH-SH Net Flux = ', num2str(EOF_DifFlux_Corr(i))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' Correlation with NH+SH Net Flux = ', num2str(EOF_SumFlux_Corr(i))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' Correlation with NAM = ', num2str(corr(U_Temp(:,i),Indices.NAM'))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' Correlation with NINO34 = ', num2str(corr(U_Temp(:,i),Indices.NINO34'))]);
    fprintf(fid, '%s\n',[FluxNames{h}, ' EOF',num2str(i), ' Correlation with SAM = ', num2str(corr(U_Temp(:,i),Indices.SAM'))]);
    fprintf(fid, '========\n');
end
fprintf(fid, '**************************\n');
end
fclose(fid);

% 
% for j=1:length(FluxNames)
% for i=1:6
%     TempFlux= V.(FluxNames{j});
%     TempEOF = sum(TempFlux(:,:,i),2);
%     EOF_NH.(FluxNames{j})(i) = sum(TempEOF(91:end));
%     EOF_SH.(FluxNames{j})(i) = sum(TempEOF(1:90));
%     EOF_HemDif.(FluxNames{j})(i) = EOF_NH.(FluxNames{j})(i)-EOF_SH.(FluxNames{j})(i);
%     EOF_Global.(FluxNames{j})(i) = sum(TempEOF);
%     EOF_NFlux_Corr.(FluxNames{j})(i)=corr(U.Net(:,i),FluxTimeSeries.Net.AllLats.N');
%     EOF_SFlux_Corr.(FluxNames{j})(i)=corr(U.Net(:,i),FluxTimeSeries.Net.AllLats.S');
%     EOF_DifFlux_Corr.(FluxNames{j})(i)=corr(U.Net(:,i),FluxTimeSeries.Net.AllLats.Dif');
%     EOF_SumFlux_Corr.(FluxNames{j})(i)=corr(U.Net(:,i),FluxTimeSeries.Net.AllLats.N'+FluxTimeSeries.Net.AllLats.S');
% end
% end

structfun(@(x) (max(sum(x))/sum(sum(x))),S,'UniformOutput',false) %none of the leading EOF explain morethan 28% of the var. most likely much lower.
%LW leading mode seems to have the best.. but still not enough..


TimeLatContour(TotalNetEOFRegression,'Net EOF TimeSeries for 1st 6 EOFs');
TimeLatContour(MonthlyFluxDepartures.Net,'Net');

plot(U.Net)
legend('EOF1','EOF2','EOF3','EOF4','EOF5','EOF6')
grid on;
%  repmat({'hi'},1,5)+int2str(1:5)
HemisphericFluxesClimSubtracted(TotalNetEOFRegression,0,90,'Net EOF TimeSeries for 1st 6 EOFs');
%should look at what % of total variance that all the EOFs explain!

NINO34Normalized = Indices.NINO34/norm(Indices.NINO34);

% load LatWeights.mat
% FluxWeighted=bsxfun(@times,MonthlyFluxDepartures.LW,sqrt(cosd(LatWeights(:,1))));
% %FluxWeighted=MonthlyFluxDepartures.LW;
% LW = reshape(FluxWeighted,180*360,156)';
% %clear MonthlyFluxDepartures;
% [U,S,V] = svds(LW);
% Regress1DTimeSeriesMap(MonthlyFluxDepartures.LW, U(:,1),'LW','EOF1 of LW')
% Regress1DTimeSeriesMap(MonthlyFluxDepartures.LW, U(:,2),'LW','EOF2 of LW')
% RescaledV = reshape(V,180,360,6); %these are the 6 principal components
% for i=1:6
%     subplot(3,3,i)
%     title(['PC ', i])
%     GlobePlot(RescaledV(:,:,i),0);
% end
% GlobePlot(RescaledV(:,:,1),0)
% GlobePlot(RescaledV(:,:,2),0)
% GlobePlot(RescaledV(:,:,3),0)
% GlobePlot(RescaledV(:,:,4),0)
% GlobePlot(RescaledV(:,:,5),0)
% GlobePlot(RescaledV(:,:,6),0)
% Regress1DTimeSeriesMap(MonthlyFluxDepartures.LW,NAMNINODif)
% Regress1DTimeSeriesMap(MonthlyFluxDepartures.LW,SAMNINODif)

cmap=jet;
for i=1:length(FluxNames)
    LatMeans = squeeze(mean(MonthlyFluxDepartures.(FluxNames{i}),2));
    for j=1:length(NormalizedIndexNames)
        plot(corr(LatMeans', NormalizedDif.(NormalizedIndexNames{j})'),'color',cmap(j*floor(64/ceil(length(NormalizedIndexNames))),:),'LineWidth',3);
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
    legend(NormalizedIndexNames)
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['AllLatIndexCorrs_',FluxNames{i},'.png']);
    hold off;
end

for i=1:length(FluxNames)
    for j=length(NormalizedIndexNames)-1:length(NormalizedIndexNames)
        Regress1DTimeSeriesMap(MonthlyFluxDepartures.(FluxNames{i}),NormalizedDif.(NormalizedIndexNames{j}),FluxNames{i},NormalizedIndexNames{j});
    end
end

for i=1:length(FluxNames)
    LatitudeTimeFlux = squeeze(mean(MonthlyFluxDepartures.(FluxNames{i}),2));
    contourf(LatitudeTimeFlux,12);colorbar
    Factor = 180/size(LatitudeTimeFlux,1);
    grid on;
    caxis([-max(max(abs(LatitudeTimeFlux))) max(max(abs(LatitudeTimeFlux)))])
    set(gca,'FontSize',20)
    xlabel('Time')
    ylabel('Latitude')
    set(gca,'xtick',11:12:size(LatitudeTimeFlux,2))
    set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
    set(gca,'ytick',10/Factor:10/Factor:180/Factor)
    set(gca,'YTickLabel',-80:10:90)
    title(['Monthly Deviation of ',FluxNames{i},' across all latitudes(Watts)'])
    set(gca,'GridLineStyle','--')
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['Monthly DeviationLatitudeTime',FluxNames{i},'.png']);
    hold off;
end

X = [NormalizedDif.(NormalizedIndexNames{1})' NormalizedDif.(NormalizedIndexNames{2})' NormalizedDif.(NormalizedIndexNames{3})'];

for i=1:length(FluxNames)
    Beta = inv(X'*X)*X'*squeeze(mean(MonthlyFluxDepartures.(FluxNames{i}),2))';   
    plot(Beta','LineWidth',3)
    grid on;
    IndexLegend = NormalizedIndexNames(1:3);
    legend(IndexLegend,'Location','NorthWestOutside')
    set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{i}),2),1))/180:180.5)
    set(gca,'xticklabel',num2cell(-90:10:90))
    set(gca,'FontSize',18)
    title(['Latitudinal Regression Coefficient of each Index for ', FluxNames{i}])
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['IndexLatRegressions_',FluxNames{i},'.png']);
    hold off;
end

contourf(corr([squeeze(mean(MonthlyFluxDepartures.Net,2))'],[squeeze(mean(MonthlyFluxDepartures.Precip,2))']))

for i=1:length(FluxNames)
    for j=i+1:length(FluxNames)
        GlobalLatCorrelationPlot(MonthlyFluxDepartures.(FluxNames{i}),MonthlyFluxDepartures.(FluxNames{j}),FluxNames{i},FluxNames{j});
    end
end

FluxNames = fieldnames(MonthlyFluxDepartures);
for j=1:length(IndexNames)
    for i=1:length(FluxNames)
        LatMeans = squeeze(mean(MonthlyFluxDepartures.(FluxNames{i}),2));
        x = 1:180/size(mean(MonthlyFluxDepartures.(FluxNames{i}),2),1):180;
        xq = 1:180/size(mean(MonthlyFluxDepartures.(FluxNames{1}),2),1):180;
        plot(interp1(x,corr(Indices.(IndexNames{j})', LatMeans'),xq),'color',cmap(i*floor(64/ceil(length(FluxNames))),:),'LineWidth',3);
        grid on;
        set(gca,'GridLineStyle','--')
        set(gca,'FontSize',18)
        set(gca,'xtick',0.5:10*(size(mean(MonthlyFluxDepartures.(FluxNames{1}),2),1))/180:180.5)
        set(gca,'xticklabel',num2cell(-90:10:90))
        xlabel('latitude')
        ylabel('correlation')
        title(['Correlation between ',IndexNames{j},' and monthly departure of flux at each latitude'])
        hold on;
    end
    legend(FluxNames,'Location','BestOutside')
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['AllLatIndexCorrs_',IndexNames{j},'.png']);
    hold off;
end


