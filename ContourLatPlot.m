function [blah] = ContourLatPlot(Flux, FluxName)
load LatWeights.mat

time = size(Flux,3);
    [X,Y]=meshgrid(1:time,sind(LatWeights(:,1)));
% contourf(X,Y,squeeze(mean(bsxfun(@times,Flux,LatWeights(:,2)),2)));colorbar
contourf(X,Y,squeeze(mean(Flux,2)),20);colorbar
grid on;
set(gca,'FontSize',20)
xlabel('Time')
ylabel('Latitude')
set(gca,'xtick',11:12:time)
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
    set(gca,'ytick',sind((0.5:10*(size(mean(Flux,2),1))/180:179.5)-90))
    set(gca,'yticklabel',num2cell(-90:10:90))
title(['Latitudinal Monthly Anomaly of Climatology-Subtracted-',FluxName,' (Watts/m^2)'])
% caxis([-max(max(abs(squeeze(mean(bsxfun(@times,Flux,LatWeights(:,2)),2))))) max(max(abs(squeeze(mean(bsxfun(@times,Flux,LatWeights(:,2)),2)))))])
caxis([-max(max(abs(squeeze(mean(Flux,2))))) max(max(abs(squeeze(mean(Flux,2)))))])
    MakeLizMap
    colormap(lizmap)
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['ContourWeighted',FluxName,'Anomaly.png']);

set(gca,'FontSize',20)
plot(sind(LatWeights(:,1)),std(squeeze(mean(Flux,2)),0,2),'LineWidth',3)
%std of latitudinal averaged flux over time. this is prolly preferable
grid on;
% set(gca,'xtick',(0.5-90:10:179.5-90))
set(gca,'xtick',sind((0.5:10*180/180:179.5)-90))
set(gca,'xticklabel',num2cell(-90:10:80))
title(['Temporal SD of ',FluxName, ' Monthly Anomaly'])
xlabel('Latitude')
ylabel('Temporal SD')
view(90,-90)
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[FluxName,'Lat_Anomaly_SD','.png']);

LatAnomaly = squeeze(mean(Flux,2));
bb = corr(flipud(LatAnomaly(1:end/2,:))',LatAnomaly(end/2+1:end,:)')
set(gca,'FontSize',20)
plot(diag(bb),'LineWidth',3)
grid on;
xlabel('Latitude (degrees away from equator)')
ylabel('Correlation Coefficient')
title(['Correlation between between Monthly Anomalies at X Latitude N and X Latitude S in ', FluxName])
view(90,-90)
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[FluxName,'Lat_Anomaly_North_vs_South_Correlation','.png']);


end