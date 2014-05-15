function RegressedFlux = TimeLatContourRegression(Flux1,TimeSeries,FluxName,IndexName)    
load LatWeights
%RegressedFlux=bsxfun(@times,Flux1, TimeSeries);
 RegressedFlux=bsxfun(@times,Flux1, permute(TimeSeries',[3 2 1]));
 

%  SSRes = 
%  SSTotal = (RegressedFlux - mean(Flux1,3)).^2;

 
    %then lat-time contour plot it
    [X,Y]=meshgrid(1:length(TimeSeries),sind(LatWeights(:,1)));
    contourf(X,Y,squeeze(mean(RegressedFlux,2)));colorbar
%     figure(2)
%     contourf(squeeze(mean(RegressedFlux,2)));
%     blah=squeeze(mean(RegressedFlux,2));
    caxis([-max(max(abs(squeeze(mean(RegressedFlux,2))))) max(max(abs(squeeze(mean(RegressedFlux,2)))))])
%     MakeLizMap;
% colormap(lizmap)
MakeLizMap
colormap(lizmap)

    grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
    set(gca,'xtick',11:12:length(TimeSeries))
    set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
    xlabel('Year');
%     plot(sind(LatWeights(:,1)))
%     plot(-1:1)
    set(gca,'ytick',sind((0.5:10*180/180:179.5)-90))
    set(gca,'yticklabel',num2cell(-90:10:90))
    ylabel('latitude')
    title(['Time-Latitudinal Contour Plot of ', FluxName,' over ',IndexName])
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['TimeLatRegressContours',FluxName,'_',IndexName,'.png']);
    hold off;
end