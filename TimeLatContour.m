function Flux1 = TimeLatContour(Flux1,FluxName)    
load LatWeights
%Flux1=bsxfun(@times,Flux1, TimeSeries);

    %then lat-time contour plot it
    [X,Y]=meshgrid(1:156,sind(LatWeights(:,1)));
    contourf(X,Y,squeeze(mean(Flux1,2)));colorbar
%     figure(2)
%     contourf(squeeze(mean(Flux1,2)));
%     blah=squeeze(mean(Flux1,2));
    caxis([-max(max(abs(squeeze(mean(Flux1,2))))) max(max(abs(squeeze(mean(Flux1,2)))))])
    caxis([-8 8])
%     MakeLizMap;
% colormap(lizmap)

    grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
    set(gca,'xtick',11:12:size(Flux1,3))
    set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
    xlabel('Year');
%     plot(sind(LatWeights(:,1)))
%     plot(-1:1)
    set(gca,'ytick',sind((0.5:10*180/180:179.5)-90))
    set(gca,'yticklabel',num2cell(-90:10:90))
    ylabel('latitude')
    title(['Time-Latitudinal Contour Plot of ', FluxName])
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',['TimeLatRegressContours',FluxName,'.png']);
    hold off;
end