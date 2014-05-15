function RegressedFlux = MWATimeLatContourRegression(Alpha,Beta,MWA,TimeSeries,FluxName,IndexName,MonthFilterSize)    
load LatWeights
lats = size(Beta,1);
longs = size(Beta,2);
time=length(TimeSeries);
Factor = 180/lats;

%RegressedFlux=bsxfun(@times,Flux1, TimeSeries);
RegressedFlux=bsxfun(@times,Beta, permute(TimeSeries',[3 2 1])) + repmat(Alpha,[1,1,length(TimeSeries)]);
TimeMean =  repmat(mean(MWA,3),[1,1,length(TimeSeries)]) ;
SSTotal = sum((MWA-TimeMean ).^2,3);
SSRes = sum((MWA-RegressedFlux).^2,3);%%why is this so high? Why is SS_Res >> SS_Total?
NumMonthsAfterMAFiltering = time-(MonthFilterSize-1);

%Flux1 is the moving average of the flux... 
%TimeSeries is the MA of the index..
%Regression of a moving average on a MA... 

%so why is SSRes so high? isn't that supposed to be mathematically
%impossible? becuase you need to input the regression coefficients into RegressedFlux...

% R_Squared = 1 - SSRes./SSTotal;
  
MakeLizMap;
colormap(lizmap)
% load geoid
% set(gca,'FontSize',20)
% ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
% geoshow(resizem(R_Squared,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% 
%     %then lat-time contour plot it
%         caxis([-1 1])
%         title([num2str(MonthFilterSize), ' Month Moving Avg. R^2 of ', FluxName,' over ',IndexName])
%     set(gcf,'paperposition',[0 0 20 10])
%     print(gcf,'-dpng','-r300',[FluxName, num2str(MonthFilterSize),'_Month_MA_R_Squared','_',IndexName,'.png']);

    hold off;

    [X,Y]=meshgrid(1:length(TimeSeries),sind(LatWeights(:,1)));
    contourf(X,Y,squeeze(mean(RegressedFlux,2)),30);colorbar
    shading flat
    caxis([-max(max(abs(squeeze(mean(RegressedFlux,2))))) max(max(abs(squeeze(mean(RegressedFlux,2)))))])

% lizmap(165:195,:) = 1;
 colormap(lizmap)
    grid on; set(gca,'FontSize',18);set(gca,'GridLineStyle','--')
set(gca,'xtick',[12-2-(MonthFilterSize-1):12:length(TimeSeries)]) %what if Filter Size is 12 months?
% [12-2-(MonthFilterSize-1):12:NumMonthsAfterMAFiltering]
set(gca,'XTickLabel',2000:2012)
xlabel('Year end');
%     plot(sind(LatWeights(:,1)))
%     plot(-1:1)
    set(gca,'ytick',sind((0.5:10*180/180:179.5)-90))
    set(gca,'yticklabel',num2cell(-90:10:90))
    ylabel('latitude')
    title([FluxName,' over ',IndexName, ' ', num2str(MonthFilterSize),'-Month Moving Average Regression'])
    set(gcf,'paperposition',[0 0 20 10])
    print(gcf,'-dpng','-r300',[FluxName,'_',IndexName,' ',num2str(MonthFilterSize),'_Month_MvgAvgTimeLat','.png']);
    hold off;
end
% 
% geoshow(resizem(mean(MWA.Net,3),Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% geoshow(resizem(mean(RegressedFlux,3),Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% geoshow(resizem(Beta,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% 
% 
% geoshow(resizem(blah(:,:,1),Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% 
% 
%   geoshow(resizem(SSRes./SSTotal,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
%   geoshow(resizem(SSTotal,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% 
%     geoshow(resizem(SSTotal-SSRes,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
