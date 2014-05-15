function Matrix = PlotIndexThreshold(Flux1,TimeSeries,Percentile,Var1,Var2)

TimePoints = find(TimeSeries > prctile(TimeSeries,Percentile));
ensoTime = ncread('enso.cdf','T');

NINO34 = ncread('enso.cdf','NINO34');
NINO34 = NINO34(find(ensoTime > 482 & ensoTime < 638));

% TimePoints = find(NINO34 < prctile(NINO34,90) & NINO34 > prctile(NINO34,10) & TimeSeries > prctile(TimeSeries,Percentile));


lats = size(Flux1,1);
longs = size(Flux1,2);
time=size(Flux1,3);
Factor = 180/lats;

Flux2Plot = mean(Flux1(:,:,TimePoints),3);
load LatWeights.mat

FluxMean = mean(mean(resizem(Flux2Plot,Factor),2).*LatWeights(:,2));

load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(resizem(Flux2Plot,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
caxis([-max(max(abs(Flux2Plot))) max(max(abs(Flux2Plot)))])
title(['Avg monthly ',Var1,' departure when index of ',Var2,' exceeds ', num2str(Percentile),' Percentile. Mean = ', num2str(FluxMean)])
grid on;
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Exceeds',num2str(Percentile),'ThresholdMonDepar_',Var1,'-',Var2,'.png']);
hold off;


end