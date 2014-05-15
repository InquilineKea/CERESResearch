function GlobalCorrTimeSeries = Corr1DTimeSeriesMap(Flux1,TimeSeries,VarName1,VarName2)
lats = size(Flux1,1);
longs = size(Flux1,2);
time=size(Flux1,3);
Factor = 180/lats;
%time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));

%TimeSeries = NetDifMonthTotal;
%Flux1 = NetClimatologySubtracted;
TimeSeries3D = ones(1,1,time);
TimeSeries3D(1,1,:) = TimeSeries;
TimeSeries3D = repmat(TimeSeries3D,[180/Factor 360/Factor 1]);

GlobalCorrs = bsxfun(@rdivide,sum(Flux1.*TimeSeries3D,3),(sqrt(sum(Flux1.^2,3)).*sqrt(sum(TimeSeries3D.^2,3))));
load LatWeights.mat
CorrMean = mean(mean(resizem(GlobalCorrs,Factor),2).*LatWeights(:,2));

load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(resizem(GlobalCorrs,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
caxis([-max(max(abs(GlobalCorrs))) max(max(abs(GlobalCorrs)))])

title(['Correlation between ',VarName1,' and Time Series of ',VarName2, '. Mean = ', num2str(CorrMean)])
grid on;
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['GlobalTimeSeriesCorr_',VarName1,'-',VarName2,'.png']);
hold off;

end