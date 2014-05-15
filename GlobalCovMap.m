function GlobalCovs = GlobalCovMap(Flux1,Flux2)
inputname(1);
inputname(2);
time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));

GlobalCovs = bsxfun(@rdivide,sum(Flux1.*Flux2,3),(time-1));

load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(GlobalCovs,geoidrefvec,'DisplayType','texturemap');colorbar
title(['Covariance between ',inputname(1),' and ',inputname(2)])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['GlobalCovs_',num2str(inputname(1)),'-',num2str(inputname(2)),'.png']);
hold off;

end

