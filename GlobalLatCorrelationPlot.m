function all = GlobalLatCorrelationPlot(Flux1,Flux2,VarName1,VarName2)

LatMeans1 = squeeze(mean(Flux1,2));
LatMeans2 = squeeze(mean(Flux2,2));
% contour3(SWClimatologySubtracted)
CorMatrix = corr([LatMeans1'],[LatMeans2']);

% contour3(CorMatrix)
%[c,h]=contourf(CorMatrix);colorbar
contourf(CorMatrix);colorbar
grid on;
set(gca,'GridLineStyle','--')
set(gca,'FontSize',18)
set(gca,'xtick',0.5:10*(size(mean(Flux2,2),1))/180:180.5)
set(gca,'xticklabel',num2cell(-90:10:90))
set(gca,'ytick',0.5:10*(size(mean(Flux1,2),1))/180:180.5)
set(gca,'yticklabel',num2cell(-90:10:90))
caxis([-1 1])
title(['Correlation Coefficients between each Latitudinal Band in ', VarName1, ' and ', VarName2])
xlabel(['Latitude Band in ', VarName2])
ylabel(['Latitude Band in ', VarName1])
grid on;
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['CorrsAllLats_',VarName1,'_',VarName2,'.png']);
hold off;
end