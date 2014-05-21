function a =PlotLatitudeSD(FluxDeparture, MonthFilterSize, VariableName)
load LatWeights.mat


set(gca,'FontSize',20)
plot(sind(LatWeights(:,1)),std(squeeze(mean(FluxDeparture,2)),0,2),'LineWidth',3)
%std of latitudinal averaged flux over time. this is prolly preferable
grid on;
% set(gca,'xtick',(0.5-90:10:179.5-90))
set(gca,'xtick',sind((0.5:10*180/180:179.5)-90))
set(gca,'xticklabel',num2cell(-90:10:80))
title(['Temporal SD of ',num2str(MonthFilterSize),' Month MA in ',VariableName])
xlabel('Latitude')
ylabel('Temporal SD')
view(90,-90)
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
end