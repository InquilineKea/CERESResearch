function [latFluxes] = LatWeightedAverage(Flux,Year,titleName)

load LatWeights.mat

colorMarker = {'gx','bx','ro','yd','c^','m-','rs','k<','ch','ko','kx','k-'};
cmap = jet;

lat = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');
lon = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lon');
time = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time');
% 180*sind(lat)

%Monthly averaged, EBAF product, 3/2000-12/2010, all at TOA: downward shortwave (rsdt), upward shortwave (rsut), upward shortwave in clear-sky (rsutcs), upward longwave (rlut), upward longwave in clear-sky (rlutcs)

%CERES Energy Balanced and Filled (EBAF) Ed2.6r data product.

[LON,LAT]=meshgrid(lon,lat);

MonthIndices = 11+12*(Year-2001):22+12*(Year-2001);

difNHSH = zeros(1,length(MonthIndices));
set(gca,'FontSize',20)
for j=MonthIndices
   %plot(LAT(:,2),sum(permuteNet(:,:,j),2).*LatWeights(:,2),colorMarker{mod(j+1:j+1,12)+1})
   %plot(LAT(:,2),sum(permuteNet(:,:,j),2).*LatWeights(:,2),'color',cmap((mod(j+1:j+1,12)+1)*5,:),'LineWidth',3)
   plot(sind(LAT(:,2)),mean(Flux(:,:,j),2),'color',cmap((mod(j+1:j+1,12)+1)*5,:),'LineWidth',3)
   latFluxes = mean(Flux(:,:,j),2).*LatWeights(:,2);
   difNHSH(j) = sum(latFluxes(91:end)-latFluxes(1:90));
   %plot(sum(permuteNet(91:end,:,j),2) - sum(flipud(permuteNet(1:90,:,j)),2))
   hold on;   
end
grid on;
title([num2str(Year),' ',titleName]);
xlabel('Latitude');
ylabel('Average Latitudinal Flux (Watts)');
set(gca,'xtick',sind(-90:30:90))
set(gca,'xticklabel',{'90S','60S','30S','0','30N','60N','90N'})
set(legend('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'),'Location','BestOutside') %NorthWestOutside

set(gcf,'paperposition',[1 1 24 12])
print(gcf,'-dpng','-r300',[titleName,'_','fig.png']);

end