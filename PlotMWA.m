function [NH,SH,HemDif] = PlotMWA(Flux,LowerLat,HigherLat,VariableName)

%time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));
time = length(Flux(1,1,:));
load LatWeights.mat
lat = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');

difNHSH = zeros(1,time);
NH = zeros(1,time);
SH = zeros(1,time);
MonthlySubtractedFluxes = zeros(3,time);

%LAND ONLY
load('topo.mat','topo','topomap1');
%Flux = bsxfun(@times, Flux,topo>=0); 
%Flux = bsxfun(@times, Flux,topo<0); 
%also this reasoning can be applied if i just want to select land tiles,
%too..

LatWeightsAlone = LatWeights(:,2);

for j=1:time
    latFluxes = bsxfun(@rdivide,mean(Flux(:,:,j),2).*LatWeights(:,2),sum(LatWeightsAlone(LatWeights(:,1) <= HigherLat & LatWeights(:,1) > LowerLat)));
    %mean(permuteNet(:,:,1:12:end),3) %march means...
    difNHSH(j) = sum(latFluxes(lat > LowerLat & lat <= HigherLat)-latFluxes(lat < -LowerLat & lat >= -HigherLat));
    NH(j) = sum(latFluxes(lat > LowerLat & lat <= HigherLat)); %sum across latitude lines
    SH(j) = sum(latFluxes(lat < -LowerLat & lat >= -HigherLat));
end

HemDif = NH-SH;

[AX,H1,H2] = plotyy(1:length(NH),[NH' ...
    SH'],1:length(NH),HemDif);
grid on;
set(gca,'xtick',[11:12:119+26])
set(AX(2),'XTickLabel',[])
set(gca,'XTickLabel',{2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year end');
set(get(AX(1),'Ylabel'),'FontSize',20,'String','Heat Flux in Latitude Band for NH, SH (Watts/m^2)') 
set(get(AX(2),'Ylabel'),'FontSize',20,'String','NH-SH Difference in Total Heat Flux (Watts/m^2)','FontSize',20) 
difSD = sqrt(var(NH) + var(SH) - 2*getfield(cov(NH,SH), {1,2}));
corrcof = getfield(corrcoef(NH,SH),{1,2});
set(legend(['NH SD = ', num2str(std(NH))],['SH SD = ', num2str(std(SH))] ,sprintf(['NH-SH SD = ', num2str(difSD), '\n Corr Cof = ', num2str(corrcof)])),'Location','BestOutside')
set(H1,'linewidth',3)
set(H2,'linewidth',3)% to change the first line
set(AX,'FontSize',20)
set(gcf, 'Units','inches', 'Position',[0 0 20 10])
set(gca, 'Units','inches', 'Position',[1 1 16 8])
title([num2str(LowerLat),'-',num2str(HigherLat),'deg MWA ', VariableName])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[VariableName,'_MWA_',num2str(LowerLat),'-',num2str(HigherLat),'.png']);
%hold off;


end