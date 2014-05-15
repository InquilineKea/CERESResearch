function [difNHSH,NH,SH,MWA,MonthlySubtractedFluxes] = NHMinusSH(Flux,LowerLat,HigherLat,FluxName)

inputname(1)

load LatWeights.mat
lat = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');
time = size(Flux,3);
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

%LatWeights is basically cos(lat)
for j=1:time
    latFluxes = bsxfun(@rdivide,mean(Flux(:,:,j),2).*LatWeights(:,2),sum(LatWeightsAlone(LatWeights(:,1) < HigherLat & LatWeights(:,1) > LowerLat)));
    %mean(permuteNet(:,:,1:12:end),3) %march means...
    %for contribution at each latitude.. shouldn't it be the sum of
    %everything rather than the mean of everything? but then mean makes the
    %units in the thousands..
    difNHSH(j) = sum(latFluxes(lat > LowerLat & lat < HigherLat)-latFluxes(lat < -LowerLat & lat > -HigherLat));
    NH(j) = sum(latFluxes(lat > LowerLat & lat < HigherLat)); %sum across latitude lines
    SH(j) = sum(latFluxes(lat < -LowerLat & lat > -HigherLat));
end

for j=1:time
    Indices = mod(j,12):12:time;
    if mod(j,12) == 0
        Indices = 12:12:time;
    end
    MonthlySubtractedFluxes(1,j) = difNHSH(j) - mean(difNHSH(Indices));
    MonthlySubtractedFluxes(2,j) = NH(j) - mean(difNHSH(Indices));
    MonthlySubtractedFluxes(3,j) = SH(j) - mean(difNHSH(Indices));     
end %only for calculating monthly subtracted fluxes, which stops mattering...

% 
% %%%%%%%%%%%%no climatology subtraction%%%%%%%%%%%%%%%%%
set(gca,'FontSize',20)
windowSize = 12;
AllRunningMeans = filter([31/365 28/365 31/365 30/365 31/365 30/365 31/365 31/365 30/365 31/365 30/365 31/365],1,difNHSH); %but wouldn't this make jan count as april for some starts?
plot(AllRunningMeans(12:end))
grid on;
set(gca,'xtick',[10:12:time])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year where mean ends on');
ylabel('NH-SH Difference in Total Heat Flux (rsdt - rsutcs - rlutcs)');

days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
allMonthWeights = repmat(days_per_month,1,12);
allMonthWeights = [31 30 31 30 31 31 30 31 30 31 allMonthWeights 31 28];
allMonthWeights(48) = 29; %leap year, http://www.wolframalpha.com/input/?i=months+between+march+2000+and+february+2004
allMonthWeights(96) = 29;

%1 is march 2000. leap years in february 2004, 2008

moving_sum = @(n, x) filter(ones(1,n), 1, x);
moving_weighted_avg = moving_sum(12, difNHSH .* allMonthWeights) ...
    ./ moving_sum(12, allMonthWeights);
NHmoving_weighted_avg = moving_sum(12, NH .* allMonthWeights) ...
    ./ moving_sum(12, allMonthWeights);
SHmoving_weighted_avg = moving_sum(12, SH .* allMonthWeights) ...
    ./ moving_sum(12, allMonthWeights);

CorrectedNHmoving_weighted_avg = NHmoving_weighted_avg(12:end) - mean(NHmoving_weighted_avg(12:end));
CorrectedSHmoving_weighted_avg = SHmoving_weighted_avg(12:end) - mean(SHmoving_weighted_avg(12:end));

CorrectedGlobalMWA = CorrectedNHmoving_weighted_avg+CorrectedSHmoving_weighted_avg ;%-mean(CorrectedNHmoving_weighted_avg'+CorrectedSHmoving_weighted_avg')

Corrected_HemDif_MWA = moving_weighted_avg(12:end); 
%index 1 is march 2000 to february 2001
%index 13 is march 2001 to february 2002
%index 11 is january 2001 to december 2001

[AX,H1,H2] = plotyy(1:length(CorrectedNHmoving_weighted_avg),[CorrectedNHmoving_weighted_avg' ...
    CorrectedSHmoving_weighted_avg' CorrectedGlobalMWA'],1:length(CorrectedNHmoving_weighted_avg),Corrected_HemDif_MWA);
grid on;
set(gca,'xtick',[11:12:119+26])
set(AX(2),'XTickLabel',[])
set(gca,'XTickLabel',{2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year end');
set(get(AX(1),'Ylabel'),'FontSize',20,'String','Heat Flux in Latitude Band for NH, SH, Global after Mean-Subtraction (Watts/m^2)') 
set(get(AX(2),'Ylabel'),'FontSize',20,'String','NH-SH Difference in Total Heat Flux (Watts/m^2)','FontSize',20) 
difSD = sqrt(var(CorrectedNHmoving_weighted_avg) + var(CorrectedSHmoving_weighted_avg) - 2*getfield(cov(CorrectedNHmoving_weighted_avg,CorrectedSHmoving_weighted_avg), {1,2}));
GlobalSD = sqrt(var(CorrectedNHmoving_weighted_avg) + var(CorrectedSHmoving_weighted_avg) + 2*getfield(cov(CorrectedNHmoving_weighted_avg,CorrectedSHmoving_weighted_avg), {1,2}));
corrcof = getfield(corrcoef(CorrectedNHmoving_weighted_avg,CorrectedSHmoving_weighted_avg),{1,2});
set(legend(['NH SD = ', num2str(std(CorrectedNHmoving_weighted_avg))],['SH SD = ', num2str(std(CorrectedSHmoving_weighted_avg))] , ['Global SD = ', num2str(GlobalSD)],sprintf(['NH-SH SD = ', num2str(difSD), '\n Corr Cof = ', num2str(corrcof)])),'Location','BestOutside')
set(H2,'linewidth',3)% to change the first line
set(H1,'linewidth',4)
set(H2,'LineStyle','--')
set(AX,'FontSize',20)
set(gcf, 'Units','inches', 'Position',[0 0 20 10])
set(gca, 'Units','inches', 'Position',[1 1 16 8])
title([num2str(LowerLat),'-',num2str(HigherLat),' degree Latitude MWA in ', FluxName])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_',FluxName,'.png']);

NH_MWA = CorrectedNHmoving_weighted_avg;
SH_MWA = CorrectedSHmoving_weighted_avg;

blahblah = xcorr(NH_MWA,SH_MWA);

MWA = [Corrected_HemDif_MWA;NH_MWA;SH_MWA];

DifMWA = moving_weighted_avg;
corrcoef(CorrectedNHmoving_weighted_avg,CorrectedSHmoving_weighted_avg)
corrcoef(CorrectedNHmoving_weighted_avg,CorrectedGlobalMWA)
corrcoef(CorrectedNHmoving_weighted_avg,DifMWA(12:end))
corrcoef(CorrectedSHmoving_weighted_avg,CorrectedGlobalMWA)
corrcoef(CorrectedSHmoving_weighted_avg,DifMWA(12:end))
corrcoef(DifMWA(12:end),CorrectedGlobalMWA)

FluxFile = [FluxName, '_MWA_', num2str(LowerLat), '_', num2str(HigherLat),'.txt'];

fid = fopen(FluxFile,'w');
fprintf(fid, '%s\n',[FluxName, ' NH SD = ', num2str(  std(CorrectedNHmoving_weighted_avg) )]);
fprintf(fid, '%s\n',[FluxName, ' SH SD = ', num2str( std(CorrectedSHmoving_weighted_avg) )]);
fprintf(fid, '%s\n',[FluxName, ' NH+SH (Global) SD = ', num2str( GlobalSD )]);
fprintf(fid, '%s\n',[FluxName, ' NH-SH SD  = ', num2str(  difSD )]);

fprintf(fid, '%s\n',[FluxName, ' NH vs. SH Corr = ', num2str(  getfield(corrcoef(CorrectedNHmoving_weighted_avg,CorrectedSHmoving_weighted_avg)  ,{1,2}))]);
fprintf(fid, '%s\n',[FluxName, ' NH vs. Global Corr = ', num2str(  getfield(corrcoef(CorrectedNHmoving_weighted_avg,CorrectedGlobalMWA),{1,2}))]);
fprintf(fid, '%s\n',[FluxName, ' NH vs. NH-SH Dif Corr = ', num2str(  getfield(corrcoef(CorrectedNHmoving_weighted_avg,DifMWA(12:end))  ,{1,2}))]);
fprintf(fid, '%s\n',[FluxName, ' SH vs. Global Corr = ', num2str(  getfield(corrcoef(CorrectedSHmoving_weighted_avg,CorrectedGlobalMWA)  ,{1,2}))]);
fprintf(fid, '%s\n',[FluxName, ' SH vs. NH-SH Dif Corr = ', num2str(  getfield(corrcoef(CorrectedSHmoving_weighted_avg,DifMWA(12:end))  ,{1,2}))]);
fprintf(fid, '%s\n',[FluxName, ' Global vs. NH-SH Dif Corr = ', num2str(  getfield(corrcoef(DifMWA(12:end),CorrectedGlobalMWA) ,{1,2}))]);
    fprintf(fid, '========\n');
fprintf(fid, '**************************\n');
fclose(fid);



end