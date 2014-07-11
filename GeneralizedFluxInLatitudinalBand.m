function [NH,SH,HemDif,HemSum,AsymmetryIndex] = GeneralizedFluxInLatitudinalBand(Flux,LowerLat,HigherLat,VariableName,MonthFilterSize,ListOfLats)
hold off;
%time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));
time = length(Flux(1,1,:));
% load LatWeights.mat
 load IndicesMWA.mat
% ListOfLats= ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');

% difNHSH = zeros(1,time);
% NH = zeros(1,time);
% SH = zeros(1,time);

%LAND ONLY
% load('topo.mat','topo','topomap1');
%Flux = bsxfun(@times, Flux,topo>=0); 
%Flux = bsxfun(@times, Flux,topo<0); 
%also this reasoning can be applied if i just want to select land tiles,
%too..

LatWeightsAlone = cosd(ListOfLats)./sum(cosd(ListOfLats));

for j=1:time
    latFluxes = bsxfun(@rdivide,mean(Flux(:,:,j),2).*LatWeightsAlone,sum(LatWeightsAlone(ListOfLats <= HigherLat& ListOfLats > LowerLat)));
    %shoudln't this be sum rather than mean? when i do sum it gets into the
    %thousands... can't work..
    %mean(permuteNet(:,:,1:12:end),3) %march means...
    AllFluxes.difNHSH(j) = sum(latFluxes(ListOfLats > LowerLat& ListOfLats<= HigherLat)-latFluxes(ListOfLats< -LowerLat& ListOfLats>= -HigherLat));
    AllFluxes.NH(j) = sum(latFluxes(ListOfLats> LowerLat& ListOfLats<= HigherLat)); %sum across latitude lines
    AllFluxes.SH(j) = sum(latFluxes(ListOfLats< -LowerLat& ListOfLats>= -HigherLat));
end

for j=1:time
    latFluxesGlobal = bsxfun(@rdivide,mean(Flux(:,:,j),2).*LatWeightsAlone,sum(LatWeightsAlone(ListOfLats <= 90 & ListOfLats > 0)));
    %shoudln't this be sum rather than mean? when i do sum it gets into the
    %thousands... can't work..
    %mean(permuteNet(:,:,1:12:end),3) %march means...
    AllFluxes.difNHSHGlobal(j) = sum(latFluxesGlobal(ListOfLats> 0 & ListOfLats<= 90)-latFluxesGlobal(ListOfLats< 0 & ListOfLats>= -90));
    AllFluxes.NHGlobal(j) = sum(latFluxesGlobal(ListOfLats> 0 & ListOfLats<= 90)); %sum across latitude lines
    AllFluxes.SHGlobal(j) = sum(latFluxesGlobal(ListOfLats< 0 & ListOfLats>= -90));
end

AllFluxes.HemDif = AllFluxes.NH-AllFluxes.SH;
AllFluxes.HemSum = (AllFluxes.NH+AllFluxes.SH)/2; %does this average to zero??? By nature this will be ~2x the magnitude of NH and SH... 
AllFluxes.AsymmetryIndex = AllFluxes.HemDif./(AllFluxes.HemSum*2);

days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
allMonthWeights = repmat(days_per_month,1,13);
allMonthWeights = [31 30 31 30 31 31 30 31 30 31 allMonthWeights 31 28];
allMonthWeights(48) = 29; %leap year, http://www.wolframalpha.com/input/?i=months+between+march+2000+and+february+2004
allMonthWeights(96) = 29;
moving_sum = @(n, x) filter(ones(1,n), 1, x);

AllFluxNames = fieldnames(AllFluxes);
for i=1:length(AllFluxNames)
AllFluxes.(AllFluxNames{i}) = ...
    bsxfun(@rdivide,moving_sum(MonthFilterSize, AllFluxes.(AllFluxNames{i}) .* allMonthWeights),moving_sum(MonthFilterSize,allMonthWeights)); AllFluxes.(AllFluxNames{i}) = AllFluxes.(AllFluxNames{i})(MonthFilterSize:end);
Mean.(AllFluxNames{i}) = mean(AllFluxes.(AllFluxNames{i}));
end
% AsymmetryIndex= bsxfun(@rdivide,moving_sum(MonthFilterSize, AsymmetryIndex .* allMonthWeights),moving_sum(MonthFilterSize,allMonthWeights)); AsymmetryIndex = AsymmetryIndex(MonthFilterSize:end);
% NH= bsxfun(@rdivide,moving_sum(MonthFilterSize, NH .* allMonthWeights),moving_sum(MonthFilterSize,allMonthWeights)); NH = NH(MonthFilterSize:end);
% SH= bsxfun(@rdivide,moving_sum(MonthFilterSize, SH .* allMonthWeights),moving_sum(MonthFilterSize,allMonthWeights)); SH = SH(MonthFilterSize:end);
% HemDif= bsxfun(@rdivide,moving_sum(MonthFilterSize, HemDif .* allMonthWeights),moving_sum(MonthFilterSize,allMonthWeights)); HemDif = HemDif(MonthFilterSize:end);
% HemSum= bsxfun(@rdivide,moving_sum(MonthFilterSize, HemSum .* allMonthWeights),moving_sum(MonthFilterSize,allMonthWeights)); HemSum = HemSum(MonthFilterSize:end);

NHMean = mean(AllFluxes.NH);
SHMean = mean(AllFluxes.SH);
GlobalMean = mean(AllFluxes.HemSum);
HemDifMean = mean(AllFluxes.HemDif);
HemDif = AllFluxes.HemDif;
difNHSHGlobal = AllFluxes.difNHSHGlobal;
NHGlobal = AllFluxes.NHGlobal;
SHGlobal = AllFluxes.SHGlobal;

if ~strcmp(VariableName,'Net') || ~(LowerLat==0) || ~(HigherLat==90)
    NH = AllFluxes.NH - mean(AllFluxes.NH);
    SH = AllFluxes.SH - mean(AllFluxes.SH);
    HemSum = AllFluxes.HemSum - mean(AllFluxes.HemSum);
else
    NH = AllFluxes.NH;
    SH = AllFluxes.SH;
    HemSum = AllFluxes.HemSum;
end
% VariableName

AsymmetryIndex = AllFluxes.AsymmetryIndex;

if strcmp(VariableName,'Precip')
        [AX,H1,H2] = plotyy(1:length(NH),[NH' SH' HemSum'],1:length(NH),[HemDif'  6*AsymmetryIndex']);
else
        [AX,H1,H2] = plotyy(1:length(NH),[NH' SH' HemSum'],1:length(NH),HemDif);
end
grid on;
set(gca,'xtick',12-2-(MonthFilterSize-1)-6:12:time)
set(AX(2),'XTickLabel',[])
set(gca,'XTickLabel',2000:2013)


if strcmp(VariableName,'Precip')
set(get(AX(1),'Ylabel'),'FontSize',20,'String','Precip in Latitude Band for NH, SH (mm/day)') 
set(get(AX(2),'Ylabel'),'FontSize',20,'String','NH-SH Difference in Precip','FontSize',20) 
% 'fadad'
else
set(get(AX(1),'Ylabel'),'FontSize',20,'String','Heat Flux in Latitude Band for NH, SH (Watts/m^2)') 
set(get(AX(2),'Ylabel'),'FontSize',20,'String','NH-SH Difference in Total Heat Flux (Watts/m^2)','FontSize',20) 
% 'SHOULD WORK'
end
difSD = sqrt(var(NH) + var(SH) - 2*getfield(cov(NH,SH), {1,2}));
GlobalSD = sqrt(var(NH) + var(SH) +2*getfield(cov(NH,SH), {1,2}));
corrcof = getfield(corrcoef(NH,SH),{1,2});
if strcmp(VariableName,'Precip')
        set(legend(['NH \mu = ', num2str(NHMean)],['SH \mu = ', num2str(SHMean)] ,['Global \mu = ' num2str(GlobalMean)], ['AsymIndex'],['NH-SH \mu = ', num2str(HemDifMean)]),'Location','BestOutside')
else
    set(legend(['NH \mu = ', num2str(NHMean)],['SH \mu = ', num2str(SHMean)] ,['Global \mu = ' num2str(GlobalMean)],sprintf(['NH-SH \\mu = ', num2str(HemDifMean), '\n Corr = ', num2str(corrcof)])),'Location','BestOutside')    
end
set(H1,'linewidth',4)
set(H2,'LineStyle','--')
set(H2,'linewidth',3)% to change the first line
set(AX,'FontSize',20)
xlabel('Year End');
uistack(H1(3), 'top') 
set(gcf, 'Units','inches', 'Position',[0 0 20 10])
set(gca, 'Units','inches', 'Position',[1 1 16 8])
title([num2str(LowerLat),'-',num2str(HigherLat),'deg ', num2str(MonthFilterSize),' Month Moving Avg ', VariableName, '. NH/SH ',sprintf(['Corr = ', num2str(corrcof)])])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',[VariableName, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.png']);
saveas(gcf,[VariableName, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.fig'],'fig')
hold off;

FluxFile = [VariableName, num2str(MonthFilterSize), '_MonthMA_', num2str(LowerLat), '_', num2str(HigherLat),'.txt'];

fid = fopen(FluxFile,'w');
fprintf(fid, '%s\n',[VariableName, ' NH SD = ', num2str(  std(NH) )]);
fprintf(fid, '%s\n',[VariableName, ' SH SD = ', num2str( std(SH) )]);
fprintf(fid, '%s\n',[VariableName, ' NH+SH (Global) SD just in latitudinal band = ', num2str( GlobalSD)]);
fprintf(fid, '%s\n',[VariableName, ' NH-SH SD  = ', num2str(  difSD )]);

fprintf(fid, '%s\n',[VariableName, ' NH vs. SH Corr = ', num2str(  getfield(corrcoef(NH,SH)  ,{1,2})), '. R-Squared = ',num2str(RSquared(NH,SH))]);
fprintf(fid, '%s\n',[VariableName, ' NH vs. Global NH+SH Corr = ', num2str(  getfield(corrcoef(NH,NHGlobal+SHGlobal),{1,2})), '. R-Squared = ',num2str(RSquared(NH,NHGlobal+SHGlobal))]);
fprintf(fid, '%s\n',[VariableName, ' NH vs. Global NH Corr = ', num2str(  getfield(corrcoef(NH,NHGlobal),{1,2})), '. R-Squared = ',num2str(RSquared(NH,NHGlobal))]);
fprintf(fid, '%s\n',[VariableName, ' SH vs. Global SH Corr = ', num2str(  getfield(corrcoef(SH,SHGlobal),{1,2})), '. R-Squared = ',num2str(RSquared(SH,SHGlobal))]);
fprintf(fid, '%s\n',[VariableName, ' NH vs. NH-SH Dif just in latitudinal band Corr = ', num2str(  getfield(corrcoef(NH,HemDif)  ,{1,2})), '. R-Squared = ',num2str(RSquared(NH,HemDif))]);
fprintf(fid, '%s\n',[VariableName, ' NH vs. Global NH-SH Dif Corr = ', num2str(  getfield(corrcoef(NH,NHGlobal-SHGlobal)  ,{1,2})), '. R-Squared = ',num2str(RSquared(NH,NHGlobal-SHGlobal))]);
fprintf(fid, '%s\n',[VariableName, ' SH vs. Global NH+SH Corr = ', num2str(  getfield(corrcoef(SH,NHGlobal+SHGlobal)  ,{1,2})), '. R-Squared = ',num2str(RSquared(SH,NHGlobal+SHGlobal))]);
fprintf(fid, '%s\n',[VariableName, ' SH vs. Global NH-SH Dif Corr = ', num2str(  getfield(corrcoef(SH,NHGlobal-SHGlobal)  ,{1,2})), '. R-Squared = ',num2str(RSquared(SH,NHGlobal-SHGlobal))]);
fprintf(fid, '%s\n',[VariableName, ' Latitude NH+SH vs. NH-SH Global Dif Corr = ', num2str(  getfield(corrcoef(HemSum,NHGlobal-SHGlobal) ,{1,2}))]);
    fprintf(fid, '========\n');
    
IndexNames = fieldnames(IndicesMvgAvg.Window12Months);
    
for i=1:length(IndexNames)
fprintf(fid, '%s\n',[VariableName, ' NH vs. ',IndexNames{i},' Corr = ', num2str(  getfield(corrcoef(NH,IndicesMvgAvg.Window12Months.(IndexNames{i}))  ,{1,2})), '. R-Squared = ', num2str(RSquared(IndicesMvgAvg.Window12Months.(IndexNames{i}),NH))]);
fprintf(fid, '%s\n',[VariableName, ' SH vs. ',IndexNames{i},' Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.(IndexNames{i}),SH)  ,{1,2})), '. R-Squared = ', num2str(RSquared(IndicesMvgAvg.Window12Months.(IndexNames{i}),SH))]);    
fprintf(fid, '%s\n',[VariableName, ' NH+SH vs. ',IndexNames{i},' Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.(IndexNames{i}),SH+NH)  ,{1,2})), '. R-Squared = ', num2str(RSquared(IndicesMvgAvg.Window12Months.(IndexNames{i}),NH+SH))]);  
fprintf(fid, '%s\n',[VariableName, ' NH-SH vs. ',IndexNames{i},' Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.(IndexNames{i}),NH-SH)  ,{1,2})), '. R-Squared = ', num2str(RSquared(IndicesMvgAvg.Window12Months.(IndexNames{i}),NH-SH))]);

end
    
fprintf(fid, '**************************\n');
fclose(fid);

% fprintf(fid, '%s\n',[VariableName, ' NH vs. NAM Corr = ', num2str(  getfield(corrcoef(NH,IndicesMvgAvg.Window12Months.NAM)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' NH vs. SAM Corr = ', num2str(  getfield(corrcoef(NH,IndicesMvgAvg.Window12Months.SAM)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' NH vs. NINO34 Corr = ', num2str(  getfield(corrcoef(NH,IndicesMvgAvg.Window12Months.NINO34)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' SH vs. NAM Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.NAM,SH)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' SH vs. SAM Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.SAM,SH)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' SH vs. NINO34 Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.NINO34,SH)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' NH+SH vs. NINO34 Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.NINO34,SH+NH)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' NH-SH vs. NINO34 Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.NINO34,NH-SH)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' NH+SH vs. SAM Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.SAM,SH+NH)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' NH-SH vs. SAM Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.SAM,NH-SH)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' NH+SH vs. NAM Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.NAM,SH+NH)  ,{1,2}))]);
% fprintf(fid, '%s\n',[VariableName, ' NH-SH vs. NAM Corr = ', num2str(  getfield(corrcoef(IndicesMvgAvg.Window12Months.NAM,NH-SH)  ,{1,2}))]);



end