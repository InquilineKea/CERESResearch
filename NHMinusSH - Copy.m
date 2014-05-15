function [difNHSH,NH,SH,MWA] = NHMinusSH(permuteNet,time)

load LatWeights.mat

difNHSH = [];
NH = [];
SH = [];

for j=1:time
    latFluxes = sum(permuteNet(:,:,j),2).*LatWeights(:,2);
    difNHSH = [difNHSH sum(latFluxes(91:end)-latFluxes(1:90))];
    NH = [NH sum(latFluxes(91:end))];
    SH = [SH sum(latFluxes(1:90))];
end

% 
% %%%%%%%%%%%%IGNORE THIS%%%%%%%%%%%%%%%%%
set(gca,'FontSize',20)
windowSize = 12;
AllRunningMeans = filter([31/365 28/365 31/365 30/365 31/365 30/365 31/365 31/365 30/365 31/365 30/365 31/365],1,difNHSH) %but wouldn't this make jan count as april for some starts?
plot(AllRunningMeans(12:end))
grid on;
set(gca,'xtick',[10:12:time])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year where mean ends on');
ylabel('NH-SH Difference in Total Heat Flux (rsdt - rsutcs - rlutcs)');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AllRunningMeans = runmean(difNHSH,12);
% plot(AllRunningMeans)
% grid on;
% set(gca,'xtick',[10:12:time])
% set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010})
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

valid_moving_weighted_avg = moving_weighted_avg(12:end); 
%index 1 is march 2000 to february 2001
%index 13 is march 2001 to february 2002
%index 11 is january 2001 to december 2001

size(valid_moving_weighted_avg)
%subplot(2,1,1)

set(gca,'FontSize',20)
plot(valid_moving_weighted_avg)
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year end');
ylabel('NH-SH Difference in Total Heat Flux (rsdt - rsutcs - rlutcs)');
title('Running Mean in NH-SH Flux Average (by year end)')

set(gca,'FontSize',20)
plot(NHmoving_weighted_avg(12:end)-mean(NHmoving_weighted_avg(12:end)),'LineWidth',3)
hold on;
plot(SHmoving_weighted_avg(12:end)-mean(SHmoving_weighted_avg(12:end)), 'g','LineWidth',3)
plot(valid_moving_weighted_avg-mean(valid_moving_weighted_avg),'r','LineWidth',3)
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year end');
ylabel('NH-SH Difference in Total Heat Flux');
set(legend('NH','SH','NH-SH'),'Location','BestOutside')

NH_MWA = NHmoving_weighted_avg(12:end);
SH_MWA = SHmoving_weighted_avg(12:end);

MWA = [valid_moving_weighted_avg;NH_MWA;SH_MWA];

%%%OLD METHOD TO COMPARE WITH- DON'T USE%%%
% subplot(2,1,2)
% set(gca,'FontSize',20)
% meanNHSH=[];
% for jj=1:119
%     meanNHSH = [meanNHSH mean(difNHSH(jj:jj+11))];
% end
% plot(meanNHSH)
% grid on;
% set(gca,'xtick',[6:12:125])
% set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010})
% xlabel('Year where mean is centered on');
% ylabel('NH-SH Difference in Total Heat Flux (rsdt - rsutcs - rlutcs)');

end