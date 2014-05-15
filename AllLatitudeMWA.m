function [MWA] = AllLatitudeMWA(Flux)

inputname(1)

lat = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');
load LatWeights.mat

time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));

MWA = zeros(180,360,time);
%Flux=Net;

days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
allMonthWeights = repmat(days_per_month,1,12);
allMonthWeights = [31 30 31 30 31 31 30 31 30 31 allMonthWeights 31 28];
allMonthWeights(48) = 29; %leap year, http://www.wolframalpha.com/input/?i=months+between+march+2000+and+february+2004
allMonthWeights(96) = 29;

%1 is march 2000. leap years in february 2004, 2008

meanLatFlux =squeeze(mean(Flux,2));
testt = bsxfun(@times,meanLatFlux,LatWeights(:,2)); %180x156
% plot(testt(:,1))
% plot(testt(90,:))


% hold on;
% plot(meanLatFlux(5,:))
% plot(meanLatFlux(90,:))

size(repmat(allMonthWeights,[180 1])'); %156x180

moving_sum = @(n, x) filter(ones(1,n), 1, x);
MWA= bsxfun(@rdivide,moving_sum(12, testt' .* repmat(allMonthWeights,[180 1])') ...
    ,moving_sum(12,repmat(allMonthWeights,[180 1])')); %or Mean Lat Flux here

%this is MWA for each latitude, but not MWA for each tile.. which is more
%difficult...

%moving_sum(12,repmat(allMonthWeights,[180 1])') is 156x180

%start: 9:23 PM

months=12;
%EachLatTimeSeries = moving_sum(months, testt(2,:));

% plot(EachLatTimeSeries(12:end)) %moving sum does poorly with multidimensions...
%like 180 x 156. it can only deal with 1 x 156...

% EachLatTimeSeries = moving_sum(months, testt'); %156x180
%moving sum does poorly with multidimensions...
cmap = jet;


%%%%%%VERSION 1%%%%%%%%%%%%%
% set(gca,'FontSize',20)
% j = 1;
% for i=1:3:180
%     plot(EachLatTimeSeries(12:end,i),'color',cmap(j,:))
%     hold on;
%     j = j + 1;
% end
% grid on;
% set(gca,'xtick',[11:12:119+26])
% set(gca,'XTickLabel',{2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
% xlabel('Year end');

%%%%%%%%%%VERSION 2 %%%%%%%%

% set(gca,'FontSize',20)
% j = 1;
% for i=1:10:171
%     plot(mean(EachLatTimeSeries(12:end,i:i+9),2)-mean(mean(EachLatTimeSeries(12:end,i:i+9),2)),'color',cmap(j*3,:),'LineWidth',3)
%     hold on;
%     j = j + 1;
% end
% grid on;
% set(gca,'xtick',[11:12:119+26])
% set(gca,'XTickLabel',{2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
% xlabel('Year end');
% ylabel('Net Flux (in Watts)');
% set(legend('90S','80S','70S','60S','50S','40S','30S','20S','10S','0S','10N','20N','30N','40N','50N','60N','70N','80N'),'Location','BestOutside') %NorthWestOutside
% title('Moving Weighted Average Net Flux per Latitude Band, after subtraction by its relative mean (Watts)')

%%%%%Version MWA%%%%%%%%%%%%%%%%%%%

set(gca,'FontSize',20)
j = 1;
for i=1:10:171
    %zonal mean subtraction below
    plot(sum(MWA(12:end,i:i+9),2)-mean(sum(MWA(12:end,i:i+9),2)),'color',cmap(j*3,:),'LineWidth',3)%zonal mean subtraction
    %plot(mean(MWA(12:end,i:i+9),2),'color',cmap(j*3,:),'LineWidth',3)
    hold on;
    j = j + 1;
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year end');
ylabel('Net Flux (in Watts)');
lh = legend('90-80S','80-70S','70-60S','60-50S','50-40S','40-30S','30-20S','20-10S','10-0S','0-10N','10-20N','20-30N','30-40N','40-50N','50-60N','60-70N','70-80N','80-90N');
set(lh,'Location','BestOutside') %NorthWestOutside
title('Moving Weighted Average Net Flux per Latitude Band, after mean-subtraction (Watts)')
   v = get(lh,'title');
   set(v,'string','Latitude','FontSize',14);


%%%%%%%%%%%%%%%%%%

% plot(sum(EachLatTimeSeries(1:89,:)))
% 
% bah2 = moving_sum(months, SH')/months;
% plot(bah2(12:end))
% 
% 
% MWA =MWA(:,12:end);
% 
% subplot(1,2,1)
% hold on;
% for i=1:30:180
%     plot(meanLatFlux(i,:))
% end
% 
% subplot(1,2,2)
% hold on;
% for i=1:30:180
%     plot(MWA(i,:))
% end

end