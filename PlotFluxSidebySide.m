function [MaxCorr, RealMonthOfMaxCorr, TimeLaggedCorr]= PlotFluxSidebySide(Flux1,Flux2,Name1,Name2,MonthFilterSize,PlotYesNo)
time = length(Flux1);

MaxLag = 20;
TimeLeadCorr = zeros(1,MaxLag);
TimeLagCorr = zeros(1,MaxLag);

for i=1:MaxLag
    TimeLeadCorr(i) = corr(Flux1(1:end-i)',Flux2(i+1:end)');
    TimeLeadRSq(i) = RSquared(Flux1(1:end-i)',Flux2(i+1:end)');
    TimeLagCorr(i)=corr(Flux1(i+1:end)',Flux2(1:end-i)');
    TimeLagRSq(i)=RSquared(Flux1(i+1:end)',Flux2(1:end-i)');
end
TimeLeadCorr = fliplr(TimeLeadCorr);
tempCorr = corr(Flux1',Flux2');
tempRSq = RSquared(Flux1',Flux2');

TimeLaggedCorr = [TimeLeadCorr tempCorr TimeLagCorr];
TimeLaggedRSq = [TimeLeadRSq tempRSq TimeLagRSq];

% TimeLaggedCorr(MaxLag+1)
MonthsBeforeAfter = [-MaxLag:MaxLag];
% plot(MonthsBeforeAfter,TimeLaggedCorr)
MonthOfMaxCorr =  find(max(abs(TimeLaggedCorr))==abs(TimeLaggedCorr));
MaxCorr = TimeLaggedCorr(MonthOfMaxCorr);
RealMonthOfMaxCorr = MonthsBeforeAfter(MonthOfMaxCorr);

if (RealMonthOfMaxCorr < 0)
LagString = [' at ', Name1, ' leading ', Name2, ' by ', num2str(-RealMonthOfMaxCorr) ' months'];
else
LagString = [' at ', Name2, ' leading ', Name1, ' by ', num2str(RealMonthOfMaxCorr) ' months'];
end


if PlotYesNo==1
	
	[AX,H1,H2] = plotyy(1:length(Flux1),[Flux1],1:length(Flux2),[Flux2]);
	% [AX,H1,H2] = plotyy(1:length(NH),[NH' SH' HemSum'],1:length(NH),HemDif);
	grid on;
	set(gca,'xtick',12-2-(MonthFilterSize-1):12:time)
	set(AX(2),'XTickLabel',[])
	set(gca,'XTickLabel',2000:2012)
	xlabel('Year End');
	set(get(AX(1),'Ylabel'),'FontSize',20,'String',Name1) 
	set(get(AX(2),'Ylabel'),'FontSize',20,'String',Name2,'FontSize',20) 
	
	corrcof = getfield(corrcoef(Flux1,Flux2),{1,2});
	set(legend([Name1],[Name2]));
	% set(legend(['NH \mu = ', num2str(NHMean)],['SH \mu = ', num2str(SHMean)] ,['Global \mu = ' num2str(GlobalMean)],sprintf(['NH-SH \mu = ', num2str(HemDifMean), '\n Corr = ', num2str(corrcof)])),'Location','BestOutside')
	
	set(H1,'linewidth',4)
	set(H2,'LineStyle','--')
	set(H2,'linewidth',3)% to change the first line
	set(AX,'FontSize',20)
	title(['Max Time-Lagged Corr is ', num2str(MaxCorr), LagString])
	set(gcf, 'Units','inches', 'Position',[0 0 20 10])
	set(gca, 'Units','inches', 'Position',[1 1 16 8])
	set(gca,'GridLineStyle','--')
	set(gcf,'paperposition',[0 0 20 10])
	print(gcf,'-dpng','-r300',['Compare',Name1,'vs',Name2,'.png']);
	saveas(gcf,['Compare',Name1,'vs',Name2,'.fig'],'fig')
	hold off;
	
	set(gca,'linewidth',4)
	set(gca,'LineStyle','--')
	plot(MonthsBeforeAfter, TimeLaggedCorr)
	title(['Time-Lagged Corr of ', Name2, ' leading ',Name1])
	xlabel('Months of Lag')
	set(gca,'GridLineStyle','--')
	set(gcf,'paperposition',[0 0 20 10])
	print(gcf,'-dpng','-r300',['TimeLagCorr',Name1,'vs',Name2,'.png']);
	saveas(gcf,['TimeLagCorr',Name1,'vs',Name2,'.fig'],'fig')
	end
end
