function [MonthVar] = MonthTrend(difNHSH,month)

%still need to offset january, february to the right by one year

%Years = [2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012];
Years = [2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011];
cmap = jet;

startMonth = month;
set(gca,'FontSize',20)
if month < 3
    startMonth = month + 12;
    plot(1:length(difNHSH(startMonth:12:end))+1,[NaN difNHSH(startMonth:12:end) - mean(difNHSH(startMonth:12:end))],'color',cmap((mod(month+1:month+1,12)+1)*5,:),'LineWidth',3)
else
    plot(1:length(difNHSH(startMonth:12:end)),difNHSH(startMonth:12:end) - mean(difNHSH(startMonth:12:end)),'color',cmap((mod(month+1:month+1,12)+1)*5,:),'LineWidth',3)
end
grid on;
%set(gca,'XTickLabel',{2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011})
set(gca,'xtick',1:length(difNHSH(startMonth:12:end)))
set(gca,'XTickLabel',{2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
xlabel('Year');
ylabel('Flux Difference');

MonthVar = var(difNHSH(startMonth:12:end));

end