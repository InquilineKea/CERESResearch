plot([corr(Net49_90N_Month',Net49_90S_Month') corr(Net49_90N_Month',Net30_49S_Month') ...
    corr(Net49_90N_Month',Net15_30S_Month') corr(Net49_90N_Month',Net0_15S_Month') ...
    corr(Net49_90N_Month',Net0_15N_Month') corr(Net49_90N_Month',Net15_30N_Month')...
    corr(Net49_90N_Month',Net30_49N_Month') corr(Net49_90N_Month',Net49_90N_Month')],'rs','MarkerFaceColor','g','markersize', 20)
set(gca,'FontSize',20)
title(['Correlation Coefficients between each Latitudinal Band and ', strrep('49_90N','_','-') ,' Difference net flux'])
xlabel('Latitudinal Band')
ylabel('Correlation Coefficient')
grid on;
set(gca,'XTickLabel',{'90S-49S','49S-30S','30S-15S','15S-0S','0N-15N','15N-30N','30N-49N','49N-90N'})
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['NetLatitudinal49_90NCorrelation','.png']);
hold off;
