function [i] = PlotDecomposedVar(Var,month)

i = 1;
hold on;
set(gca,'FontSize',20)
plot(Var(1,:),'LineWidth',3)
plot(Var(2,:),'g','LineWidth',3)
plot(Var(3,:),'r','LineWidth',3)

set(legend('NH-SH Var','NH Var','SH Var'),'Location','BestOutside')
title('Variances of each NH-SH component by month')
set(gca,'xtick',1:12)


end