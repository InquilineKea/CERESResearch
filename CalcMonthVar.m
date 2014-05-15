function [MonthVar] = CalcMonthVar(difNHSH,month)

startMonth = month;
if month < 3
    startMonth = month + 12;
end
MonthVar = var(difNHSH(startMonth:12:end));


end