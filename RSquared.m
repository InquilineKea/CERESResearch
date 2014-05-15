function RSquare = RSquared(x,y)

SSTotal = sum((y-mean(y)).^2);
Beta = (sum(x.*y)-sum(x)*sum(y)/length(y))/(sum(x.^2) - sum(x)^2/length(y));
Alpha = mean(y)-Beta*mean(x);
RegressedY = Alpha+Beta*x;
SSRes = sum((y-RegressedY).^2);
RSquare= 1 - SSRes./SSTotal;
end