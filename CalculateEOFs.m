function [U,S,V] = CalculateEOFs(Flux,VarName,EOF)

load LatWeights.mat
Factor = 180/size(Flux,1);
x=1:1:180;
xq = 1:180/size(Flux,1):180;

FluxWeighted=bsxfun(@times,Flux,sqrt(cosd(interp1(x,LatWeights(:,1),xq)))');
%FluxWeighted=MonthlyFluxDepartures.LW;
ToSVD = reshape(FluxWeighted,size(Flux,1)*size(Flux,2),size(Flux,3))';
%clear MonthlyFluxDepartures;
[U,S,V] = svds(ToSVD);
Regress1DTimeSeriesMapWithTTestContours(Flux, U(:,1),VarName,['EOF1 of ', VarName])
hold off;
Regress1DTimeSeriesMapWithTTestContours(Flux, U(:,2),VarName,['EOF2 of ', VarName])
Regress1DTimeSeriesMapWithTTestContours(Flux, U(:,3),VarName,['EOF3 of ', VarName])


end