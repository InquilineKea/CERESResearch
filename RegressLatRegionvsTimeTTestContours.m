function [Beta,pValues] = RegressLatRegionvsTimeTTestContours(Flux1,TimeSeries,Name1,Name2,LowerLat,HigherLat)
load LatWeights.mat
LatWeightsAlone = LatWeights(:,2);

%time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));
lats = size(Flux1,1);
longs = size(Flux1,2);
time=size(Flux1,3);
Factor = 180/lats;
lat = ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','lat');

%TimeSeries = NetDifMonthTotal;
%Flux1 = NetClimatologySubtracted;
TimeSeries3D = ones(1,time);
TimeSeries3D(1,:) = TimeSeries;

Flux1 = squeeze(sum(Flux1,2));
Flux1 = Flux1.*repmat(LatWeightsAlone,[1 156]);
TimeSeries3D = repmat(TimeSeries3D,[size(Flux1(LatWeights(:,1) < HigherLat & LatWeights(:,1) > LowerLat,:),1) 1]);
Flux1 = sum(Flux1(LatWeights(:,1) < HigherLat & LatWeights(:,1) > LowerLat,:));
% Flux1 = Flux1';
TimeSeries3D = TimeSeries;
%Beta = bsxfun(@rdivide,time*(sum(Flux1.*TimeSeries3D,3))-sum(Flux1,3).*sum(TimeSeries3D,3),time*sum(TimeSeries3D.*TimeSeries3D,3)-sum(TimeSeries3D,3).*sum(TimeSeries3D,3))*std(TimeSeries);
%multiply by std(TimeSeries) in end so that you come out with units of
Sxy = sum(Flux1.*TimeSeries3D,2);
Sxx = sum(TimeSeries3D.*TimeSeries3D,2);
Syy = sum(Flux1.*Flux1,2);
SSE = Syy-Sxy.*Sxy./Sxx;
Sy = sum(Flux1,2);
Sx = sum(TimeSeries);
Beta = bsxfun(@rdivide,time*(Sxy)-Sx*Sy,time*Sxx-Sx.*Sx);
[x,y] = autocorr(TimeSeries);
AutoCorrTime = length(x(x > 1/exp(1))) + 0.5;
DOF = time/(2*AutoCorrTime);
tScore = Beta*sqrt(DOF-2)./sqrt(SSE/sum((TimeSeries-mean(TimeSeries)).^2));
tTestAlpha = 0.05;
pValues = tcdf(tScore,DOF-1);
end