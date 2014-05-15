function [Beta,pValues, Alpha] = RegressLatvsTimeTTestContours(Flux1,TimeSeries,Name1,Name2)
load LatWeights

%time = length(ncread('rlutcs_CERES-EBAF_L3B_Ed2-7_200003-201302.nc','time'));
lats = size(Flux1,1);
longs = size(Flux1,2);
time=size(Flux1,3);
Factor = 180/lats;

%TimeSeries = NetDifMonthTotal;
%Flux1 = NetClimatologySubtracted;
TimeSeries3D = ones(1,time);
TimeSeries3D(1,:) = TimeSeries;
TimeSeries3D = repmat(TimeSeries3D,[180/Factor 1]);

Flux1 = squeeze(mean(Flux1,2));

%Beta = bsxfun(@rdivide,time*(sum(Flux1.*TimeSeries3D,3))-sum(Flux1,3).*sum(TimeSeries3D,3),time*sum(TimeSeries3D.*TimeSeries3D,3)-sum(TimeSeries3D,3).*sum(TimeSeries3D,3))*std(TimeSeries);
%multiply by std(TimeSeries) in end so that you come out with units of
Sxy = sum(Flux1.*TimeSeries3D,2);
Sxx = sum(TimeSeries3D.*TimeSeries3D,2);
Syy = sum(Flux1.*Flux1,2);
SSE = Syy-Sxy.*Sxy./Sxx;
Sy = sum(Flux1,2);
Sx = sum(TimeSeries);
Beta = bsxfun(@rdivide,time*(Sxy)-Sx*Sy,time*Sxx-Sx.*Sx);
TotalBeta = sum(Beta);
Alpha = 1/time*(Sy-Beta*Sx);
% Lats2Divide = [0,30,50,70,90];
% BetaSection = zeros(length(Lats2Divide)-1,1);
% for i=1:length(Lats2Divide)-1
%     BetaSection(i) = sum(Beta(LatWeights(:,1) > Lats2Divide(i) & LatWeights(:,1) < Lats2Divide(i+1)));
% end
% 
[x,y] = autocorr(TimeSeries);
AutoCorrTime = length(x(x > 1/exp(1))) + 0.5;
DOF = time/(2*AutoCorrTime);
tScore = Beta*sqrt(DOF-2)./sqrt(SSE/sum((TimeSeries-mean(TimeSeries)).^2));
tTestAlpha = 0.05;
pValues = tcdf(tScore,DOF-1);

% tScoreSection = BetaSection*sqrt(DOF-2)./sqrt(SSE/sum((TimeSeries-mean(TimeSeries)).^2));
% pValuesSection = tcdf(tScoreSection,DOF-1);


end