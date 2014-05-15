% 3 datasets:
% 
% Argo
% Interpolated Argo
% Seal

%switch directory first!

seals = catpad(3,comet,blitzer,Soern,undin,viken,dasher,donner,bernt,aspasia);

[numDepths,numFields,numSeals] = size(seals);

color=[':g',':b',':r',':y',':c',':m',':r',':k',':k'];
colorMarker = {'gx','bx','ro','yd','c^','mv','rs','k<','ch'};

for i=1:numSeals
    plot3(seals(:,3,i),seals(:,4,i),seals(:,5,i),colorMarker{i})
    hold on;
    grid on;
end

set(legend('comet','blitzer','Soern','undin','viken','dasher','donner','bernt','aspasia'),'Location','BestOutside') %NorthWestOutside

%logical array indexing is amazing
 seals(seals(:,3,1) < 10,:,1)
  
 Pdifs = [];
 Tdifs = [];
 Sdifs = [];
 
  SealDisplacements = [];
  SealIndexFirst = [];
  SealIndexLast = [];
  SealSurfaceTDifference = [];
  SealSurfaceSDifference = [];
 
 for j=1:numSeals
     [a,SealIndexFirst] = unique(seals(:,1,j)*999999999+seals(:,2,j),'first');
     [a,SealIndexLast] = unique(seals(:,1,j)*999999999+seals(:,2,j),'last');
     SealDisplacements = catpad(2,sqrt(power(abs(diff(seals(SealIndexFirst,1,j))),2) + power(abs(diff(seals(SealIndexFirst,2,j))),2)),SealDisplacements);
     SealSurfaceTDifference = catpad(2,diff(seals(SealIndexFirst,4,j)),SealSurfaceTDifference);
     SealSurfaceSDifference = catpad(2,diff(seals(SealIndexFirst,5,j)),SealSurfaceSDifference);
     for i=1:size(SealIndexFirst)*[1,0]'
        Pdifs = catpad(2, Pdifs, diff(seals(SealIndexFirst(i):SealIndexLast(i),3,j)));
        Tdifs = catpad(2, Tdifs, diff(seals(SealIndexFirst(i):SealIndexLast(i),4,j)));
        Sdifs = catpad(2, Sdifs, diff(seals(SealIndexFirst(i):SealIndexLast(i),5,j)));
     end

 end
 
 SealTSurfaceGrad = SealSurfaceTDifference./SealDisplacements;
 SealSSurfaceGrad = SealSurfaceSDifference./SealDisplacements;
 
 hist(SealTSurfaceGrad(abs(SealTSurfaceGrad) < 5),100)
 hist(SealSSurfaceGrad(abs(SealSSurfaceGrad) < 5),100)
 
 AllTGrads = Tdifs./Pdifs;
 AllSGrads = Sdifs./Pdifs;
 AllTGrads = AllTGrads(isfinite(AllTGrads));
 AllSGrads = AllSGrads(isfinite(AllSGrads));
 hist(AllTGrads,100)
 [mean(AllTGrads) var(AllTGrads)]
 
 SealDisplacements = [];
 
 sqrt(power(diff(unique(seal(:,1,1))),2)+power(diff(unique(seal(:,2,1))),2))
 
 for j=1:numSeals 
    SealDisplacements = catpad(2,sqrt(power(diff(unique(seal(:,1,j))),2)+power(diff(unique(seal(:,1,j))),2)),SealDisplacements);
 end
 
%       [a,SealIndexFirst] = unique(seals(:,1,1)*999999999+seals(:,2,1),'first');
%      [a,SealIndexLast] = unique(seals(:,1,1)*999999999+seals(:,2,1),'last');
% 
%      for i=1:size(SealIndexFirst)*[1,0]'
%         Pdifs = catpad(2, diff(seals(SealIndexFirst(i):SealIndexLast(i),3,1)),Pdifs);
%         Tdifs = catpad(2, diff(seals(SealIndexFirst(i):SealIndexLast(i),4,1)),Tdifs);
%      end

 
dirls=dir('*prof.nc')

ArgoTemps = [];
ArgoPressures = [];
ArgoSalinities = [];
ArgoLongLats = [];

for filenum=1:length(dirls)
  ArgoTemps = catpad(3, ncread(dirls(filenum).name,'TEMP'),ArgoTemps);
  ArgoPressures = catpad(3, ncread(dirls(filenum).name,'PRES'),ArgoPressures);
  ArgoSalinities= catpad(3, ncread(dirls(filenum).name,'PSAL'),ArgoSalinities);
  ArgoLongLats= catpad(3, [ncread(dirls(filenum).name,'LONGITUDE') ncread(dirls(filenum).name,'LATITUDE')],ArgoLongLats);
end

%"Conversion to double from cell is not possible." if no {}

for i=1:length(dirls)
    plot(ArgoLongLats(:,1,i),ArgoLongLats(:,2,i))
    %scatter(ArgoLongLats(:,1,i),ArgoLongLats(:,2,i),[],ArgoTemps(1,:,i)');colorbar %surface Ts
    %scatter(ArgoLongLats(:,1,i),ArgoLongLats(:,2,i),[],ArgoSalinities(1,:,i)');colorbar %surface Ts
    hold all;
    plot(seals(:,2,i),seals(:,1,i),colorMarker{i})
    grid on;
end

xlabel('longitude')
ylabel('latitude')
set(gca,'FontSize',20)
set(get(gca,'xlabel'),'string','Longitude','fontsize',20)
set(get(gca,'ylabel'),'string','Latitude','fontsize',20)
title('Surface Distribution of Seals (Points) Compared with Argo Floats (Lines)')


%plot(ArgoPressures(:,1,1),ArgoTemps(:,1,1),'x') %a single vert profile
%plot(ArgoPressures(:,1:10,1),ArgoTemps(:,1:10,1),'x') %10 vert profiles

for i=1:length(dirls)
    plot(ArgoPressures(:,1:10,i),ArgoTemps(:,1:10,i),'x') %10 vert profiles
    subplot(3,3,i)
end

%vertical temperature gradients: diff(ArgoTemps(:,:,1))./diff(ArgoPressures(:,:,1))

hist(ArgoTemps(1,:,1))

hist(diff(ArgoTemps(1,:,1))')

%horizontal surface temperature gradient in absolute displacement direction

AllArgoSurfaceTGrads = [];
AllArgoSurfaceSGrads = [];

for i=1:length(dirls)
    AllArgoSurfaceTGrads = catpad(2,diff(ArgoTemps(1,:,i))'./sqrt(abs(diff(ArgoLongLats(:,1,i))).^2+abs(diff(ArgoLongLats(:,2,i))).^2),AllArgoSurfaceTGrads);
    AllArgoSurfaceSGrads = catpad(2,diff(ArgoSalinities(1,:,i))'./sqrt(abs(diff(ArgoLongLats(:,1,i))).^2+abs(diff(ArgoLongLats(:,2,i))).^2),AllArgoSurfaceSGrads);
    %subplot(3,3,i)
    %hist(diff(ArgoTemps(1,:,i))'./sqrt(abs(diff(ArgoLongLats(:,1,i))).^2+abs(diff(ArgoLongLats(:,2,i))).^2))
    %hist(diff(ArgoSalinities(1,:,i))'./sqrt(abs(diff(ArgoLongLats(:,1,i))).^2+abs(diff(ArgoLongLats(:,2,i))).^2))
   
end

[n1, xout1] =  hist(SealTSurfaceGrad(abs(SealTSurfaceGrad) < 5),100);
bar(xout1,n1,'r'); grid; hold
[n2, xout2] = hist(AllArgoSurfaceTGrads(abs(AllArgoSurfaceTGrads) < 5),100);
bar(xout2,n2,'b');

set(legend('Seals','Floats'),'Location','BestOutside')
set(gca,'FontSize',20)
set(get(gca,'xlabel'),'string','Surface Displacement Temperature Gradients (dT/(d\phi^2+d\lambda^2)^{0.5})','fontsize',20)
set(get(gca,'ylabel'),'string','Frequency','fontsize',20)
title('Surface Temperature Gradients of Seals vs Floats')

[n1, xout1] =  hist(SealSSurfaceGrad(abs(SealSSurfaceGrad) < 1),100);
bar(xout1,n1,'r'); grid; hold
[n2, xout2] = hist(AllArgoSurfaceSGrads(abs(AllArgoSurfaceSGrads) < 1),100);
bar(xout2,n2,'b');

set(legend('Seals','Floats'),'Location','BestOutside')
set(gca,'FontSize',20)
set(get(gca,'xlabel'),'string','Surface Displacement Salinity Gradients (dS/(d\phi^2+d\lambda^2)^{0.5})','fontsize',20)
set(get(gca,'ylabel'),'string','Frequency','fontsize',20)
title('Surface Salinity Gradients of Seals vs Floats')


errorbar(nanmean(SealTSurfaceGrad),nansem(SealTSurfaceGrad))
errorbar(nanmean(AllArgoSurfaceTGrads),nansem(AllArgoSurfaceTGrads))


%you can actually make a movie going from the top layer of ArgoTemps to the
%bottom layer! or try multiple subplots

set(legend('comet','blitzer','Soern','undin','viken','dasher','donner','bernt','aspasia'),'Location','BestOutside') %NorthWestOutside

% a = ArgoTemps(ArgoPressures < 50)./ArgoPressures(ArgoPressures < 50);
% a(isinf(a))=[];
% hist(a)
% b = ArgoTemps(50 < ArgoPressures < 100)./ArgoPressures(50 < ArgoPressures < 100);
% b(isinf(b))=[];
% hist(b)
% c = ArgoTemps(100 < ArgoPressures < 200)./ArgoPressures(100 < ArgoPressures < 200);
% c(isinf(c))=[];
% hist(c)
% d = ArgoTemps(200 < ArgoPressures < 500)./ArgoPressures(200 < ArgoPressures < 500);
% d(isinf(d))=[];
% hist(d)

AllTDiff = diff(ArgoTemps)./diff(ArgoPressures);
AllSDiff = diff(ArgoSalinities)./diff(ArgoPressures);
AllTDiff = AllTDiff(isfinite(AllTDiff));
AllSDiff = AllSDiff(isfinite(AllSDiff));

hist(AllTDiff(abs(AllTDiff) < .05),100)

[n2, xout2] = histnorm(AllTDiff(abs(AllTDiff) < .05),100);
bar(xout2,n2,'b'); grid; hold
[n1, xout1] =  histnorm(AllTGrads(abs(AllTGrads) < .05),100);
bar(xout1,n1,'r'); 
parentHandle = bar(xout2,n2,'b');
childHandle = get(parentHandle,'Children');
set(childHandle,'FaceAlpha',0.7); % 0 = transparent, 1 = opaque.

set(legend('Floats','Seals'),'Location','BestOutside')
set(gca,'FontSize',20)
set(get(gca,'xlabel'),'string','Vertical Displacement Temperature Gradients (dT/dP)','fontsize',20)
set(get(gca,'ylabel'),'string','Normalized Frequency','fontsize',20)
title('Vertical Displacement Temperature Gradients of Seals vs Floats')


[n2, xout2] = histnorm(AllSDiff(abs(AllSDiff) < 0.1),100);
bar(xout2,n2,'b'); grid; hold
[n1, xout1] =  histnorm(AllSGrads(abs(AllSGrads) < 0.1),100);
bar(xout1,n1,'r');

set(legend('Floats','Seals'),'Location','BestOutside')
set(gca,'FontSize',20)
set(get(gca,'xlabel'),'string','Vertical Displacement Salinity Gradients (dS/dP)','fontsize',20)
set(get(gca,'ylabel'),'string','Normalized Frequency','fontsize',20)
title('Vertical Displacement Salinity Gradients of Seals vs Floats')

StepSize = 50;
sealIndexCount = [];
argoIndexCount = [];

for x=StepSize:StepSize:700
    grid on;
    sealIndex = Tdifs(cumsum(Pdifs)+5 < x & cumsum(Pdifs)+5 > x-StepSize)./Pdifs(cumsum(Pdifs)+5 < x & cumsum(Pdifs)+5 > x-StepSize);
    sealIndex(abs(sealIndex)>10)=[];
    sealIndexCount = [sealIndexCount length(sealIndex)];
    errorbar(x,nanmedian(sealIndex),nanstd(sealIndex),'<r','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
    hold all;
    %length(~isnan(diff(ArgoTemps(x > ArgoPressures & ArgoPressures > x-StepSize))./diff(ArgoPressures(x > ArgoPressures & ArgoPressures > x-StepSize))))
    argoIndex = diff(ArgoTemps(x > ArgoPressures & ArgoPressures > x-StepSize))./diff(ArgoPressures(x > ArgoPressures & ArgoPressures > x-StepSize));
    argoIndex = argoIndex(isfinite(argoIndex));
    argoIndex(abs(argoIndex)>10)=[];
    argoIndexCount = [argoIndexCount length(argoIndex)];
    [nanmedian(argoIndex), nansem(argoIndex)];
    errorbar(x,nanmedian(argoIndex), nansem(argoIndex),'ob','MarkerFaceColor','b','MarkerSize',10)
end

set(legend('Seals','Floats'),'Location','BestOutside')
set(gca,'FontSize',20)
set(get(gca,'xlabel'),'string','Pressure (decibars)','fontsize',20)
set(get(gca,'ylabel'),'string','\Delta Temperature (C)/\Delta Pressure (decibars)','fontsize',20)
title('Pressure Binning the Vertical Temperature Gradients of Seals vs Floats')



plot(StepSize:StepSize:2000,sealIndexCount)
hold all;
grid on;
plot(StepSize:StepSize:2000,argoIndexCount)

for x=StepSize:StepSize:1200
    grid on;
    sealIndex = Sdifs(cumsum(Pdifs)+5 < x & cumsum(Pdifs)+5 > x-StepSize)./Pdifs(cumsum(Pdifs)+5 < x & cumsum(Pdifs)+5 > x-StepSize);
    sealIndex(abs(sealIndex)>10)=[];
    sealIndexCount = [sealIndexCount length(sealIndex)];
    errorbar(x,nanmedian(sealIndex),nanstd(sealIndex),'<r','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
    hold all;
    %length(~isnan(diff(ArgoTemps(x > ArgoPressures & ArgoPressures > x-StepSize))./diff(ArgoPressures(x > ArgoPressures & ArgoPressures > x-StepSize))))
    argoIndex = diff(ArgoSalinities(x > ArgoPressures & ArgoPressures > x-StepSize))./diff(ArgoPressures(x > ArgoPressures & ArgoPressures > x-StepSize));
    argoIndex = argoIndex(isfinite(argoIndex));
    argoIndex(abs(argoIndex)>10)=[];
    argoIndexCount = [argoIndexCount length(argoIndex)];
    [nanmedian(argoIndex), nansem(argoIndex)];
    errorbar(x,nanmedian(argoIndex), nansem(argoIndex),'ob','MarkerFaceColor','b','MarkerSize',10)
end

set(legend('Seals','Floats'),'Location','BestOutside')
set(gca,'FontSize',20)
set(get(gca,'xlabel'),'string','Pressure (decibars)','fontsize',20)
set(get(gca,'ylabel'),'string','\Delta Salinity (psu)/\Delta Pressure (decibars)','fontsize',20)
title('Pressure Binning the Vertical Salinity Gradients of Seals vs Floats')



% hist(diff(ArgoTemps)./diff(ArgoPressures))

T50 = ArgoTemps;
P50 = ArgoPressures;
T50(P50 > 1000) = NaN;
P50(P50 > 1000) = NaN;
d50 = diff(T50)./diff(P50);

hist(d50(abs(d50)<5),100)

hist(AllTGrads,100)
hold on;
hist(d50(abs(d50)<.1),100)

T100 = ArgoTemps;
P100 = ArgoPressures;
T100(P100 < 50 | P100 > 100) = NaN;
P100(P100 < 50 | P100 > 100) = NaN;
d100 = diff(T100)./diff(P100);
hist(d100(isfinite(d100)),100)
hist(d100(abs(d100)<20),100)



 ncdisp(dirls(2).name)