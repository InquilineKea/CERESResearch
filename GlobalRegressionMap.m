function [Beta,NHMinusSHIndividualLatitudes,FluxSumsAllLatitudeLines] = GlobalRegressionMap(Flux1,Flux2, Var1Name,Var2Name)
lats = min(size(Flux1,1),size(Flux2,1));
longs = min(size(Flux1,2),size(Flux2,2));
time=size(Flux1,3);
Factor = 180/lats;

%GlobalBetaCorr = bsxfun(@rdivide,sum(Flux1.*Flux2,3),(sqrt(sum(Flux1.^2,3)).*sqrt(sum(Flux2.^2,3))));
if size(Flux1,1) < size(Flux2,1)
    Matrix2Resize = Flux2;
    OtherMatrix = Flux1;
elseif size(Flux2,1) < size(Flux1,1)
    Matrix2Resize = Flux1;
    OtherMatrix = Flux2;
end

if size(Flux2,1) ~= size(Flux1,1)
    E=zeros(size(Matrix2Resize,1)/Factor,size(Matrix2Resize,2)/Factor,size(Matrix2Resize,3));
    for depth=1:size(Matrix2Resize,3)
      E(:,:,depth)=imresize(Matrix2Resize(:,:,depth),1/Factor);
    end
    Matrix2Resize = E;
    if size(Flux1,1) < size(Flux2,1)
        Flux2 = Matrix2Resize;
    else
        Flux1 = Matrix2Resize;
    end
end

Beta = bsxfun(@rdivide,time*(sum(Flux1.*Flux2,3))-sum(Flux1,3).*sum(Flux2,3),time*sum(Flux2.*Flux2,3)-sum(Flux2,3).*sum(Flux2,3));
load LatWeights.mat
RegressMean = mean(mean(resizem(Beta,Factor),2).*LatWeights(:,2));

FluxSumsAllLatitudeLines = sum(resizem(Beta,Factor),2).*LatWeights(:,2);
NHPlusSHTotalFluxContribution = sum(FluxSumsAllLatitudeLines);
NHTotalFluxContribution = sum(FluxSumsAllLatitudeLines(end/2+1:end));
SHTotalFluxContribution = sum(FluxSumsAllLatitudeLines(1:end/2));
NHMinusSHTotalFluxContribution = NHTotalFluxContribution - SHTotalFluxContribution;
NHMinusSHIndividualLatitudes = FluxSumsAllLatitudeLines(end/2+1:end)-flipud(FluxSumsAllLatitudeLines(1:end/2));

HemisphericFluxes = [NHTotalFluxContribution, SHTotalFluxContribution, NHMinusSHTotalFluxContribution,NHPlusSHTotalFluxContribution];

load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(resizem(Beta,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
% caxis([-max(max(Beta)) max(max(Beta))])
MakeLizMap
colormap(lizmap)
title(['Regression Coefficient of ',Var1Name,' Regressed Over ',Var2Name, '. Mean = ', num2str(RegressMean)])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['GlobalRegress_',Var1Name,'-RegressedOver-',Var2Name,'.png']);
hold off;

end