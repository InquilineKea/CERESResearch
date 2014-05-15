function GlobalCorrs = GlobalCorrMap(Flux1,Flux2)
lats = min(size(Flux1,1),size(Flux2,1));
longs = min(size(Flux1,2),size(Flux2,2));
time=size(Flux1,3);
Factor = 180/lats;

if size(Flux1,1) < size(Flux2,1)
    Matrix2Resize = Flux2;
    OtherMatrix = Flux1;
elseif size(Flux2,1) < size(Flux1,1)
    Matrix2Resize = Flux1;
    OtherMatrix = Flux2;
else
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

%take correlation across ALL time here...
GlobalCorrs = bsxfun(@rdivide,sum(Flux1.*Flux2,3),(sqrt(sum(Flux1.^2,3)).*sqrt(sum(Flux2.^2,3))));
load LatWeights.mat
CorrMean = mean(mean(resizem(GlobalCorrs,Factor),2).*LatWeights(:,2));

load geoid
set(gca,'FontSize',20)
ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(resizem(GlobalCorrs,Factor),geoidrefvec,'DisplayType','texturemap');colorbar
caxis([-1 1])

title(['Correlation between ',inputname(1),' and ',inputname(2), '. Mean = ', num2str(CorrMean)])
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['GlobalCorr_',num2str(inputname(1)),'-',num2str(inputname(2)),'.png']);
hold off;

end