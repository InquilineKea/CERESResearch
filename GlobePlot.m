function [test] = GlobePlot(Flux,Month)

load geoid; ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true); geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
if Month==0
    geoshow(mean(Flux,3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
else
    geoshow(mean(Flux(:,:,Month:12:end),3), geoidrefvec, 'DisplayType', 'texturemap');colorbar
end

test = 1;


end