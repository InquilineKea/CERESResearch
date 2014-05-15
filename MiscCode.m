varfun = @(func, varargin)cellfun(...
    @(x)evalin('base', sprintf('%s(%s, y)', func, x)), ...
    arrayfun(@inputname, 2:nargin, 'UniformOutput', false), ...
    'UniformOutput', false);

varfun(myfun(Net,LAT))

for i=1:length(patient)
    sum(patient(i).test)
    varname(patient(i).name)
end

AllFluxes(1).name = 'blah';
AllFluxes(1).Flux = 5;

for i=1:length(AllFluxes)
    sum(patient(i).test)
    varname(patient(i).name)
end


a = 1:18
plot(xcorr(1:18,1:18))
plot(xcorr(1:18,-1:18))
plot(xcorr(1:18,rand(18,1)))
a = 1:18
b = [5 3 0 9 a];
plot(xcorr(a,b))
plot(xcorr(1:18,rand(18,1),'coeff'))
plot(xcorr(1:18,1:18,'coeff'))
plot(xcorr(sind(1:180),cosd(1:180),'coeff'))
plot(xcorr(sind(1:1000),cosd(1:1000),'coeff'))
plot(xcorr(sind(1:10000),cosd(1:10000),'coeff'))
plot(xcorr(sind(1:1000),cosd(1:1000),'coeff'))
plot(xcorr(sind(1:1000),sind(1:1000),'coeff'))
plot(xcorr(sind(1:1000),cosd(1:1000),'coeff'))
plot(xcorr(sind(1:1000),sind(1:1000),'coeff'))
plot(xcorr(sind(1:1000),cosd(1:1000),'coeff'))
plot(xcorr((1:1000)^2,cosd(1:1000),'coeff'))
plot(xcorr((1:1000).^2,cosd(1:1000),'coeff'))
plot(xcorr((1:1000).^2,(1:1000),'coeff'))
plot(xcorr((1:1000).^2,rand(1:1000),'coeff'))
plot(xcorr((1:1000).^2,rand(1000,1),'coeff'))
plot(xcorr(rand(1:1000),rand(1000,1),'coeff'))
plot(xcorr(rand(1000,1),rand(1000,1),'coeff'))
plot(xcorr(rand(1000,1),rand(1000,1)))