
ensoTime = ncread('enso.cdf','T');
NINO34 = ncread('enso.cdf','NINO34')
NINO34 = NINO34(find(ensoTime > 228 & ensoTime < 638));

temp = ncread('t2m.nc','t2m');
temptime = ncread('t2m.nc','time');
tempLats= ncread('t2m.nc','latitude');
tempLongs = ncread('t2m.nc','longitude');
temp = temp(:,:,1:find(temptime==991296));
temp = permute(temp,[2 1 3]);
E=zeros(180,360,156);
for depth=1:size(temp,3)
  temp(:,:,depth) = flipud(temp(:,:,depth));
  E(:,:,depth)=imresize(temp(:,:,depth),[180 360]);
end
RawFlux.Temp = E;
clear E;
