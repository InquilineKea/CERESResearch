

%%%%%%%%%%

ncdisp('tavg_vvel.nc')
lat=ncread('tavg_vvel.nc','Latitude_v');
lon=ncread('tavg_vvel.nc','Longitude_t');
dep=ncread('tavg_vvel.nc','Depth_c');
depth_bottom = ncread('tavg_vvel.nc','Depth_l');
time=ncread('tavg_vvel.nc','Time');
curltau = ncread('tavg_curltau.nc','curltau');
etan = ncread('tavg_etan.nc','etan');


diffZ = repmat(diff(depth_bottom),1,360);
diffZ = diffZ';

diffDepth = [10;diff(depth_bottom)];

AllDepths = repmat(diffDepth,1,360);
AllDepths = AllDepths';

[LON,LAT]=meshgrid(lon,lat);
contourf(LON,LAT,ncread('tavg_curltau.nc','curltau')');colorbar
contourf(LON,LAT,ncread('tavg_taux.nc','taux')');colorbar


LAT2USE = -32;
Lon2Use = lon > 27.5 & lon < 115.6;
Omega=7.29e-5; 
f=2*Omega*sin(LAT2USE*pi/180);
beta = (2*Omega*sind(LAT2USE+1)-2*Omega*sind(LAT2USE-1))/(2*1113194);
dx = abs(110000*cosd(LAT2USE));

% Heres where you select your latitude
ndx=find(lat==LAT2USE)

% Here's the critical loading step where you grab only 1 latitude
v=ncread('tavg_vvel.nc','v',[1,ndx,1,1],[length(lon),1,length(dep),length(time)]);
u=ncread('tavg_uvel.nc','u',[1,ndx,1,1],[length(lon),1,length(dep),length(time)]);
w=ncread('tavg_wvel.nc','w',[1,ndx,1,1],[length(lon),1,length(dep),length(time)]);
phihyd = ncread('tavg_phihyd.nc','phyhydid',[1,ndx,1,1],[length(lon),1,length(dep),length(time)]);
rhoanom = ncread('tavg_rhoanoma.nc','rhoanoma',[1,ndx,1,1],[length(lon),1,length(dep),length(time)]);

% It is easier to interpret rhoanom, if you subtract off the average of rhoanom from each z level.  For example (sorry for the approximate coding here):
% 

v=squeeze(v);

rhoanom2=rhoanom;
for zz=1:length(rhoanom(1,1,:))
  rhoanom2(:,1,zz)=rhoanom(:,1,zz)-nanmean(rhoanom(Lon2Use,1,zz));
end

phihyd2=phihyd;
for zz=1:length(phihyd(1,1,:))
  phihyd2(:,1,zz)=phihyd(:,1,zz)-nanmean(phihyd(Lon2Use,1,zz));
end

phihyd = squeeze(phihyd);
phihyd2 = squeeze(phihyd2);
rhoanom = squeeze(rhoanom);
rhoanom2 = squeeze(rhoanom2);
CurlTauLat = curltau(:,ndx);
etanLat = etan(:,ndx);
taux = ncread('tavg_taux.nc','taux');
taux = taux(:,ndx);
H = 1500; %%%THIS PROBABLY NEEDS TO BE FIXED, EVERYONE HAS DIFFERENT VALUES
SSHGeo = 9.81*H*diff(etanLat)/(f*dx); 
Sverdrup = CurlTauLat/(1035*beta);
EkmanLayerV = taux/(1035*f);
EkmanLayerVSverdrup = EkmanLayerV(Lon2Use)*dx/10^6;

DiffV = diff(v,1,2);
vflux = double(v.*AllDepths*111000*cosd(LAT2USE));
vDz = double(v.*AllDepths);
dVdZ = double(DiffV./diffZ);
% 
% plot(CurlTauLat,'r')
% hold on;
% plot(curltau(:,47))
% plot(curltau(:,55),'g')

% figure(1)
% pcolor(lon,-dep,rhoanom');shading('flat');colorbar
% figure(2)

% For the thermal wind equation, we only need d/dx(rhoanom), so subtracting out any function of z doesn't affect the result.  In this case, the function of z is the along-section average of rhoanom.
% 
% What's left over will show only the density anomailies relevant to the thermal wind equation.

% You might want to squeeze out the extra dimension in v


% Then you can plot
subplot(3,2,1)
set(gca,'FontSize',20)
pcolor(lon(Lon2Use),-dep,v(Lon2Use,:)');shading('flat');colorbar
title('v along constant latitude line')
xlabel('longitude')
ylabel('depth')
subplot (3,2,2)
pcolor(lon(Lon2Use),-dep,vflux(Lon2Use,:)'/10^6);shading('flat');colorbar
set(gca,'FontSize',20)
title('vflux (in Sverdrups?)')
xlabel('longitude')
ylabel('depth')
subplot(3,2,3)
pcolor(lon(Lon2Use),-dep,rhoanom2(Lon2Use,:)');shading('flat');colorbar
set(gca,'FontSize',20)
title('rhoanom along constant latitude line')
xlabel('longitude')
ylabel('depth')
subplot(3,2,4)
pcolor(lon(Lon2Use),-dep,phihyd2(Lon2Use,:)');shading('flat');colorbar
% caxis([-20 50])
%pcolor(lon(Lon2Use),-dep(dep < 1500),phihyd(Lon2Use,dep < 1500)');shading('flat');colorbar
set(gca,'FontSize',20)
title('phihyd along constant latitude line')
xlabel('longitude')
ylabel('depth')
MeanLayer = nanmean(dVdZ);
subplot(3,2,5)
pcolor(lon(Lon2Use),-cumsum(diff(dep)),dVdZ(Lon2Use,:)');shading('flat');colorbar
set(gca,'FontSize',20)
title('dV/dZ along constant latitude line')
xlabel('longitude')
ylabel('depth')
caxis([-MeanLayer(2) MeanLayer(2)])
str=sprintf('Plot with latitude = %d',LAT2USE);
[ax4,h3]=suplabel(str,'t');
set(h3,'FontSize',30)

%%%%%%%%%

subplot (2,2,1)
pcolor(lon(Lon2Use),-dep,vflux(Lon2Use,:)'/10^6);shading('flat');colorbar
set(gca,'FontSize',20)
title('Meridional Volume Transport (in Sv)')
xlabel('longitude')
ylabel('depth')
subplot(2,2,2)
pcolor(lon(Lon2Use),-dep,rhoanom2(Lon2Use,:)');shading('flat');colorbar
set(gca,'FontSize',20)
title('Density Anomalies (in kg/m^3)')
xlabel('longitude')
ylabel('depth')
subplot(2,2,3)
pcolor(lon(Lon2Use),-dep,phihyd2(Lon2Use,:)');shading('flat');colorbar
set(gca,'FontSize',20)
title('Hydrostatic Potential Anomalies (in m^2/s^2)')
xlabel('longitude')
ylabel('depth')
subplot(2,2,4)
plot(lon(Lon2Use),etanLat(Lon2Use,:)')
grid on;
set(gca,'FontSize',20)
title('Sea Surface Height Anomalies')
xlabel('longitude')
ylabel('Deviation from Sea Level (meters)')


%%%%%%%%%%


dep2 = cumsum(diff(dep));

pcolor(lon(Lon2Use),-dep(dep < 1500),phihyd(Lon2Use,dep < 1500)');shading('flat');colorbar

pcolor(lon(Lon2Use),-dep2(dep < 1500),dVdZ(Lon2Use,dep < 1500)');shading('flat');colorbar
pcolor(lon(Lon2Use),-dep2(dep < 500),dVdZ(Lon2Use,dep < 500)');shading('flat');colorbar


%%%%%%u & w
subplot(2,2,1)
pcolor(lon(Lon2Use),-dep,v(Lon2Use,:)');shading('flat');colorbar
subplot(2,2,2)
pcolor(lon(Lon2Use),-dep,u(Lon2Use,:)');shading('flat');colorbar
subplot(2,2,3)
pcolor(lon(Lon2Use),-dep,w(Lon2Use,:)');shading('flat');colorbar

%%%%all longs
subplot(2,2,1)
pcolor(lon(:),-dep,v(:,:)');shading('flat');colorbar
subplot(2,2,2)
pcolor(lon(:),-dep,u(:,:)');shading('flat');colorbar
subplot(2,2,3)
pcolor(lon(:),-dep,w(:,:)');shading('flat');colorbar


%%%%%%%%%%%%%%2D Volume Transports


%%%%%%%%%%%%%

vDz = double(v.*AllDepths);
PressureGeoV = diff(phihyd)/(dx*f);
ThermalShear = -9.81*diff(rhoanom2)/(dx*f*1035);
PressureGeoVDz= double(PressureGeoV.*AllDepths(1:end-1,:));
ThermalShearV = double(ThermalShear.*AllDepths(1:end-1,:));
ThermalShearV = [-PressureGeoV(1,:);ThermalShearV];
ThermalShearV2 = double(nancumsum(-ThermalShearV,2));
ThermalShearV2 = ThermalShearV2(1:end-1,:);
%GroundV = repmat(v(1,:),359,1);
ThermalShearV2(isnan(PressureGeoVDz)) = NaN;
%ThermalShearV2 = fliplr(ThermalShearV2);
ThermalShearV2Sverdrup = double(ThermalShearV2.*AllDepths(2:end,:)*111000*cosd(LAT2USE))/10^6;

subplot(2,1,1)
pcolor(lon(Lon2Use),-dep,ThermalShear(Lon2Use,:)');shading('flat');colorbar
title('thermal shear')
subplot(2,1,2)
pcolor(lon(Lon2Use),-dep,ThermalShearV2Sverdrup(Lon2Use,:)');shading('flat');colorbar
title('thermal shear sverdrup')

PressureGeoVflux = double(PressureGeoV.*AllDepths(1:end-1,:)*111000*cosd(LAT2USE));
ThermalShearVFlux = double(ThermalShearV.*AllDepths(1:end-1,:)*111000*cosd(LAT2USE));
ThermalShearVFlux2 = double(ThermalShearV2.*AllDepths(1:end-1,:)*111000*cosd(LAT2USE));

% dVdZ2 = double(nancumsum(fliplr(dVdZ),2).*fliplr(AllDepths(:,2:end)));
% dVdZ2(dVdZ2==0) = NaN;
% dVdZ2 = fliplr(dVdZ2);
% dVdZ2Sverdrup = double(dVdZ2.*AllDepths(:,2:end)*111000*cosd(LAT2USE))/10^6;
% 
% pcolor(lon(Lon2Use),-cumsum(diff(dep)),dVdZ2Sverdrup(Lon2Use,:)');shading('flat');colorbar
% 

PressureGeoSverdrup = nansum(PressureGeoVflux(Lon2Use,:),2)'/10^6;
SSHGeoSverdrup = 111000*cosd(LAT2USE)*SSHGeo(Lon2Use)/10^6;
WindSverdrup = 111000*cosd(LAT2USE)*Sverdrup(Lon2Use)/10^6;
ThermalWindSverdrup = nansum(ThermalShearV2Sverdrup(Lon2Use,:),2)';

%%%%%%%%%%%%

subplot(2,1,1)
pcolor(lon(Lon2Use),-dep,v(Lon2Use,:)');shading('flat');colorbar
%caxis([0 0.1])
subplot(2,1,2)
pcolor(lon(Lon2Use),-dep,PressureGeoV(Lon2Use,:)');shading('flat');colorbar
set(gca,'FontSize',20)
title('V meridional transport along constant latitude line')
xlabel('longitude')
ylabel('depth')

%%%%%

subplot(2,1,1)
pcolor(lon(Lon2Use),-dep,vflux(Lon2Use,:)'/10^6);shading('flat');colorbar
subplot(2,1,2)
pcolor(lon(Lon2Use),-dep,PressureGeoVflux(Lon2Use,:)'/10^6);shading('flat');colorbar
set(gca,'FontSize',20)
title('v along constant latitude line (in Sverdrups)')
xlabel('longitude')
ylabel('depth')

%%%%

subplot(2,1,1)
pcolor(lon(Lon2Use),-dep,vflux(Lon2Use,:)'/10^6);shading('flat');colorbar
subplot(2,1,2)
%pcolor(lon(Lon2Use),-dep,ThermalShearVFlux(Lon2Use,:)'/10^6);shading('flat');colorbar
pcolor(lon(Lon2Use),-dep,ThermalShearV2Sverdrup(Lon2Use,:)');shading('flat');colorbar
set(gca,'FontSize',20)
title('v along constant latitude line (in Sverdrups) for Thermal Shear')
xlabel('longitude')
ylabel('depth')
pcolor(lon(Lon2Use),-dep,ThermalShearVFlux(Lon2Use,:)'/10^6-vflux(Lon2Use,:)'/10^6);shading('flat');colorbar
caxis([-.5 0])

%is diff(dep) really the best thing to use? maybe try interpolating between
%the gridpoints instead...

%also, look at the v plots. dVdZ should be negative at a WBC since it
%becomes more negative as you go up in depth

subplot(2,1,1)
pcolor(lon(Lon2Use,:),-cumsum(diff(dep)),-dVdZ(Lon2Use,:)');shading('flat');colorbar
caxis([-5e-4 5e-4])
set(gca,'FontSize',20)
set(colorbar,'fontsize',20);
set(get(colorbar,'ylabel'),'string','dV/dZ (m/s)','FontSize',20)
title('dV/dZ, calculated from ECCO 3.73 v')
xlabel('longitude')
ylabel('depth')
subplot(2,1,2)
pcolor(lon(Lon2Use,:),-dep,ThermalShear(Lon2Use,:)');shading('flat');colorbar
caxis([-5e-4 5e-4])
set(gca,'FontSize',20)
set(get(colorbar,'ylabel'),'string','dV/dZ (m/s)','FontSize',20)
title('dV/dZ, calculated from ECCO 3.73 v')
%ylabel(colorbar, 'dV/dZ (m/s)')
title('estimated dV/dZ from thermal wind relation')
xlabel('longitude')
ylabel('depth')



% contourf(lon,-dep,vflux');shading('flat');colorbar
% 
% contourf(lon,-dep,rhoanom');shading('flat');colorbar
% 
% contourf(lon,-dep,phihyd');shading('flat');colorbar
% 



plot(lon(Lon2Use),etanLat(Lon2Use,:)')

%These are all Ekman transports, m^2/s

% plot(lon(Lon2Use),nansum(vDz(Lon2Use,:)'))
% grid on;
% hold on;
% plot(lon(Lon2Use),Sverdrup(Lon2Use),'g')%should NOT agree in the WBC
% plot(lon(Lon2Use)+0.5,SSHGeo(Lon2Use),'r')
% plot(lon(Lon2Use)+0.5,nansum(PressureGeoVDz(Lon2Use,:)'),'m')
% set(gca,'FontSize',20)
% str2=sprintf('Integrated vertical flux at latitude = %d',LAT2USE);
% xlabel('longitude')
% ylabel('Integrated vertical flux')
% title(str2)
% set(legend('Integrated vDz','Sverdrup, should not agree for WBC','SSHGeo','INtegrated vDz from phihyd relation'),'Location','NorthWestOutside')
% 
% sum(nansum(vDz'))

%%%ALL FLOWS IN SVERDRUP BELOW

VSverdrup = 111000*cosd(LAT2USE)*nansum(vDz(Lon2Use,:)'/10^6);
AllVSverdrupLayers = 111000*cosd(LAT2USE)*vDz(Lon2Use,:)/10^6;
VSverdrupByZ = nansum(AllVSverdrupLayers,1);
AbsVSverdrupByZ = nansum(abs(AllVSverdrupLayers),1);

subplot(2,1,1)
plot(dep,VSverdrupByZ,'-*','LineWidth',3)
hold on;
grid on;
plot(dep,nansum(PressureGeoVflux(Lon2Use,:))/10^6,'r','LineWidth',3)
%plot(dep,nansum(ThermalShearV2Sverdrup(Lon2Use,:)),'color',[0 .5 0],'LineWidth',3)
set(legend('ECCO 3.73-Calculated v', 'Hydrostatic Geostrophic Relation'),'Location','NorthWestOutside')
set(gca,'FontSize',20)
title('Zonally-Integrated Net Meridional Flows Categorized by Depth')
xlabel('Depth (m)')
ylabel('Volume Transport (Sv)')

subplot(2,1,2)

%%%%
plot(dep,AbsVSverdrupByZ,'-*','LineWidth',3)
hold on;
grid on;
plot(dep,nansum(abs(PressureGeoVflux(Lon2Use,:)))/10^6,'r','LineWidth',3)
%plot(dep,nansum(abs(ThermalShearV2Sverdrup(Lon2Use,:))),'color',[0 .5 0],'LineWidth',3)
%set(legend('ECCO 3.73-Calculated v', 'Hydrostatic Geostrophic Relation', 'Thermal Wind Relation'),'Location','NorthWestOutside')
set(gca,'FontSize',20)
title('Zonally-Integrated Absolute Value Meridional Flow Categorized by Depth')
xlabel('Depth (m)')
ylabel('Volume Transport (Sv)')

%%%%

sum(VSverdrupByZ)
sum(AbsVSverdrupByZ)

%OldThermalWindSverdrup = nansum(ThermalShearVFlux(Lon2Use,:),2)'/10^6;


% plot(lon(Lon2Use),EkmanLayerVSverdrup,'LineWidth',2)


plot(lon(Lon2Use),VSverdrup,'-*','LineWidth',3)
grid on;
hold on;
plot(lon(Lon2Use),WindSverdrup,'g','LineWidth',3)%should NOT agree in the WBC
plot(lon(Lon2Use)+0.5,SSHGeoSverdrup,'r','LineWidth',3)
plot(lon(Lon2Use)+0.5,PressureGeoSverdrup,'m','LineWidth',3)
plot(lon(Lon2Use)+0.5,ThermalWindSverdrup,'k','LineWidth',3)
%plot(lon(Lon2Use)+0.5,OldThermalWindSverdrup,'--x','LineWidth',2)
%plot(lon(Lon2Use)+0.5,SSHGeoSverdrup'+ThermalWindSverdrup,'k','LineWidth',2)
set(gca,'FontSize',20)
str2=sprintf('Estimated Meridional Volume Transports at latitude = %d',LAT2USE);
xlabel('longitude')
ylabel('Meridional Volume Transport (Sv)')
title(str2)
set(legend('ECCO 3.73 Vertically-Integrated Meridional Transport','Sverdrup relation','Sea Surface Height Geostrophic', 'Hydrostatic Geostrophic', 'Thermal Wind'),'Location','NorthWestOutside')
%set(legend('Integrated Meridional Volume Transport from ECCO 3.73','Geostrophic (from SSH)', 'Geostrophic (from Pressure)', 'ThermalWind','SSHGeo + Thermal Wind'),'Location','NorthWestOutside')
%set(legend('Integrated Meridional Volume Transport from ECCO 3.73','Geostrophic (from SSH)', 'Geostrophic (from Pressure)', 'ThermalWind'),'Location','NorthWestOutside')
ylim([-30,20])

semilogy(lon(Lon2Use),VSverdrup,'-*','LineWidth',3)
grid on;
hold on;
semilogy(lon(Lon2Use),WindSverdrup,'g','LineWidth',3)%should NOT agree in the WBC
semilogy(lon(Lon2Use)+0.5,SSHGeoSverdrup,'r','LineWidth',3)
semilogy(lon(Lon2Use)+0.5,PressureGeoSverdrup,'m','LineWidth',3)
semilogy(lon(Lon2Use)+0.5,ThermalWindSverdrup,'k','LineWidth',3)
set(gca,'FontSize',20)
str2=sprintf('Estimated Meridional Volume Transports at latitude = %d',LAT2USE);
xlabel('longitude')
ylabel('Meridional Volume Transport (Sv)')
title(str2)
set(legend('ECCO 3.73 Vertically-Integrated Meridional Transport','Sverdrup relation','Sea Surface Height Geostrophic', 'Hydrostatic Geostrophic', 'Thermal Wind'),'Location','NorthWestOutside')
%set(legend('Integrated Meridional Volume Transport from ECCO 3.73','Geostrophic (from SSH)', 'Geostrophic (from Pressure)', 'ThermalWind','SSHGeo + Thermal Wind'),'Location','NorthWestOutside')
%set(legend('Integrated Meridional Volume Transport from ECCO 3.73','Geostrophic (from SSH)', 'Geostrophic (from Pressure)', 'ThermalWind'),'Location','NorthWestOutside')


nancov(SSHGeoSverdrup'+ThermalWindSverdrup,PressureGeoSverdrup)
nancov(VSverdrup,PressureGeoSverdrup)
nancov(VSverdrup,SSHGeoSverdrup'+ThermalWindSverdrup)

LonEastofWBC = Lon2Use & lon > 40;


nancorrcoef(ThermalWindSverdrup(10:end-6),SSHGeoSverdrup(10:end-6)')
nancorrcoef(ThermalWindSverdrup(10:end-6),VSverdrup(10:end-6))
nancorrcoef(ThermalWindSverdrup(10:end-6),PressureGeoSverdrup(10:end-6))
nancorrcoef(SSHGeoSverdrup(10:end-6)',PressureGeoSverdrup(10:end-6))

nancorrcoef(SSHGeoSverdrup(10:end-6)',VSverdrup(10:end-6))
nancorrcoef(SSHGeoSverdrup(10:end-6)',PressureGeoSverdrup(10:end-6))
nancorrcoef(WindSverdrup(10:end-6)',VSverdrup(10:end-6))
nancorrcoef(VSverdrup(10:end-6),PressureGeoSverdrup(10:end-6))
nancorrcoef(SSHGeoSverdrup(10:end-6)'+ThermalWindSverdrup(10:end-6),PressureGeoSverdrup(10:end-6))
nancorrcoef(VSverdrup(10:end-6),SSHGeoSverdrup(10:end-6)'+ThermalWindSverdrup(10:end-6))

nancorrcoef(OldThermalWindSverdrup(10:end-6),VSverdrup(10:end-6))
nancorrcoef(OldThermalWindSverdrup(10:end-6),SSHGeoSverdrup(10:end-6)')
nancorrcoef(ThermalWindSverdrup(10:end-6),OldThermalWindSverdrup(10:end-6))

[sum(nansum(VSverdrup')) sum((abs(VSverdrup')))]
sum(nansum(VSverdrup'))/nansum((abs(VSverdrup')))
sum(nansum(WindSverdrup'))/nansum((abs(WindSverdrup')))
sum(nansum(SSHGeoSverdrup'))/nansum((abs(SSHGeoSverdrup')))
sum(nansum(PressureGeoSverdrup'))/nansum((abs(PressureGeoSverdrup')))
sum(nansum(ThermalWindSverdrup'))/nansum((abs(ThermalWindSverdrup')))
sum(nansum(ThermalWindSverdrup+SSHGeoSverdrup'))/nansum((abs(ThermalWindSverdrup+SSHGeoSverdrup')))

%%%TESTING SOME FEATURES BELOW

% plot(lon(Lon2Use),111000*cosd(LAT2USE)*vDz(Lon2Use,[1 5 9 13 17])'/10^6)
% set(legend('1','2','3','4','5','6','7'),'Location','NorthWestOutside')

subplot(3,1,1)

plot(lon(Lon2Use),nansum(AllVSverdrupLayers(:,:),2),'-*','LineWidth',2)
hold on;
grid on;
plot(lon(Lon2Use),nansum(AllVSverdrupLayers(:,depth_bottom < 1000),2),'c','LineWidth',3)
plot(lon(Lon2Use),nansum(AllVSverdrupLayers(:,depth_bottom > 1000 & depth_bottom < 2000),2),'color',[0 .2 0],'LineWidth',3)
plot(lon(Lon2Use),nansum(AllVSverdrupLayers(:,depth_bottom > 2000 & depth_bottom < 3000),2),'r','LineWidth',3)
plot(lon(Lon2Use),nansum(AllVSverdrupLayers(:,depth_bottom > 3000),2),'-xm','LineWidth',3)
set(gca,'FontSize',20)
str2=sprintf('Estimated Meridional Volume Transports at latitude = %d, calculated from ECCO 3.73 v values',LAT2USE);
xlabel('longitude')
ylabel('Meridional Volume Transport (Sv)')
title(str2)
set(legend('all depths','0-985m','985-1750m','1750-2700m','below 2700m'),'Location','NorthWestOutside')

BottomSverdrup = nansum(AllVSverdrupLayers(:,depth_bottom > 3000),2);
nancorrcoef(BottomSverdrup(10:end-6)',VSverdrup(10:end-6))

%%%%%%%%%%

%%%%%%%%%%%
hold off;

subplot(3,1,2)

plot(lon(Lon2Use),111000*cosd(LAT2USE)*nansum(PressureGeoVDz(Lon2Use,:)'/10^6))
hold on;
grid on;
plot(lon(Lon2Use),111000*cosd(LAT2USE)*nansum(PressureGeoVDz(Lon2Use,1:6)'/10^6),'c')
plot(lon(Lon2Use),111000*cosd(LAT2USE)*nansum(PressureGeoVDz(Lon2Use,7:11)'/10^6),'g')
plot(lon(Lon2Use),111000*cosd(LAT2USE)*nansum(PressureGeoVDz(Lon2Use,12:17)'/10^6),'r')
plot(lon(Lon2Use),111000*cosd(LAT2USE)*nansum(PressureGeoVDz(Lon2Use,17:22)'/10^6),'m')
set(legend('Total','0-100','100-510','510-2700','2700-5200'),'Location','NorthWestOutside')

hold off;

subplot(3,1,3)


%%%%%

subplot(2,1,1)
plot(lon(Lon2Use),nanmean(dVdZ(Lon2Use,1:6)'))
hold on;
grid on;
plot(lon(Lon2Use),nanmean(dVdZ(Lon2Use,7:11)'),'c')
plot(lon(Lon2Use),nanmean(dVdZ(Lon2Use,12:17)'),'g')
plot(lon(Lon2Use),nanmean(dVdZ(Lon2Use,17:22)'),'r')
set(legend('0-100','100-510','510-2700','2700-5200'),'Location','NorthWestOutside')
title('Thermal Shear')

subplot(2,1,2)
plot(lon(Lon2Use),nanmean(ThermalShear(Lon2Use,1:6)'))
hold on;
grid on;
plot(lon(Lon2Use),nanmean(ThermalShear(Lon2Use,7:11)'),'c')
plot(lon(Lon2Use),nanmean(ThermalShear(Lon2Use,12:17)'),'g')
plot(lon(Lon2Use),nanmean(ThermalShear(Lon2Use,17:22)'),'r')
set(legend('0-100','100-510','510-2700','2700-5200'),'Location','NorthWestOutside')

title('Thermal Shear from density anomalies')

%%%%%%%%%%%%% rhoanoms - are they governed by SSH?

hold off;
subplot(2,1,1)
plot(lon(Lon2Use),nanmean(rhoanom2(Lon2Use,1:6)'),'m')
hold on;
grid on;
plot(lon(Lon2Use),nanmean(rhoanom2(Lon2Use,7:11)'),'c')
%plot(lon(Lon2Use),nanmean(rhoanom2(Lon2Use,17:22)'),'--r')
[ax,h1,h2] = plotyy(lon(Lon2Use),nanmean(rhoanom2(Lon2Use,12:13)'),lon(Lon2Use)+0.5, etanLat(Lon2Use,:)','plot')
set(legend('0-100','100-510','SSH'),'Location','NorthOutside')
set(h2,'marker','s','color','green')
set(get(ax(2),'ylabel'), 'String', 'New Left Label')
set(get(ax(2),'ylabel'), 'color', 'green')
set(ax(2), 'YColor', 'green')
set(ax(1), 'YColor', 'k')
%axis([30 120 -2 2])
axis(ax(1), [29 116 -2 2])
axis(ax(2), [29 116 -1 1])

plot(diff(etanLat(Lon2Use,:)))

subplot(2,1,2)

pcolor(lon(Lon2Use),-dep,phihyd2(Lon2Use,:)');shading('flat');colorbar('horiz')
set(gca,'FontSize',20)
title('phihyd along constant latitude line')
xlabel('longitude')
ylabel('depth')



pcolor(lon(Lon2Use),-dep,rhoanom2(Lon2Use,:)');shading('flat');colorbar
set(gca,'FontSize',20)
title('rhoanom along constant latitude line')
xlabel('longitude')
ylabel('depth')



plot(lon(Lon2Use),nanmean(phihyd2(Lon2Use,1:6)'))
hold on;
grid on;
plot(lon(Lon2Use),nanmean(phihyd2(Lon2Use,7:11)'),'c')
plot(lon(Lon2Use),nanmean(phihyd2(Lon2Use,12:17)'),'g')
plot(lon(Lon2Use),nanmean(phihyd2(Lon2Use,17:22)'),'r')
%plot(lon(Lon2Use),nanmean(phihyd(Lon2Use,:)'))


%%%%%%%%

subplot(2,2,1)
plot(lon(Lon2Use),nanmean(rhoanom2(Lon2Use,12:13)'))
subplot(2,2,2)
plot(lon(Lon2Use),111000*cosd(LAT2USE)*nansum(PressureGeoVDz(Lon2Use,12:13)'/10^6))
subplot(2,2,3)
plot(lon(Lon2Use),111000*cosd(LAT2USE)*nansum(vDz(Lon2Use,12:13)'/10^6))

subplot(2,2,1)
plot(lon(Lon2Use),rhoanom2(Lon2Use,16)')
title('rhoanom')
subplot(2,2,2)
plot(lon(Lon2Use),111000*cosd(LAT2USE)*PressureGeoVDz(Lon2Use,16)'/10^6)
title('Pressure Geo v')
subplot(2,2,3)
plot(lon(Lon2Use),111000*cosd(LAT2USE)*vDz(Lon2Use,16)'/10^6)
title('actual v')

% VSverdrup = 111000*cosd(LAT2USE)*nansum(vDz(Lon2Use,:)'/10^6);
% PressureGeoSverdrup = nansum(PressureGeoVflux(Lon2Use,:),2)'/10^6;
% SSHGeoSverdrup = 111000*cosd(LAT2USE)*SSHGeo(Lon2Use)/10^6;
% WindSverdrup = 111000*cosd(LAT2USE)*Sverdrup(Lon2Use)/10^6;
% ThermalWindSverdrup = nansum(ThermalShearVFlux(Lon2Use,:),2)'/10^6;

%%%%%%CUMSUM BELOW

plot(lon(Lon2Use),nancumsum(VSverdrup),'-*','LineWidth',2)
grid on;
hold on;
%plot(lon(Lon2Use),nancumsum(WindSverdrup),'g','LineWidth',2)%should NOT agree in the WBC
plot(lon(Lon2Use)+0.5,nancumsum(SSHGeoSverdrup),'r','LineWidth',3)
plot(lon(Lon2Use)+0.5,cumsum(PressureGeoSverdrup),'m','LineWidth',3)
plot(lon(Lon2Use)+0.5,nancumsum(ThermalWindSverdrup),'color',[0 .5 0],'LineWidth',3)
set(gca,'FontSize',20)
str2=sprintf('Cumulative Meridional Volume Transport at latitude = %d',LAT2USE);
xlabel('longitude')
ylabel('Cumulative Volume Transport (Sv) added up to a given longitude')
title(str2)
set(legend('Integrated Meridional Volume Transport from ECCO 3.73 Velocity','Sea Surface Height Geostrophic','Hydrostatic Geostrophic', 'Thermal Wind'),'Location','NorthWestOutside')



cumsum(111000*cosd(LAT2USE)*nansum(vDz(Lon2Use,:)'/10^6))

%this is the total volume of the WBC?
min(cumsum(111000*cosd(LAT2USE)*nansum(vDz(Lon2Use,:)'/10^6)))

%compare net Sverdrup transport with total Sverdrup transport

sum(111000*cosd(LAT2USE)*nansum(vDz(Lon2Use,:)'/10^6))
sum(abs(111000*cosd(LAT2USE)*nansum(vDz(Lon2Use,:)'/10^6)))

nansum(111000*cosd(LAT2USE)*Sverdrup(Lon2Use)/10^6)
nansum(abs(111000*cosd(LAT2USE)*Sverdrup(Lon2Use)/10^6))

nansum(111000*cosd(LAT2USE)*SSHGeo(Lon2Use)/10^6)
nansum(abs(111000*cosd(LAT2USE)*SSHGeo(Lon2Use)/10^6))

plotyy(lon(Lon2Use)+0.5,111000*cosd(LAT2USE)*nansum(PressureGeoVDz(Lon2Use,:)'/10^6),lon(Lon2Use)+0.5, etanLat(Lon2Use,:)','plot') 

plotyy(lon(Lon2Use),nanmean(rhoanom2(Lon2Use,:)'),lon(Lon2Use)+0.5, etanLat(Lon2Use,:)','plot') 


plotyy(lon,nansum(vDz'),lon, Sverdrup,'plot') %SSH Geostrophic, should agree everywhere


%    [ax1,h1]=suplabel('super X label');
%    [ax2,h2]=suplabel('super Y label','y');
%    [ax3,h2]=suplabel('super Y label (right)','yy');


text(0.5, 1,'\bf Do you like this title?','HorizontalAlignment' ,'center','VerticalAlignment', 'top')

