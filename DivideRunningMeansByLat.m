NHMinusSH(Net,time,30,90)
title('30deg-90deg Running Mean in net')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','net30_90.png']);
NHMinusSH(Net,time,0,30)
title('0-30deg Running Mean in net')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','net0_30.png']);

NHMinusSH(LWCF,time,30,90)
title('30deg-90deg Running Mean in LWCF')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LWCF30_90.png']);
NHMinusSH(LWCF,time,0,30)
title('0-30deg Running Mean in LWCF')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LWCF0_30.png']);

NHMinusSH(SWCF,time,30,90)
title('30deg-90deg Running Mean in SWCF')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SWCF30_90.png']);
NHMinusSH(SWCF,time,0,30)
title('0-30deg Running Mean in SWCF')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SWCF0_30.png']);

NHMinusSH(LW,time,30,90)
title('30deg-90deg Running Mean in LW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LW30_90.png']);
NHMinusSH(Net,time,0,30)
title('0-30deg Running Mean in LW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LW0_30.png']);

NHMinusSH(LWclear,time,30,90)
title('30deg-90deg Running Mean in LWClear')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LWClear30_90.png']);
NHMinusSH(LWclear,time,0,30)
title('0-30deg Running Mean in LWClear')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LWClear0_30.png']);

NHMinusSH(SWclear,time,30,90)
title('30deg-90deg Running Mean in SWclear')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SWclear30_90.png']);
NHMinusSH(SWclear,time,0,30)
title('0-30deg Running Mean in SWclear')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SWclear0_30.png']);

NHMinusSH(SW,time,30,90)
title('30deg-90deg Running Mean in SW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SW30_90.png']);
NHMinusSH(SW,time,0,30)
title('0-30deg Running Mean in SW')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SW0_30.png']);

NHMinusSH(netclear,time,30,90)
title('30deg-90deg Running Mean in netclear')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','netclear30_90.png']);
NHMinusSH(netclear,time,0,30)
title('0-30deg Running Mean in netclear')
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','netclear0_30.png']);

NHMinusSH(netclear,time,0,30)
NHMinusSH(SWCF,time,0,30)
cellfun(@(x) NHMinusSH(x,time,0,30),{SWCF,netclear},'UniformOutput',false);

cellfun(@(x) AllLatitudeMWA(x),{Net,netclear,SWCF},'UniformOutput',false);

AllLatitudeMWA(Net)
title('Moving Weighted Average Net Flux per Latitude Band, after mean-subtraction (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','net_All_Lat.png']);
hold off;
AllLatitudeMWA(netclear)
title('Moving Weighted Average Net CS Flux per Latitude Band, after mean-subtraction (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','netclear_All_Lat.png']);
hold off;
AllLatitudeMWA(SWCF)
title('Moving Weighted Average SWCF per Latitude Band, after mean-subtraction (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SWCF_AllLat.png']);
hold off;
AllLatitudeMWA(LWCF)
title('Moving Weighted Average LWCF per Latitude Band, after mean-subtraction (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LWCF_AllLat.png']);
hold off;
AllLatitudeMWA(SW)
title('Moving Weighted Average SW per Latitude Band, after mean-subtraction (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SW_AllLat.png']);
hold off;
AllLatitudeMWA(LW)
title('Moving Weighted Average LW per Latitude Band, after mean-subtraction (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LW_AllLat.png']);
hold off;
AllLatitudeMWA(SWclear)
title('Moving Weighted Average SWclear per Latitude Band, after mean-subtraction (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','SWclear_AllLat.png']);
hold off;
AllLatitudeMWA(LWclear)
title('Moving Weighted Average LWClear per Latitude Band, after mean-subtraction (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','LWclear_AllLat.png']);
hold off;
AllLatitudeMWA(TotalCloudForcing)
title('Moving Weighted Average TotalCloudForcing per Latitude Band, after mean-subtraction (Watts)')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['RunningMean_','TotalCloudForcing_AllLat.png']);


%land only!
AllLatitudeMWA(bsxfun(@times, Net,topo>=0))

%%%%

NetMWA = AllTileMWA(Net,'Net');
SWMWA = AllTileMWA(SW,'SW');
LWMWA = AllTileMWA(LW,'LW');
netclearMWA = AllTileMWA(netclear,'netclear');
SWCFMWA = AllTileMWA(SWCF,'SWCF');
LWCFMWA =AllTileMWA(LWCF,'LWCF');
SWclearMWA = AllTileMWA(SWclear,'SWclear');
LWclearMWA = AllTileMWA(LWclear,'LWclear');
TotalCloudForcingMWA = AllTileMWA(TotalCloudForcing,'TotalCloudForcing');

load MWAs.mat
cmap = jet;

[NetNHMWAAllLats, NetSHMWAAllLats, NetNHMinusSHMWAAllLats] = HemisphericFluxesClimSubtracted(NetMWA,0,90,'Net');
[Net0_15N_MWA, Net0_15S_MWA, Net0_15Dif_MWA] =PlotMWA(NetMWA,0,14.5,'Net');
[Net15_30N_MWA, Net15_30S_MWA, Net15_30Dif_MWA] = PlotMWA(NetMWA,14.5,30,'Net');
[Net30_49N_MWA, Net30_49S_MWA, Net30_49Dif_MWA] = PlotMWA(NetMWA,30,48.5,'Net');
[Net49_90N_MWA, Net49_90S_MWA, Net49_90Dif_MWA] = PlotMWA(NetMWA,48.5,90,'Net');

%plot([Net49_90S_MWA' Net30_49S_MWA' Net15_30S_MWA' Net0_15S_MWA' Net0_15N_MWA' Net15_30N_MWA'  Net30_49N_MWA' Net49_90N_MWA'])
NetMWA2Plot = [detrend(Net49_90S_MWA',0) detrend(Net30_49S_MWA',0) detrend(Net15_30S_MWA',0) detrend(Net0_15S_MWA',0) detrend(Net0_15N_MWA',0) detrend(Net15_30N_MWA',0)  detrend(Net30_49N_MWA',0) detrend(Net49_90N_MWA',0)];
for i=1:8
   hold on;
   plot(NetMWA2Plot(:,i),'color',cmap(8*i,:),'LineWidth',3)
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'FontSize',20)
title('Mean-Subtracted Net MWA Fluxes for Various Latitude Bands')
xlabel('Year end')
ylabel('Mean-Subtracted Flux (Watts/m^2)')
grid on;
set(legend(['90S-49S'],['49S-30S'],['30S-15S'],['15S-0S'],['0N-15N'],['15N-30N'],['30N-49N'],['49N-90N']))
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['Net_MWA_Flux_All_Lats','.png']);
hold off;

[SWNHMWAAllLats, SWSHMWAAllLats, SWNHMinusSHMWAAllLats] = HemisphericFluxesClimSubtracted(SWMWA,0,90,'SW');
[SW0_15N_MWA, SW0_15S_MWA, SW0_15Dif_MWA] =PlotMWA(SWMWA,0,14.5,'SW');
[SW15_30N_MWA, SW15_30S_MWA, SW15_30Dif_MWA] = PlotMWA(SWMWA,14.5,30,'SW');
[SW30_49N_MWA, SW30_49S_MWA, SW30_49Dif_MWA] = PlotMWA(SWMWA,30,48.5,'SW');
[SW49_90N_MWA, SW49_90S_MWA, SW49_90Dif_MWA] = PlotMWA(SWMWA,48.5,90,'SW');

%plot([SW49_90S_MWA' SW30_49S_MWA' SW15_30S_MWA' SW0_15S_MWA' SW0_15N_MWA' SW15_30N_MWA'  SW30_49N_MWA' SW49_90N_MWA'])
SWMWA2Plot = [detrend(SW49_90S_MWA',0) detrend(SW30_49S_MWA',0) detrend(SW15_30S_MWA',0) detrend(SW0_15S_MWA',0) detrend(SW0_15N_MWA',0) detrend(SW15_30N_MWA',0)  detrend(SW30_49N_MWA',0) detrend(SW49_90N_MWA',0)];
for i=1:8
   hold on;
   plot(SWMWA2Plot(:,i),'color',cmap(8*i,:),'LineWidth',3)
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'FontSize',20)
title('Mean-Subtracted SW MWA Fluxes for Various Latitude Bands')
xlabel('Year end')
ylabel('Mean-Subtracted Flux (Watts/m^2)')
grid on;
set(legend(['90S-49S'],['49S-30S'],['30S-15S'],['15S-0S'],['0N-15N'],['15N-30N'],['30N-49N'],['49N-90N']))
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['SW_MWA_Flux_All_Lats','.png']);
hold off;

[LWNHMWAAllLats, LWSHMWAAllLats, LWNHMinusSHMWAAllLats] = HemisphericFluxesClimSubtracted(LWMWA,0,90,'LW');
[LW0_15N_MWA, LW0_15S_MWA, LW0_15Dif_MWA] =PlotMWA(LWMWA,0,14.5,'LW');
[LW15_30N_MWA, LW15_30S_MWA, LW15_30Dif_MWA] = PlotMWA(LWMWA,14.5,30,'LW');
[LW30_49N_MWA, LW30_49S_MWA, LW30_49Dif_MWA] = PlotMWA(LWMWA,30,48.5,'LW');
[LW49_90N_MWA, LW49_90S_MWA, LW49_90Dif_MWA] = PlotMWA(LWMWA,48.5,90,'LW');

%plot([LW49_90S_MWA' LW30_49S_MWA' LW15_30S_MWA' LW0_15S_MWA' LW0_15N_MWA' LW15_30N_MWA'  LW30_49N_MWA' LW49_90N_MWA'])
LWMWA2Plot = [detrend(LW49_90S_MWA',0) detrend(LW30_49S_MWA',0) detrend(LW15_30S_MWA',0) detrend(LW0_15S_MWA',0) detrend(LW0_15N_MWA',0) detrend(LW15_30N_MWA',0)  detrend(LW30_49N_MWA',0) detrend(LW49_90N_MWA',0)];
for i=1:8
   hold on;
   plot(LWMWA2Plot(:,i),'color',cmap(8*i,:),'LineWidth',3)
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'FontSize',20)
title('Mean-Subtracted LW MWA Fluxes for Various Latitude Bands')
xlabel('Year end')
ylabel('Mean-Subtracted Flux (Watts/m^2)')
grid on;
set(legend(['90S-49S'],['49S-30S'],['30S-15S'],['15S-0S'],['0N-15N'],['15N-30N'],['30N-49N'],['49N-90N']))
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['LW_MWA_Flux_All_Lats','.png']);
hold off;

[SWCFNHMWAAllLats, SWCFSHMWAAllLats, SWCFNHMinusSHMWAAllLats] = HemisphericFluxesClimSubtracted(SWCFMWA,0,90,'SWCF');
[SWCF0_15N_MWA, SWCF0_15S_MWA, SWCF0_15Dif_MWA] =PlotMWA(SWCFMWA,0,14.5,'SWCF');
[SWCF15_30N_MWA, SWCF15_30S_MWA, SWCF15_30Dif_MWA] = PlotMWA(SWCFMWA,14.5,30,'SWCF');
[SWCF30_49N_MWA, SWCF30_49S_MWA, SWCF30_49Dif_MWA] = PlotMWA(SWCFMWA,30,48.5,'SWCF');
[SWCF49_90N_MWA, SWCF49_90S_MWA, SWCF49_90Dif_MWA] = PlotMWA(SWCFMWA,48.5,90,'SWCF');

%plot([SWCF49_90S_MWA' SWCF30_49S_MWA' SWCF15_30S_MWA' SWCF0_15S_MWA' SWCF0_15N_MWA' SWCF15_30N_MWA'  SWCF30_49N_MWA' SWCF49_90N_MWA'])
SWCFMWA2Plot = [detrend(SWCF49_90S_MWA',0) detrend(SWCF30_49S_MWA',0) detrend(SWCF15_30S_MWA',0) detrend(SWCF0_15S_MWA',0) detrend(SWCF0_15N_MWA',0) detrend(SWCF15_30N_MWA',0)  detrend(SWCF30_49N_MWA',0) detrend(SWCF49_90N_MWA',0)];
for i=1:8
   hold on;
   plot(SWCFMWA2Plot(:,i),'color',cmap(8*i,:),'LineWidth',3)
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'FontSize',20)
title('Mean-Subtracted SWCF MWA Fluxes for Various Latitude Bands')
xlabel('Year end')
ylabel('Mean-Subtracted Flux (Watts/m^2)')
grid on;
set(legend(['90S-49S'],['49S-30S'],['30S-15S'],['15S-0S'],['0N-15N'],['15N-30N'],['30N-49N'],['49N-90N']))
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['SWCF_MWA_Flux_All_Lats','.png']);
hold off;

[SWclearNHMWAAllLats, SWclearSHMWAAllLats, SWclearNHMinusSHMWAAllLats] = HemisphericFluxesClimSubtracted(SWclearMWA,0,90,'SWclear');
[SWclear0_15N_MWA, SWclear0_15S_MWA, SWclear0_15Dif_MWA] =PlotMWA(SWclearMWA,0,14.5,'SWclear');
[SWclear15_30N_MWA, SWclear15_30S_MWA, SWclear15_30Dif_MWA] = PlotMWA(SWclearMWA,14.5,30,'SWclear');
[SWclear30_49N_MWA, SWclear30_49S_MWA, SWclear30_49Dif_MWA] = PlotMWA(SWclearMWA,30,48.5,'SWclear');
[SWclear49_90N_MWA, SWclear49_90S_MWA, SWclear49_90Dif_MWA] = PlotMWA(SWclearMWA,48.5,90,'SWclear');

%plot([SWclear49_90S_MWA' SWclear30_49S_MWA' SWclear15_30S_MWA' SWclear0_15S_MWA' SWclear0_15N_MWA' SWclear15_30N_MWA'  SWclear30_49N_MWA' SWclear49_90N_MWA'])
SWclearMWA2Plot = [detrend(SWclear49_90S_MWA',0) detrend(SWclear30_49S_MWA',0) detrend(SWclear15_30S_MWA',0) detrend(SWclear0_15S_MWA',0) detrend(SWclear0_15N_MWA',0) detrend(SWclear15_30N_MWA',0)  detrend(SWclear30_49N_MWA',0) detrend(SWclear49_90N_MWA',0)];
for i=1:8
   hold on;
   plot(SWclearMWA2Plot(:,i),'color',cmap(8*i,:),'LineWidth',3)
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'FontSize',20)
title('Mean-Subtracted SWclear MWA Fluxes for Various Latitude Bands')
xlabel('Year end')
ylabel('Mean-Subtracted Flux (Watts/m^2)')
grid on;
set(legend(['90S-49S'],['49S-30S'],['30S-15S'],['15S-0S'],['0N-15N'],['15N-30N'],['30N-49N'],['49N-90N']))
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['SWclear_MWA_Flux_All_Lats','.png']);
hold off;

[LWCFNHMWAAllLats, LWCFSHMWAAllLats, LWCFNHMinusSHMWAAllLats] = HemisphericFluxesClimSubtracted(LWCFMWA,0,90,'LWCF');
[LWCF0_15N_MWA, LWCF0_15S_MWA, LWCF0_15Dif_MWA] =PlotMWA(LWCFMWA,0,14.5,'LWCF');
[LWCF15_30N_MWA, LWCF15_30S_MWA, LWCF15_30Dif_MWA] = PlotMWA(LWCFMWA,14.5,30,'LWCF');
[LWCF30_49N_MWA, LWCF30_49S_MWA, LWCF30_49Dif_MWA] = PlotMWA(LWCFMWA,30,48.5,'LWCF');
[LWCF49_90N_MWA, LWCF49_90S_MWA, LWCF49_90Dif_MWA] = PlotMWA(LWCFMWA,48.5,90,'LWCF');

%plot([LWCF49_90S_MWA' LWCF30_49S_MWA' LWCF15_30S_MWA' LWCF0_15S_MWA' LWCF0_15N_MWA' LWCF15_30N_MWA'  LWCF30_49N_MWA' LWCF49_90N_MWA'])
LWCFMWA2Plot = [detrend(LWCF49_90S_MWA',0) detrend(LWCF30_49S_MWA',0) detrend(LWCF15_30S_MWA',0) detrend(LWCF0_15S_MWA',0) detrend(LWCF0_15N_MWA',0) detrend(LWCF15_30N_MWA',0)  detrend(LWCF30_49N_MWA',0) detrend(LWCF49_90N_MWA',0)];
for i=1:8
   hold on;
   plot(LWCFMWA2Plot(:,i),'color',cmap(8*i,:),'LineWidth',3)
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'FontSize',20)
title('Mean-Subtracted LWCF MWA Fluxes for Various Latitude Bands')
xlabel('Year end')
ylabel('Mean-Subtracted Flux (Watts/m^2)')
grid on;
set(legend(['90S-49S'],['49S-30S'],['30S-15S'],['15S-0S'],['0N-15N'],['15N-30N'],['30N-49N'],['49N-90N']))
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['LWCF_MWA_Flux_All_Lats','.png']);
hold off;

[LWclearNHMWAAllLats, LWclearSHMWAAllLats, LWclearNHMinusSHMWAAllLats] = HemisphericFluxesClimSubtracted(LWclearMWA,0,90,'LWclear');
[LWclear0_15N_MWA, LWclear0_15S_MWA, LWclear0_15Dif_MWA] =PlotMWA(LWclearMWA,0,14.5,'LWclear');
[LWclear15_30N_MWA, LWclear15_30S_MWA, LWclear15_30Dif_MWA] = PlotMWA(LWclearMWA,14.5,30,'LWclear');
[LWclear30_49N_MWA, LWclear30_49S_MWA, LWclear30_49Dif_MWA] = PlotMWA(LWclearMWA,30,48.5,'LWclear');
[LWclear49_90N_MWA, LWclear49_90S_MWA, LWclear49_90Dif_MWA] = PlotMWA(LWclearMWA,48.5,90,'LWclear');

%plot([LWclear49_90S_MWA' LWclear30_49S_MWA' LWclear15_30S_MWA' LWclear0_15S_MWA' LWclear0_15N_MWA' LWclear15_30N_MWA'  LWclear30_49N_MWA' LWclear49_90N_MWA'])
LWclearMWA2Plot = [detrend(LWclear49_90S_MWA',0) detrend(LWclear30_49S_MWA',0) detrend(LWclear15_30S_MWA',0) detrend(LWclear0_15S_MWA',0) detrend(LWclear0_15N_MWA',0) detrend(LWclear15_30N_MWA',0)  detrend(LWclear30_49N_MWA',0) detrend(LWclear49_90N_MWA',0)];
for i=1:8
   hold on;
   plot(LWclearMWA2Plot(:,i),'color',cmap(8*i,:),'LineWidth',3)
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'FontSize',20)
title('Mean-Subtracted LWclear MWA Fluxes for Various Latitude Bands')
xlabel('Year end')
ylabel('Mean-Subtracted Flux (Watts/m^2)')
grid on;
set(legend(['90S-49S'],['49S-30S'],['30S-15S'],['15S-0S'],['0N-15N'],['15N-30N'],['30N-49N'],['49N-90N']))
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['LWclear_MWA_Flux_All_Lats','.png']);
hold off;

[netclearNHMWAAllLats, netclearSHMWAAllLats, netclearNHMinusSHMWAAllLats] = HemisphericFluxesClimSubtracted(netclearMWA,0,90,'netclear');
[netclear0_15N_MWA, netclear0_15S_MWA, netclear0_15Dif_MWA] =PlotMWA(netclearMWA,0,14.5,'netclear');
[netclear15_30N_MWA, netclear15_30S_MWA, netclear15_30Dif_MWA] = PlotMWA(netclearMWA,14.5,30,'netclear');
[netclear30_49N_MWA, netclear30_49S_MWA, netclear30_49Dif_MWA] = PlotMWA(netclearMWA,30,48.5,'netclear');
[netclear49_90N_MWA, netclear49_90S_MWA, netclear49_90Dif_MWA] = PlotMWA(netclearMWA,48.5,90,'netclear');

%plot([netclear49_90S_MWA' netclear30_49S_MWA' netclear15_30S_MWA' netclear0_15S_MWA' netclear0_15N_MWA' netclear15_30N_MWA'  netclear30_49N_MWA' netclear49_90N_MWA'])
netclearMWA2Plot = [detrend(netclear49_90S_MWA',0) detrend(netclear30_49S_MWA',0) detrend(netclear15_30S_MWA',0) detrend(netclear0_15S_MWA',0) detrend(netclear0_15N_MWA',0) detrend(netclear15_30N_MWA',0)  detrend(netclear30_49N_MWA',0) detrend(netclear49_90N_MWA',0)];
for i=1:8
   hold on;
   plot(netclearMWA2Plot(:,i),'color',cmap(8*i,:),'LineWidth',3)
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'FontSize',20)
title('Mean-Subtracted netclear MWA Fluxes for Various Latitude Bands')
xlabel('Year end')
ylabel('Mean-Subtracted Flux (Watts/m^2)')
grid on;
set(legend(['90S-49S'],['49S-30S'],['30S-15S'],['15S-0S'],['0N-15N'],['15N-30N'],['30N-49N'],['49N-90N']))
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['netclear_MWA_Flux_All_Lats','.png']);
hold off;

[TotalCloudForcingNHMWAAllLats, TotalCloudForcingSHMWAAllLats, TotalCloudForcingNHMinusSHMWAAllLats] = HemisphericFluxesClimSubtracted(TotalCloudForcingMWA,0,90,'TotalCloudForcing');
[TotalCloudForcing0_15N_MWA, TotalCloudForcing0_15S_MWA, TotalCloudForcing0_15Dif_MWA] =PlotMWA(TotalCloudForcingMWA,0,14.5,'TotalCloudForcing');
[TotalCloudForcing15_30N_MWA, TotalCloudForcing15_30S_MWA, TotalCloudForcing15_30Dif_MWA] = PlotMWA(TotalCloudForcingMWA,14.5,30,'TotalCloudForcing');
[TotalCloudForcing30_49N_MWA, TotalCloudForcing30_49S_MWA, TotalCloudForcing30_49Dif_MWA] = PlotMWA(TotalCloudForcingMWA,30,48.5,'TotalCloudForcing');
[TotalCloudForcing49_90N_MWA, TotalCloudForcing49_90S_MWA, TotalCloudForcing49_90Dif_MWA] = PlotMWA(TotalCloudForcingMWA,48.5,90,'TotalCloudForcing');

%plot([TotalCloudForcing49_90S_MWA' TotalCloudForcing30_49S_MWA' TotalCloudForcing15_30S_MWA' TotalCloudForcing0_15S_MWA' TotalCloudForcing0_15N_MWA' TotalCloudForcing15_30N_MWA'  TotalCloudForcing30_49N_MWA' TotalCloudForcing49_90N_MWA'])
TotalCloudForcingMWA2Plot = [detrend(TotalCloudForcing49_90S_MWA',0) detrend(TotalCloudForcing30_49S_MWA',0) detrend(TotalCloudForcing15_30S_MWA',0) detrend(TotalCloudForcing0_15S_MWA',0) detrend(TotalCloudForcing0_15N_MWA',0) detrend(TotalCloudForcing15_30N_MWA',0)  detrend(TotalCloudForcing30_49N_MWA',0) detrend(TotalCloudForcing49_90N_MWA',0)];
for i=1:8
   hold on;
   plot(TotalCloudForcingMWA2Plot(:,i),'color',cmap(8*i,:),'LineWidth',3)
end
grid on;
set(gca,'xtick',[11:12:119+26])
set(gca,'XTickLabel',{2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013})
set(gca,'FontSize',20)
title('Mean-Subtracted TotalCloudForcing MWA Fluxes for Various Latitude Bands')
xlabel('Year end')
ylabel('Mean-Subtracted Flux (Watts/m^2)')
grid on;
set(legend(['90S-49S'],['49S-30S'],['30S-15S'],['15S-0S'],['0N-15N'],['15N-30N'],['30N-49N'],['49N-90N']))
set(gca,'GridLineStyle','--')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['TotalCloudForcing_MWA_Flux_All_Lats','.png']);
hold off;

NHMinusSH2(Net,0,14.5,'Net')
NHMinusSH2(Net,14.5,30,'Net')
NHMinusSH2(Net,30,48.5,'Net')
NHMinusSH2(Net,48.5,90,'Net')

NHMinusSH2(SW,0,14.5,'SW')
NHMinusSH2(SW,14.5,30,'SW')
NHMinusSH2(SW,30,48.5,'SW')
NHMinusSH2(SW,48.5,90,'SW')

NHMinusSH2(LW,0,14.5,'LW')
NHMinusSH2(LW,14.5,30,'LW')
NHMinusSH2(LW,30,48.5,'LW')
NHMinusSH2(LW,48.5,90,'LW')

NHMinusSH2(netclear,0,14.5,'netclear')
NHMinusSH2(netclear,14.5,30,'netclear')
NHMinusSH2(netclear,30,48.5,'netclear')
NHMinusSH2(netclear,48.5,90,'netclear')

NHMinusSH2(SWCF,0,14.5,'SWCF')
NHMinusSH2(SWCF,14.5,30,'SWCF')
NHMinusSH2(SWCF,30,48.5,'SWCF')
NHMinusSH2(SWCF,48.5,90,'SWCF')

NHMinusSH2(LWCF,0,14.5,'LWCF')
NHMinusSH2(LWCF,14.5,30,'LWCF')
NHMinusSH2(LWCF,30,48.5,'LWCF')
NHMinusSH2(LWCF,48.5,90,'LWCF')

NHMinusSH2(SWclear,0,14.5,'SWclear')
NHMinusSH2(SWclear,14.5,30,'SWclear')
NHMinusSH2(SWclear,30,48.5,'SWclear')
NHMinusSH2(SWclear,48.5,90,'SWclear')

NHMinusSH2(LWclear,0,14.5,'LWclear')
NHMinusSH2(LWclear,14.5,30,'LWclear')
NHMinusSH2(LWclear,30,48.5,'LWclear')
NHMinusSH2(LWclear,48.5,90,'LWclear')

NHMinusSH2(TotalCloudForcing,0,14.5,'TotalCloudForcing')
NHMinusSH2(TotalCloudForcing,14.5,30,'TotalCloudForcing')
NHMinusSH2(TotalCloudForcing,30,48.5,'TotalCloudForcing')
NHMinusSH2(TotalCloudForcing,48.5,90,'TotalCloudForcing')

NHMinusSH2(TotalCloudForcing,0,90,'TotalCloudForcing')

NetClimatologySubtracted = GlobalValuesMinusClimatology(Net,'Net');
SWClimatologySubtracted = GlobalValuesMinusClimatology(SW,'SW');
LWClimatologySubtracted = GlobalValuesMinusClimatology(LW,'LW');
SWCFClimatologySubtracted = GlobalValuesMinusClimatology(SWCF,'SWCF');
LWCFClimatologySubtracted = GlobalValuesMinusClimatology(LWCF,'LWCF');
netclearClimatologySubtracted = GlobalValuesMinusClimatology(netclear,'netclear');
LWclearClimatologySubtracted = GlobalValuesMinusClimatology(LWclear,'LWclear');
SWclearClimatologySubtracted = GlobalValuesMinusClimatology(SWclear,'SWclear');
TotalCloudForcingClimatologySubtracted = GlobalValuesMinusClimatology(TotalCloudForcing,'TotalCloudForcing');

HemisphericFluxesClimSubtracted(NetClimatologySubtracted,0,90,'Net') %do I really expect monthly-subtracted to be different...?
HemisphericFluxesClimSubtracted(NetClimatologySubtracted,0,30,'Net') %do I really expect monthly-subtracted to be different...?
HemisphericFluxesClimSubtracted(NetClimatologySubtracted,30,90,'Net') %do I really expect monthly-subtracted to be different...?

[NetNHMonthTotal, NetSHMonthTotal, NetDifMonthTotal] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,0,90,'Net');
[NetNHMonth1, NetSHMonth1, NetDifMonth1] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,0,14.5,'Net');
[NetNHMonth2, NetSHMonth2, NetDifMonth2] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,14.5,30,'Net');
[NetNHMonth3, NetSHMonth3, NetDifMonth3] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,30,48.5,'Net');
[NetNHMonth4, NetSHMonth4, NetDifMonth4] = HemisphericFluxesClimSubtracted(NetClimatologySubtracted,48.5,90,'Net');

HemisphericFluxesClimSubtracted(SWClimatologySubtracted,0,90,'SW')
HemisphericFluxesClimSubtracted(SWClimatologySubtracted,0,14.5,'SW')
HemisphericFluxesClimSubtracted(SWClimatologySubtracted,14.5,30,'SW')
HemisphericFluxesClimSubtracted(SWClimatologySubtracted,30,48.5,'SW')
HemisphericFluxesClimSubtracted(SWClimatologySubtracted,48.5,90,'SW')

HemisphericFluxesClimSubtracted(LWClimatologySubtracted,0,90,'LW')
HemisphericFluxesClimSubtracted(LWClimatologySubtracted,0,14.5,'LW')
HemisphericFluxesClimSubtracted(LWClimatologySubtracted,14.5,30,'LW')
HemisphericFluxesClimSubtracted(LWClimatologySubtracted,30,48.5,'LW')
HemisphericFluxesClimSubtracted(LWClimatologySubtracted,48.5,90,'LW')

HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,0,90,'netclear')
HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,0,14.5,'netclear')
HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,14.5,30,'netclear')
HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,30,48.5,'netclear')
HemisphericFluxesClimSubtracted(netclearClimatologySubtracted,48.5,90,'netclear')

HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,0,90,'SWCF')
HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,0,14.5,'SWCF')
HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,14.5,30,'SWCF')
HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,30,48.5,'SWCF')
HemisphericFluxesClimSubtracted(SWCFClimatologySubtracted,48.5,90,'SWCF')

HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,0,90,'LWCF')
HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,0,14.5,'LWCF')
HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,14.5,30,'LWCF')
HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,30,48.5,'LWCF')
HemisphericFluxesClimSubtracted(LWCFClimatologySubtracted,48.5,90,'LWCF')

HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,0,90,'SWclear')
HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,0,14.5,'SWclear')
HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,14.5,30,'SWclear')
HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,30,48.5,'SWclear')
HemisphericFluxesClimSubtracted(SWclearClimatologySubtracted,48.5,90,'SWclear')

HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,0,90,'LWclear')
HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,0,14.5,'LWclear')
HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,14.5,30,'LWclear')
HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,30,48.5,'LWclear')
HemisphericFluxesClimSubtracted(LWclearClimatologySubtracted,48.5,90,'LWclear')

HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,0,90,'TotalCloudForcing')
HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,0,14.5,'TotalCloudForcing')
HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,14.5,30,'TotalCloudForcing')
HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,30,48.5,'TotalCloudForcing')
HemisphericFluxesClimSubtracted(TotalCloudForcingClimatologySubtracted,48.5,90,'TotalCloudForcing')

%accuracy check - should be close to 0
%max(max(mean(LWClimatologySubtracted,3)))
%mean(max(mean(LWClimatologySubtracted,3)))

LWSWCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*SWClimatologySubtracted,3),(std(LWClimatologySubtracted,0,3).*std(LWClimatologySubtracted,0,3)));
LWSWCorrMap = bsxfun(@rdivide,155^2*sum(LWClimatologySubtracted.*SWClimatologySubtracted,3),(std(LWClimatologySubtracted,0,3).*std(LWClimatologySubtracted,0,3)));
LWSWCorrMap = bsxfun(@rdivide,mean(LWClimatologySubtracted.*SWClimatologySubtracted,3),(std(LWClimatologySubtracted,0,3).*std(LWClimatologySubtracted,0,3)));
%LWSWCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*SWClimatologySubtracted,3),(moment(LWClimatologySubtracted,2,3).*moment(LWClimatologySubtracted,2,3)));


LWSWCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*SWClimatologySubtracted,3),(sqrt(sum(LWClimatologySubtracted.^2,3)).*sqrt(sum(SWClimatologySubtracted.^2,3))));
LWNetCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*NetClimatologySubtracted,3),(sqrt(sum(LWClimatologySubtracted.^2,3)).*sqrt(sum(NetClimatologySubtracted.^2,3))));
SWNetCorrMap = bsxfun(@rdivide,sum(SWClimatologySubtracted.*NetClimatologySubtracted,3),(sqrt(sum(SWClimatologySubtracted.^2,3)).*sqrt(sum(NetClimatologySubtracted.^2,3))));
LW_LWClearCorrMap = bsxfun(@rdivide,sum(LWclearClimatologySubtracted.*LWClimatologySubtracted,3),(sqrt(sum(LWclearClimatologySubtracted.^2,3)).*sqrt(sum(LWClimatologySubtracted.^2,3))));
LW_LWCFCorrMap = bsxfun(@rdivide,sum(LWClimatologySubtracted.*LWCFClimatologySubtracted,3),(sqrt(sum(LWCFClimatologySubtracted.^2,3)).*sqrt(sum(LWClimatologySubtracted.^2,3))));
LWclear_SWclearCorrMap = bsxfun(@rdivide,sum(LWclearClimatologySubtracted.*SWclearClimatologySubtracted,3),(sqrt(sum(LWclearClimatologySubtracted.^2,3)).*sqrt(sum(SWclearClimatologySubtracted.^2,3))));

%9*8/2 = 36 possible pairwise comparisons

ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);land = shaperead('landareas', 'UseGeoCoords', true);geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(LWSWCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(LWNetCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(SWNetCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(LW_LWClearCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(LW_LWCFCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar
geoshow(LWclear_SWclearCorrMap,geoidrefvec,'DisplayType','texturemap');colorbar

GlobalCorrMap(LW,SW)
GlobalCorrMap(LWClimatologySubtracted,SWClimatologySubtracted)
GlobalCorrMap(LWclearClimatologySubtracted,SWClimatologySubtracted)
GlobalCorrMap(SWclearClimatologySubtracted,LWClimatologySubtracted);
GlobalCorrMap(SWCFClimatologySubtracted,LWCFClimatologySubtracted);
GlobalCorrMap(SWCFClimatologySubtracted,LWClimatologySubtracted);
GlobalCorrMap(NetClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(SWclearClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(LWclearClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(SWclearClimatologySubtracted,SWClimatologySubtracted);
GlobalCorrMap(LWClimatologySubtracted,LWCFClimatologySubtracted);
GlobalCorrMap(SWCFClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(LWCFClimatologySubtracted,TotalCloudForcingClimatologySubtracted);
GlobalCorrMap(NetClimatologySubtracted,netclearClimatologySubtracted);
GlobalCorrMap(TotalCloudForcingClimatologySubtracted,netclearClimatologySubtracted);



GlobalCorrMap(LWClimatologySubtracted,LW)
GlobalCorrMap(LW,LW);


%could also try for i=0:.25:0.75, arccos(i)to arccos(i+.25)
