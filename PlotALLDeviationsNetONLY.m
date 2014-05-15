% f = figure('Visible','off')
% subplot(2,2,1)
% plot(1:10)
% subplot(2,2,2)
% plot(5:15)
% 
% print(gcf,'-dpng','-r300',['TEST.png']);


% Net12MnthVid = VideoWriter('Net12MonthMA3.avi');
% Net12MnthVid.FrameRate = 5;
% open(Net12MnthVid);


load('MovingAvg.mat', 'MovingAvg');
FluxNames = fieldnames(MovingAvg.Window12Months);

%2000-03 to 2001-02


Month1=3;Year1=2000;
Month2=2;Year2=2001;
for k=1:size(MovingAvg.Window12Months.Net,3)
f = figure('Visible','off')
% f = figure('Visible','on')
    MakeLizMap;
colormap(lizmap)
set(gca,'FontSize',20)
    load geoid;ax = worldmap('World');setm(ax, 'Origin', [0 180 0]);
    curPlot = MovingAvg.Window12Months.Net(:,:,k);
    geoshow(curPlot,geoidrefvec,'DisplayType','texturemap');colorbar
    load coast;plotm(lat,long,'black')
    MaxBeta = max(max(abs(curPlot)));
%     caxis([-MaxBeta MaxBeta])
caxis([-20 20])
    title(['12-Month Net Anomaly on ' num2str(Year1), '-', num2str(Month1), ' to ', num2str(Year2),'-',num2str(Month2)]);

grid on;
% tightfig
 set(gcf,'Renderer','painters')
set(gcf,'paperposition',[0 0 20 10])
print(gcf,'-dpng','-r300',['NetONLY',num2str(k),'.png']);

% 
%     currFrame = getframe;
%     writeVideo(Net12MnthVid,currFrame);
    hold off;
    Month2 = Month2+1;Month1 = Month1+1;
    if eq(Month2,13)
        Month2 = 1;
    Year2 = Year2+1;
    end
    if eq(Month1,13)
        Month1=1;
        Year1=Year1+1;
    end
end

% close(Net12MnthVid);