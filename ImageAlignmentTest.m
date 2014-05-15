[x,y,z] = peaks(50);
figure;
ah1 = subplot(10,1,1:9); %# capture handle of first axes
pcolor(x,y,z);
shading flat;
colorbar;
ah2 = subplot(10,1,10); %# capture handle of second axes
plot(x(end/2,:), z(end/2,:));

 subplots = get(gcf,'Children');
 AllPositions= get(subplots,'Position');


%# find current position [x,y,width,height]
pos2 = get(ah2,'Position');
pos1 = get(ah1,'Position');

%# set width of second axes equal to first
pos2(3) = pos1(3);
set(ah2,'Position',pos2)

set(ah1,'XTickLabel','')
pos2(2) = pos1(2) - pos2(4);
set(ah2,'Position',pos2)
