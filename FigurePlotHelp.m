


close all; clear all; clc
% Read png


[f1]=imread('f1.png'); [f2]=imread('f2.png');
% Plot the two images in 2x1 array, without using subimage or subplot
% function. We brute force it using the axes function. 
% % Note gcf background color is set to yellow, and axis are turned off
[f1]=imread('aaa.png');[f2]=imread('aaa.png');
figure(1), hold on; set(gcf,'Color',[1 1 0]), set(gca,'Visible','off'), set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1])
axes('position',[.25 0 .5 .5]), image(f1), set(gca,'Visible','off')
axes('position',[.25 .5 .5 .5]), image(f2), set(gca,'Visible','off'), axis manual
% Create a 2x2 array of images. Also using axes. Bckgrnd not visible here.
figure(2), hold on; set(gca,'Visible','off'), set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1])
axes('position',[0 0 .5 .5]), image(f1), set(gca,'Visible','off')
axes('position',[.5 0 .5 .5]), image(f2), set(gca,'Visible','off')
axes('position',[.5 .5 .5 .5]), image(f1), set(gca,'Visible','off')
axes('position',[0 .5 .5 .5]), image(f2), set(gca,'Visible','off'), axis manual
% % If you want to use scaling/shifting
close all; clear all; clc
% Read png
[f1]=imread('f1.png'); [f2]=imread('f2.png');
R=0.6; m=1.3;
% Plot the two images in 2x1 array, without using subimage or subplot
% function. We brute force it using the axes function. 
% % Note gcf background color is set to yellow, and axis are turned off
figure(1), hold on; set(gcf,'Color',[1 1 0]), set(gca,'Visible','off'), set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1])
axes('position',[.1 0 m*R m*(1-R)]), image(f1), set(gca,'Visible','off')
axes('position',[.1 .5 m*R m*(1-R)]), image(f2), set(gca,'Visible','off'), axis manual
m=1;
% Create a 2x2 array of images. Also using axes. Bckgrnd not visible here.
figure(2), hold on; set(gca,'Visible','off'), set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1])
axes('position',[0   0 m*R m*(1-R)]), image(f1), set(gca,'Visible','off')
axes('position',[.5  0 m*R m*(1-R)]), image(f2), set(gca,'Visible','off')
axes('position',[.5 .5 m*R m*(1-R)]), image(f1), set(gca,'Visible','off')
axes('position',[0  .5 m*R m*(1-R)]), image(f2), set(gca,'Visible','off'), axis manual
% Update the axes array here to increase the offset between the images and reduce m



close all;
% Create png for demo
t=0:1e-3:5; plot(t,sin(2*pi*3*t),'r','LineWidth',2), saveas(gca,'f1.png'), plot(t,cos(2*pi*3*t),'b','LineWidth',2), saveas(gca,'f2.png')
close all; clear all; clc
% Read png
LowerLat = 0;
HigherLat = 90;
MonthFilterSize = 12;
load('2014-04-30.mat', 'FluxNames')

[f1]=imread([FluxNames{11}, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.png']); 
[f2]=imread([FluxNames{1}, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.png']);
% Plot the two images in 2x1 array, without using subimage or subplot
% function. We brute force it using the axes function. 
% % Note gcf background color is set to yellow, and axis are turned off
figure(1), hold on; set(gcf,'Color',[1 1 0]), set(gca,'Visible','off')
axes('position',[.25 0 .5 .5]), image(f1), set(gca,'Visible','off')
axes('position',[.25 .5 .5 .5]), image(f2), set(gca,'Visible','off'), axis manual

print(gcf,'-dpng','-r300',['2By1.png']);

% Create a 2x2 array of images. Also using axes. Bckgrnd not visible here.
figure(2), hold on; set(gca,'Visible','off')
axes('position',[0 0 .5 .5]), image(f1), set(gca,'Visible','off')
axes('position',[.5 0 .5 .5]), image(f2), set(gca,'Visible','off')
axes('position',[.5 .5 .5 .5]), image(f1), set(gca,'Visible','off')
axes('position',[0 .5 .5 .5]), image(f2), set(gca,'Visible','off'), axis manual

print(gcf,'-dpng','-r300',['blah.png']);


LowerLat = 0;
HigherLat = 90;
MonthFilterSize = 12;
load('2014-04-30.mat', 'FluxNames')
h1 = openfig([FluxNames{11}, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.fig'],'reuse');
ax1 = gca;
gcf1 = gcf;
Lines1 = findobj('type','line');
h2 = openfig([FluxNames{1}, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.fig'],'reuse')
ax2 = gca;
gcf2 = gcf;
% Lines2 = findobj('type','line'); %%not sure about this
Lines2 = findobj(ax2,'type','line'); %%not sure about this

h3 = figure; %// create new figure
s1 = subplot(2,1,1); % // create and get handle to the subplot axes
s2 = subplot(2,1,2);
% /// Step #4
fig1 = get(ax1,'children'); %// get handle to all the children in the figure
fig2 = get(ax2,'children');
 copyobj(Lines1,s1)
set(s1,'XTick',get(ax1,'XTick'),'XTickLabel',get(ax1,'XTickLabel'),'YTick',get(ax1,'YTick'),'YTickLabel',...
    get(ax1,'YTickLabel'),'YLim',get(ax1,'YLim'))
ax = findobj(gcf1,'Type','axes','Tag','');
y1Label = get(get(ax(1),'YLabel'),'String')
y2Label = get(get(ax(2),'YLabel'),'String')
xLabel  = get(get(ax(2),'XLabel'),'String')
%%%%%
h = get(gcf1,'children');
hLeg = [];
for k = 1:length(h)
    if strcmpi(get(h(k),'Tag'),'legend')
        hLeg = h(k);
        break;
    end
end
legend1=legend(s1,'show')
set(legend1,get(hLeg,'String')) %this is the legend!
 copyobj(Lines1,s2)
legend2=legend(s2,'show')
set(legend2,get(hLeg,'String')) %this is the legend!
set(legend2,{'hihi'}) %this is the legend!

%%%%%%%%%%%%

junk=figure, plot(nan,nan);

openfig('Net12MonthMA_HemisphericDifs_0-90.fig','new','invisible');
gcfNew = gcf;
a=get(gcf,'Children');
leg = findobj(gcfNew,'Tag','legend')
for i=1:length(a),
    if isequal(get(a(i),'Tag'),'legend'),
        copyobj(a(i),junk)
%         copyobj(a(i),get(s1,'parent'))
    end;
end

copyobj(a(1),get(s1,'parent')) %THIS WORKS but only for one of the plots
copyobj(leg,get(s1,'parent')) %THIS WORKS but only for one of the plots

set(legend1,a(i))
set(get(legend1,'parent'),a(i))

% make a figure with several axes

fig1 = figure;

xx = 0:pi/10:2*pi;

sp(1) = subplot(3,1,1);

plot(xx, 10*sin(xx));

sp(2) = subplot(3,1,2);

plot(xx, cos(xx));

sp(3) = subplot(3,1,3);

plot(xx, tan(xx));

% Create a Legend for the first axes

hLeg = legend(sp(1),'Signal')

% create a new figure for saving and printing

fig2 = figure('visible','off');

% copy axes into the new figure

newax = copyobj(sp(1),fig2); 

% Since you have a LEGEND associated with the figure in the first axes 

% you can also copy the legend to the new figure:

newLeg = copyobj(hLeg,fig2); 

% If you would like to have the axes cover a larger area of the figure window rather

% than the original size as in the subplot, then change it's Position value:

set(newax, 'units', 'normalized', 'position', [0.13 0.11 0.775 0.815]);

% print and/or save the figure

print(fig2)             % print it

hgsave(fig2,'myfig')    % save it

close(fig2)             % clean up by closing it



hFig = figure;

%# create temporary subplots as template
for i=1:2, h(i) = subplot(2,1,i); end       %# create subplots
pos = get(h, 'Position');                   %# record their positions
delete(h)                                   %# delete them

%# load the .fig files inside the new figure
fileNames = {'a.fig' 'b.fig'};              %# saved *.fig file names
for i=1:2
    %# load fig
    hFigFile = hgload( fileNames{i} );

    %# move/copy axis from old fig to new fig
    hAx = get(hFigFile, 'Child');           %# hAx = gca;
    set(hAx, 'Parent',hFig)
    %#hAx = copyobj(hAx,hFig);

    %# resize it to match subplot position
    set(hAx, 'Position',pos{i});

    %# delete old fig
    delete(hFigFile)
end



 figure(1);
clf(1);
axsource=axes('Parent',1);
plot(axsource,rand(10,1));
hgsave(1,'sourcefig.fig');
close all;

% Load the fig file
h=hgload([FluxNames{11}, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.fig'])
set(h,'Visible','off');
tmpaxes=findobj(h,'Type','axes');

% Copy the axe
figure(2)
destaxes=subplot(2,2,1,'Parent',2)
copyobj(get(tmpaxes,'children'),destaxes);

copyobj(allchild(tmpaxes),destaxes);

% Clean up
close(h)




























% create a fig
     fnam='foo.fig';
     fh=figure;
     surfl(peaks(32));
     shading interp;
     camlight right;
     lighting phong;
     title('TEST');
     saveas(gcf,fnam);
     delete(fh);
% the engine
% ...a subplot template
     fhs=figure;
for i=1:2
     sh(i)=subplot(2,2,i);
end
% ...reload the figs (here we use the same fig!)
for i=1:2
     fh=open([FluxNames{11}, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.fig']);
     ah=gca;
% ...move it to the subplot
     ch=copyobj(ah,fhs);
% ...resize it
     set(ch,'position',get(sh(i),'position'));
% ...and delete the fig's canvas
     delete(fh);
end
% ...delete the template
     delete(sh);
figs2subplots([a1 a2 a3],[2 2],{[1 3],2,4})
     
     
ohf = figure
ohp= uipanel( ohf,'units','normal','pos',[0,0,1/2,1/2] )
ohs=hgload([FluxNames{11}, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.fig']);
set( get(ohs,'child'), 'parent', ohp )
close(ohs)
       fh=open([FluxNames{11}, num2str(MonthFilterSize),'MonthMA_HemisphericDifs_',num2str(LowerLat),'-',num2str(HigherLat),'.fig']);

unsubplot(2,1,1,fh)



hFig = figure;

%# create temporary subplots as template
for i=1:2, h(i) = subplot(2,1,i); end       %# create subplots
pos = get(h, 'Position');                   %# record their positions
delete(h)                                   %# delete them

%# load the .fig files inside the new figure
fileNames = {'a.fig' 'b.fig'};              %# saved *.fig file names
for i=1:2
    %# load fig
    hFigFile = hgload( fileNames{i} );

    %# move/copy axis from old fig to new fig
    hAx = get(hFigFile, 'Child');           %# hAx = gca;
    set(hAx, 'Parent',hFig)
    %#hAx = copyobj(hAx,hFig);

    %# resize it to match subplot position
    set(hAx, 'Position',pos{i});

    %# delete old fig
    delete(hFigFile)
end



%%%%%

copyobj(findobj(gcf1,'Type','axes','Tag','legend'),s1)
copyobj(findobj(ax1,'Type','axes','Tag','legend'),s1)

copyobj(findobj('type','axes'),s1)

 axes_handle = findobj(h1, 'Type', 'Axes'); % find handle to axes in figure
axes_children_handle = get(axes_handle, 'Children');
copyobj(axes_children_handle, s1); % children of original axes are copied to new axes

s1Handle = get(s1,'children');
get(s1Handle,'type')

 figure_children = get(gcf1,'Children');
children_axes = findall(figure_children,'Type','axes');
copyobj(children_axes,s1)

% Identify axes to be copied 
axes_to_be_copied = findobj(h1,'type','axes'); 
% Identify the children of this axes 
chilred_to_be_copied = get(axes_to_be_copied,'children'); 
% Identify orientation of the axes 
[az,el] = view; 
% Copy the children of the axes 
copyobj(chilred_to_be_copied,s1); 

% /// Step #5
copyobj(gcf1,s1);

copyobj(fig1,s1); %// copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);

