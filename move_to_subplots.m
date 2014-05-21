function [hFigure]= move_to_subplots(ax,a,b)
% This is a function to automatically move a selection of axis in
% preexisting plots to one subplot figure. This is helpful for comparing
% large amounts of data.
%
% Inputs:
%		inputname: ax - array of axes handles.
%                   a (optional) - number of rows of figures in subplot.
%                   b (optional) - number of columns of figures in subplot.
% Outputs:
%		name:  figurewindow handle
%		plots: new figure window
%
% Standard call:
% for i=1:n_i
% figure(i)
% ax_h(i)=gca;
% end
% move_to_subplots(ax_h,2,ceil(n_i/2))
%
% Can also be used to copy a figure by  move_to_subplots(gca,1,1)
%
% Written by C. Hogg Date 2012_06_01
% Modified 20130808 Added the handle return
% 
debugmode=0;

hFigure=figure();

if ~exist('a')
        a=ceil(sqrt(length(ax)));
end

if ~exist('b')
        b=1;
    end    

if a*b<length(ax)|~exist('a')|~exist('b')
    disp('Auto subplot sizing')
    
    b=ceil(length(ax)/a);
end

for i=1:length(ax)

hTemp = subplot(a,b,i,'Parent',hFigure);         %# Make a new temporary subplot
newPos = get(hTemp,'Position');                  %# Get its position
delete(hTemp);

thisaxes=findobj(ax(i),'Type','axes');          % get the axes handle, even if ax(i) is a figure handle
hNew = copyobj(thisaxes,hFigure);
set(hNew,'Position',newPos)
end

%% Debug. Break point here.
if debugmode==1; dbstop tic; tic; dbclear all;end

end
