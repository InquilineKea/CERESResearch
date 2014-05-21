close all; clear all; clc; warning off % get rid of java depreciation warning
%% Usar defined properties
% Dummy block for search path/directory declaration
fdir='http://www.atmos.uw.edu/~akchen0/CERES_Project/';
flist={'Net12MonthMA_HemisphericDifs_0-90.fig', 'Precip12MonthMA_HemisphericDifs_0-90.fig','NetClear12MonthMA_HemisphericDifs_0-90.fig'};
FigNumber=11; % Initialize figure numbering variable, should be >=2
plotsPerWindow=3; % Define subplots per window. Assumption, subplots are columns
subplotNumber=1; % Initialize subplotNumber. No need to change this.

% Location of legends wrt gcf. Determined manually. Will be a pain in butt for >500 images, can be automated, but I am bored now.
% To automate this, get OuterPos and Pos of each legend handle in set of
% figures. Uncomment line 68 and run code on all figs you need plotted.
% Then analyze the numbers and their trends. You should be able to predict
% those, w/ formula or interpolation (use splines for exact fit).
legendLocation={[17.4 7.4 3 2], ...
                [17.4 4.5 3 2],...
                [17.4 1.6 3 2]};

% Define the properties attached to axes objects that need extracting
axisParameters={'GridLineStyle','LineWidth','Title','XLabel','YLabel','YAxisLocation','XLim','XTick','XTickLabel','YLim','YTick','YTickLabel'};
%     get(gca,'GridLineStyle')           %char    % Equivalent manual implementation
%     get(gca,'LineWidth')               %double
%     get(get(gca,'Title'),'String')     %char
%     get(get(gca,'XLabel'),'String')    %char
%     get(get(gca,'YLabel'),'String')    %char
%     get(gca,'YAxisLocation')           %char <--- Marker for dataset on the right axis
%     get(gca,'XLim')                    %double array
%     get(gca,'XTick')                   %double array
%     get(gca,'XTickLabel')              %char array
%     get(gca,'YLim')                    %double array
%     get(gca,'YTick')                   %double array
%     get(gca,'YTickLabel')              %char array

% Define properties attached to lineseries objects that need extracting
lineParameters={'XData','YData','Color','DisplayName','LineStyle'};
%     get(b(1),'XData')                   %double array
%     get(b(1),'YData')                   %double array
%     get(b(1),'Color')                   %double array
%     get(b(1),'DisplayName')             %char
%     get(b(1),'LineStyle')               %char

%% Dummy block for file i/o and management
% Fetch files;
for i=1:length(flist)
    fn{i}=['F' num2str(i,'%.4d')]; % New file name, also gives us a var
    
    % Check if file exists. If it does, skip; if it dont, fetch. urlwrite
    % reads from the url and writes it to pwd
    if exist([fn{i},'.fig'],'file')==0, urlwrite([fdir, flist{i}],[fn{i} '.fig']); end
end

%% Main loop
for j=1:length(flist)
    % Open/Raise figure
    h1=openfig(fn{j},'new','invisible');
    
    % Get child objects attached to figure
    a=get(gcf,'Children');
    
    % Identify which one is legend object, the other two are assumed to be
    % axes objects from plotyy. Assumption that there's only 1 legend
    % object/orig. plot. There is no restriction on number of axes objects.
    for i=1:length(a),if isequal(get(a(i),'Tag'),'legend'), iLegend=i;end;end
    iAxis=setdiff(1:length(a),iLegend); % Separate legend from non-legend entries.
    
    % %     % Get data from legend. We only need names.
    % %     strLegend=get(a(iLegend),'String');
    % get(a(iLegend),'OuterPosition'), get(a(iLegend),'Position') 
    
    %% Get the object properties
    for u=1:length(iAxis) % Cycle through axes objects
        % Extract properties from current axes and store in cell container called axisProp
        for i=1:length(axisParameters) % Cycle through axes properties that need extraction
            pstr=axisParameters{i};
            % Get axis parameters
            if isequal('Title',pstr) || isequal('YLabel',pstr) || isequal('XLabel',pstr) % Check if the property is a child to axes
                axisProp(u,i)={get(get(a(iAxis(u)),pstr),'String')};
            else
                axisProp(u,i)={get(a(iAxis(u)),pstr)};
            end
        end
        
        b(u)={get(a(iAxis(u)),'Children')}; % Get the lineseries children attached to uth axis
        
        % Extract properties from lineseries objects and store in cell container called lineProp
        for i=1:length(b{u}) % Cycle through lineseries objects
            for k=1:length(lineParameters) % Cycle through lineseries properties that need extraction
                lineProp(u,i,k)={get(b{u}(i),lineParameters{k})};
            end
        end
    end
    
    %% Variables to control plot process
    % Create conditionals that prevent axes objects from being declared multiple times
    junkCreateRightAxis=true; junkCreateLeftAxis=true; junkExistAxes=false;
    
    %% Start the plot process
    f1=figure(FigNumber); vv=subplot(plotsPerWindow,1,subplotNumber); grid on;
    axes1=gca; % Define current subplot position to be the one to track
    for u=1:length(iAxis) % Cycle through axes objects
        for i=1:length(b{u}) % Cycle through lineseries attached to each axes object
            if isequal(axisProp{u,6},'right') % Distinguish between RHS or LHS axes
                if junkCreateRightAxis
                    set(axes1,'Color','none','YAxisLocation','right', 'XLim',axisProp{u,7}, 'YColor',lineProp{u,i,3}, 'XTick',axisProp{u,8}, 'XTickLabel',axisProp{u,9}, 'YLim',axisProp{u,10}, 'YTick',axisProp{u,11}, 'YTickLabel',axisProp{u,12});
                    junkCreateRightAxis=false; junkExistAxes=true;
                end
                p1=line(lineProp{u,i,1},lineProp{u,i,2},'Color',lineProp{u,i,3},'LineStyle',lineProp{u,i,5},'LineWidth',3,'Parent',axes1);
                %                 set()
                title(axisProp{u,3}),xlabel(axisProp{u,4}),ylabel(axisProp{u,5}),
            else % Plot on left axis
                if junkCreateLeftAxis
                    axes2=axes('Position',get(axes1,'Position'),'Color','none','YAxisLocation','left', 'XLim',axisProp{u,7},'YLim',axisProp{u,10},'XTick',axisProp{u,8}, 'XTickLabel',axisProp{u,9}, 'YTick',axisProp{u,11}, 'YTickLabel',axisProp{u,12});
                    junkCreateLeftAxis=false;
                end;
                line(lineProp{u,i,1},lineProp{u,i,2},'Color',lineProp{u,i,3},'LineWidth',3,'LineStyle',lineProp{u,i,5},'Parent',axes2);
                title(axisProp{u,3}),xlabel(axisProp{u,4}),ylabel(axisProp{u,5}),
            end
        end
    end % End of plot process
    
    % Convert the lowest axes background to white from transparent
    set(axes1,'Color',[1 1 1]);
    
    % Remove properties from legend object, then copy to current plot
    set(a(iLegend), 'FontSize',10, 'OuterPosition',legendLocation{j}, 'Position',legendLocation{j});
    copyobj(a(iLegend),f1);
    close(h1); % Close the original figure
    
    %% Update counters for figure & subplots
    subplotNumber=subplotNumber+1;
    if subplotNumber>plotsPerWindow
        subplotNumber=1;        % Step subplot number
        FigNumber=FigNumber+1;  % Step figure
        set(get(f1,'JavaFrame'),'Maximized',true); % Resize the plot window to fit the screen. Does not require localization.
    end
end
% Use cyclefigs to scroll through the mess you just created.
%eof-ssh-May 20, 2014
