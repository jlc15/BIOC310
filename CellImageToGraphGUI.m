function varargout = CellImageToGraphGUI(varargin)
% CELLIMAGETOGRAPHGUI MATLAB code for CellImageToGraphGUI.fig
%      CELLIMAGETOGRAPHGUI, by itself, creates a new CELLIMAGETOGRAPHGUI or raises the existing
%      singleton*.
%
%      H = CELLIMAGETOGRAPHGUI returns the handle to a new CELLIMAGETOGRAPHGUI or the handle to
%      the existing singleton*.
%
%      CELLIMAGETOGRAPHGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLIMAGETOGRAPHGUI.M with the given input arguments.
%
%      CELLIMAGETOGRAPHGUI('Property','Value',...) creates a new CELLIMAGETOGRAPHGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellImageToGraphGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellImageToGraphGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellImageToGraphGUI

% Last Modified by GUIDE v2.5 21-Feb-2017 18:01:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CellImageToGraphGUI_OpeningFcn, ...
    'gui_OutputFcn',  @CellImageToGraphGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CellImageToGraphGUI is made visible.
function CellImageToGraphGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellImageToGraphGUI (see VARARGIN)

% Choose default command line output for CellImageToGraphGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CellImageToGraphGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CellImageToGraphGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function colonyObj = standardizeMATFILE(colonies)
%Input - empty first colony struct of type dynCell
%Output - removed empty cells populate 'Colony' object

%             if isempty(handles.colonies(colony).cells.onframes)
%                 continue
%             end
colonyObj = Colony();
for colony = 1 : length(colonies)
    numCells = length(colonies(colony).cells);
    for eachCell = 1 : numCells
        if isempty(colonies(colony).cells(eachCell).position);
            continue
        else
            %build deep copy of original colony structure without empty
            %cells
            colonyObj = Colony.addCell(colonyObj, colonies(colony).cells(eachCell));
        end
    end
end

% Executes on button press in loadMATFILEButton and constructs Pop-Up window and
% populates handles structure with .mat analysis information. 
% Plots .mat data on GRAPHAXES
function loadMATFILEButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadMATFILEButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Output - handles structure with updated information from .mat file and
% rendered graph on GUI

%TODO maybe syntax to make empty string
inputWindow.Mat_File = { {'uigetfile(''.'')'} };
inputWindow = StructDlg(inputWindow, 'Data To Plot');
handles.matFileName = inputWindow.Mat_File;
load(handles.matFileName);

%new data contains colony information
if exist('colonies','var')
    %maybe have different init graph calls here depends on if
    %colonies/cells provided
    %nonEmptyCells = standardizeMATFILE(colonies);
    handles.colonies = colonies;
elseif exist('cells', 'var')
    %     handles.cells = cells;
end

%possibly use isLoaded to tell user that the .mat file and the images
%don't match up
handles.loadedMAT = 1;
handles.colonyColors = [];
%Plots in GUI the MATFILE contents
colonyColors = initGraph(hObject, handles);
handles.colonyColors = colonyColors;
set(handles.GraphAxes, 'HandleVisibility', 'off'); % AN
guidata(hObject,handles);


function plotSpecificColony(handles, colony, colonyColors)
axes(handles.GraphAxes);
currentColony = handles.colonies(colony);
totalCells = length(currentColony.cells);
%only start plotting at idx <- 2
for eachCell = 2 : totalCells                                       %%% AN
    currentCell = currentColony.cells(eachCell);
    cytoToNuclearFluor = currentCell.fluorData(:, 2)./ ...
        currentCell.fluorData(:, 3);
    plot(handles.GraphAxes, currentCell.onframes, cytoToNuclearFluor,...
        '-o', 'Color', colonyColors(colony, :), 'MarkerSize', 5,...
        'MarkerEdgeColor', colonyColors(colony,:),...
        'MarkerFaceColor', colonyColors(colony, :) );
    hold on; 
    %store cell traces
    %cellTrace = [(ones(1, length(currentCell.onframes)).*colony)' (ones(1, length(currentCell.onframes)).*eachCell)' currentCell.onframes cytoToNuclearFluor]; %% AN currentCell.onframes'
    %originalCellTraces = [originalCellTraces; cellTrace];
    
end

%Plots the given colonies
function colonyColors = helpPlot(hObject, handles, colN)
%xlim MAX
%globalMaxT = -Inf;
%ylim MAX
%globalMaxCytoToNuclearFluor = -Inf;
    %     currentMaxT = max(currentCell.onframes);
    %update global xlim and ylim                                    %% AN, don;t update the x and
    %y lims, keep them constant at max time point overall and max
    %possible value of signaling
    %     if currentMaxT > globalMaxT
    %         globalMaxT = currentMaxT;
    %     end
    
    %     if max(cytoToNuclearFluor) > globalMaxCytoToNuclearFluor
    %         globalMaxCytoToNuclearFluor = max(cytoToNuclearFluor);
    %     end
if isempty(handles.colonyColors)
    %init colonyColors once
    colonyColors = rand(3, length(handles.colonies))';
else
    %don't want to keep updating colonyColors once initialized
    colonyColors = handles.colonyColors;
end

if ~isempty(colN)                                               % AN plot only specific colony
    %call plot colN
    plotSpecificColony(handles, colN, colonyColors)
else
    for colony = 1 : length(handles.colonies)
        plotSpecificColony(handles, colony, colonyColors)
        hold on;
    end
end
%TO KEEP TRACK OF WHICH COLONY IS WHICH
%     text(currentMaxT - 20, 0 + globalMaxCytoToNuclearFluor / 10 , ...
%         ['Colony: ' num2str(colony)], 'Color', colonyColors{colony});
% text(1, 1.5 , ...                                                                  %%% AN
%  ['Colony: ' num2str(colony)], 'Color', colonyColors(colony,:));

xlim(handles.GraphAxes, [0 99]);                                                       % AN, possibly set when the person uploads the mat file, since it is known how many time points there is in advance
ylim(handles.GraphAxes, [0 2]);                                                        % AN
xlabel(handles.GraphAxes, 'Frames');
ylabel(handles.GraphAxes, 'Nuclear to Cytosolic Fluorescence Ratio');
set(gca, 'fontsize', 16);
box on;
guidata(hObject, handles);

%Initialize graph view to visualize all colonies
function colonyColors = initGraph(hObject, handles)
if isfield(handles,'colonies')
    %originally set graphvis to view all colonies
    handles.coloniesToPlot = ones(1 : length(handles.colonies));
    colonyColors = helpPlot(hObject, handles, []);            % AN initially want to look at all colonies
end

if (length(handles.colonies) > 1)
    %now that you have plotted all colonies you can select specific colonies if
    %more than one colony exists
    set(handles.selectColoniesButton, 'Visible', 'on');
end
guidata(hObject,handles);
%set(gcf,'WindowButtonMotionFcn', {@getPixelOnMouseMove, handles.GraphAxes});


%Generates the max intensity projection for each channel once and stores
function maxIntensityProjections = getMaxIntensityProjections(images, hObject, handles)
%parallelize?? <- parfor()
%cellstr(images.ordering')'
%---------AN          commented out to troubleshoot the image displaying
%part without loading the .mat file
% totalTime = 0;
% for colony = 1 : length(handles.colonies)
%     cells = cat(2, handles.colonies(colony).cells.onframes);
%     maxColonyTime = max(cells');
%     if totalTime < maxColonyTime
%         totalTime = maxColonyTime;
%     end
% end
%-------------
totalTime = 2;% for the purposes of troubleshooting look at the first two time points
%Init max intensity projection storage for all timepoints for all channels
maxIntensityProjections = cell(1, totalTime);
numberOfChannels = length(images.w);
for timept = 1 : totalTime % parfor
    %assume all channels and timepoints are the same length
    maxIntensityProjections{timept} = ...
        zeros(1024, 1024, numberOfChannels); %TODO: get dim without bfgetreader
    for channel = 1 : numberOfChannels
        %get max intensity z slice
        currentChannel = ...
            andorMaxIntensityBF(images, handles.microscopePosition, timept - 1, images.w(channel));
        maxIntensityProjections{timept}(:, :, channel) = currentChannel;
    end
end


%Initializes Image Viewer Panel based on user input
function maxIntensityProjections = initImageViewer(hObject, eventdata, handles)
axes(handles.ImageAxes);
images = readAndorDirectory(handles.imageDirectory);
maxIntensityProjections = getMaxIntensityProjections(images, hObject, handles);
handles.maxIntensityProjections = maxIntensityProjections;
set(handles.scrollImageSet, 'Min', 1);                                       % handles is empty, so returs the error when try to set
set(handles.scrollImageSet, 'Max', length(handles.maxIntensityProjections));
set(handles.scrollImageSet, 'Value', 1); %init at timept 1
set(handles.scrollImageSet, 'SliderStep', [1 / (length(handles.maxIntensityProjections)) ,...
    10 / (length(handles.maxIntensityProjections)) ]);
Nuclear_Callback(hObject, eventdata, handles);
Cytosolic_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in loadTIFFILE. Occurs after loading .mat
% file
function loadTIFFILE_Callback(hObject, eventdata, handles)
% hObject    handle to loadTIFFILE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image.Directory = { {'uigetdir(''.'')'} };
image.MicroscopePosition = 0;
image = StructDlg(image);
handles.microscopePosition = image.MicroscopePosition;
handles.imgTimept = 1; %arbitrarily picked to start time from 1
handles.imageDirectory = image.Directory;
maxIntensityProjections = initImageViewer(hObject, eventdata, handles);
handles.maxIntensityProjections = maxIntensityProjections;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function scrollImageSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scrollImageSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on slider movement.
function scrollImageSet_Callback(hObject, eventdata, handles)
% hObject    handle to scrollImageSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%axes(handles.ImageAxes);
%hold on;
imageNumberSelected = int32(get(handles.scrollImageSet, 'Value'));
%get(handles.imadjust(img{channelIdx},...
%stretchlim(img{channelIdx}, sclimits));
%imshow(handles.maxIntensityProjections{imageNumberSelected}(:, :, 2), []);

%Update Image Display based on channel selection
handles.imgTimept = imageNumberSelected;
Nuclear_Callback(hObject, eventdata, handles)
Cytosolic_Callback(hObject, eventdata, handles)
cellsToEnlarge = getXYColonyTime(handles, handles.coloniesToPlot);
delete(allchild(handles.GraphAxes));
helpPlot(hObject, handles, handles.coloniesToPlot);

plot(handles.GraphAxes, cellsToEnlarge(:, 1), cellsToEnlarge(:, 4), ...
            'Marker', '.','Markersize', 50, 'Color', handles.colonyColors(handles.coloniesToPlot, :));
%CellLabelling_Callback(hObject, eventdata, handles)

guidata(hObject, handles);

% --- Executes on button press in Nuclear.
function Nuclear_Callback(hObject, eventdata, handles)
% hObject    handle to Nuclear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ImageAxes);
hold on;
nuc = handles.maxIntensityProjections{handles.imgTimept}(:, :, 1);
cyto = handles.maxIntensityProjections{handles.imgTimept}(:, :, 2);
displayNuc = get(handles.Nuclear,'Value');
set(handles.Nuclear, 'Value', displayNuc);
if displayNuc
    if get(handles.Cytosolic, 'Value')
        imshow(imfuse(nuc, cyto), []);
    else
        imshow(nuc, []);
    end
else
    %just display cyto and not nucl bc can't show nothing
    set(handles.Cytosolic, 'Value', 1);
    imshow(cyto, []);
end
if get(handles.CellLabelling, 'Value')
    %displays cell labels on ImageAxes
    CellLabellingTrue(hObject, eventdata, handles)
end
%guidata(hObject, handles);

% --- Executes on button press in Cytosolic.
function Cytosolic_Callback(hObject, eventdata, handles)
% hObject    handle to Cytosolic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ImageAxes);
nuc = handles.maxIntensityProjections{handles.imgTimept}(:, :, 1);
cyto = handles.maxIntensityProjections{handles.imgTimept}(:, :, 2);
displayCyt = get(handles.Cytosolic, 'Value');
set(handles.Cytosolic, 'Value', displayCyt);
if displayCyt
    if get(handles.Nuclear, 'Value')
        imshow(imfuse(nuc, cyto), []);
    else
        imshow(cyto, []);
    end
else
    %just display nuc and not cyto
    set(handles.Nuclear, 'Value', 1);
    imshow(nuc, []);
end
if get(handles.CellLabelling, 'Value')
    %displays cell labels on ImageAxes
    CellLabellingTrue(hObject, eventdata, handles)
end
guidata(handles.output, handles);

function CellLabellingTrue(hObject, eventdata, handles)
%Called if User wants Cell Labelling to Occur
axes(handles.ImageAxes);
if handles.loadedMAT  %check if there is a relevant "lookup table" a.k.a. the .mat file
    for colony = 1 : length(handles.colonies)
        cellInfo = getXYColonyTime(handles, colony);
        plot(cellInfo(:, 2), cellInfo(:, 3), 'LineStyle', 'none',...
            'Marker', '.','Markersize', 30, 'Color', handles.colonyColors(colony,:));
        %         if colony == 1
        %             %init array for plotting updater
        %             colors = ones(1, length(cellInfo)).*colony;
        %         else
        %             colors = [colors  ones(1, length(cellInfo)).*colony];
        %         end
    end
end
%updates the graph for the current timept by enlarging the dots
%EnlargeDots_Update(hObject, eventdata, handles, colors);                                   % AN
%don't need
guidata(hObject, handles);

% --- Executes on button press in CellLabelling. Updates ImageAxes and
% Calls Fxn to handle GraphAxes Update.
function CellLabelling_Callback(hObject, eventdata, handles)
% hObject    handle to CellLabelling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%TODO: redundant calculations here!!! cache calc values
displayCellLabel = get(handles.CellLabelling, 'Value');
set(handles.CellLabelling, 'Value', displayCellLabel);
if displayCellLabel %add cells labels if desired
    CellLabellingTrue(hObject, eventdata, handles)
else
    %you want the cell labelling display off
    %h = findobj(gca,'Type','line');
    %[x, y]= findall(handles.ImageAxes.Children,'Type','Image')
    delete(allchild(handles.ImageAxes));
    initImageViewer(hObject, eventdata, handles);
end
%update Graph as well
%getXYColonyTime(handles)

function cellInfo = getXYColonyTime(handles, colony)
%Called When cell labelling is selected in slider
%Output - returns x, y coordinates to plot for a given timept for a given
%colony
times = cat(1, handles.colonies(colony).cells.onframes);%% AN
coordinates = cat(1, handles.colonies(colony).cells.position);
allSignal = cat(1, handles.colonies(colony).cells.fluorData);
signalData = allSignal(:, 2) ./ allSignal(:, 3);
search = [times, coordinates, signalData];% AN
cellsToPlot = find(search(:, 1) == (handles.imgTimept));
cellInfo = search(cellsToPlot, :);




% --- Executes on button press in selectColoniesButton.
%Creates StructDlg for selection of specific colonies to plot
function selectColoniesButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectColoniesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(allchild(handles.GraphAxes));
col = size(handles.colonies, 2);                                                                 %  AN
userSelectedColonies.SelectedColonies = 0;% have the box empty before selection of colony to plot

%TODO: throw exception if user tries to not put in units/etc.
userSelectedColonies = StructDlg(userSelectedColonies, 'Select Colonies');
%handles.coloniesToPlot = handles.colonies(userSelectedColonies.SelectedColonies);                  %AN
handles.coloniesToPlot = userSelectedColonies.SelectedColonies;
colN = userSelectedColonies.SelectedColonies;
helpPlot(hObject, handles, colN);                                                                % AN which colony to plot, comes from user input
guidata(hObject,handles);%%%

function GraphAxes_CreateFcn(hObject, eventdata, handles)
disp('called before plotting anything');
function EnlargeDots_CreateFcn(hObject, eventdata, handles)
set(gca, 'xticklabel', [], 'yticklabel', []);
