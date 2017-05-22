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

% Last Modified by GUIDE v2.5 23-Mar-2017 12:41:38

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
handles.loadedMAT = 0;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = CellImageToGraphGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Executes on button press in loadMATFILEButton and constructs Pop-Up window and
% populates handles structure with .mat analysis information.
% Plots .mat data on GRAPHAXES
function loadMATFILEButton_Callback(hObject, eventdata, handles)
% Inputs -
% hObject    handle to loadMATFILEButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Output - handles structure with updated information from .mat file and
% rendered graph on GUI
inputWindow.Mat_File = { {'uigetfile(''.'')'} };
inputWindow.MicroscopePosition = 0;
inputWindow = StructDlg(inputWindow, 'Data To Plot');
handles.matFileName = inputWindow.Mat_File;
handles.MATMicroscopePosition = inputWindow.MicroscopePosition;
load(handles.matFileName);
if exist('colonies','var')
    handles.colonies = colonies;
    handles.withoutColonyInfo = 0;
elseif exist('positions', 'var')
    %there are no unique colony groupings; everything is in one "colony"
    handles.matFileNameForPosition = positions(handles.MATMicroscopePosition + 1).filename;
    handles.colonies = positions(handles.MATMicroscopePosition + 1);
    handles.withoutColonyInfo = 1;
end

%TODO: possibly use isLoaded to tell user that the .mat file and the images don't match up
handles.loadedMAT = 1;
handles.colonyColors = [];
handles.coloniesToPlot = [];
%Plots in GUI the MATFILE contents
colonyColors = initGraph(hObject, handles);
handles.colonyColors = colonyColors;
set(handles.GraphAxes, 'HandleVisibility', 'off');
guidata(hObject,handles);

function [globalXMax, globalYMax] = plotSpecificColonyHelperColony(handles, colony, colonyColors, globalXMax, globalYMax)
currentColony = handles.colonies(colony);
totalCells = length(currentColony.cells);
%only start plotting at idx <- 2 because some files have an empty first
%cell
for eachCell = 2 : totalCells
    currentCell = currentColony.cells(eachCell);
    cytoToNuclearFluor = currentCell.fluorData(:, 2) ./ ...
        currentCell.fluorData(:, 3);
    if max(currentCell.onframes) > globalXMax
        globalXMax = max(currentCell.onframes);
    end
    if max(cytoToNuclearFluor) > globalYMax
        globalYMax = max(cytoToNuclearFluor);
    end
    plot(handles.GraphAxes, currentCell.onframes, cytoToNuclearFluor,...
        '-o', 'Color', colonyColors(colony, :), 'MarkerSize', 5,...
        'MarkerEdgeColor', colonyColors(colony, :),...
        'MarkerFaceColor', colonyColors(colony, :) );
    text(2, (- colony / 15 + 0.25), ['Colony: ' num2str(colony)], ...
        'Color', colonyColors(colony, :), 'Parent', handles.GraphAxes);
end

function [globalXMax, globalYMax] = plotSpecificColonyHelperNoColony(handles, colony, colonyColors, globalXMax, globalYMax)
dasCells = handles.colonies.cellData;
toPlotMean = (cat(1, dasCells.nucLevelAvg) - cat(1, dasCells.background))...
    ./ (cat(1, dasCells.cytLevelAvg) - cat(1, dasCells.background));
err = zeros(1, handles.colonies.tPerFile);
sampleSize = zeros(1, handles.colonies.tPerFile);
for timept = 1 : handles.colonies.tPerFile
    someValNAN = dasCells(timept).nucLevel ./ dasCells(timept).cytLevel;
    nanidx = isnan(someValNAN);
    toTakeSTD = someValNAN(nanidx ~= 1);
    if max(toTakeSTD) > globalYMax
        globalYMax = max(toTakeSTD);
    end
    %display sample size on Image Axes
    sampleSize(timept) = length(toTakeSTD);
    err(timept) = std(toTakeSTD);
end
errorbar(handles.GraphAxes, 1 : handles.colonies.tPerFile, toPlotMean,...
    err, 'Color', colonyColors(colony, :), 'MarkerSize', 5,...
    'MarkerEdgeColor', colonyColors(colony, :),...
    'MarkerFaceColor', colonyColors(colony, :), 'LineWidth', 1.5);
hold on;
globalXMax = handles.colonies.tPerFile;

function [globalXMax, globalYMax] = plotSpecificColony(handles, colony, colonyColors, globalXMax, globalYMax)
%TODO: Refactor into two separate helper functions
axes(handles.GraphAxes);
hold on;
if handles.withoutColonyInfo
    %No Specific Colonies therefore plot averages
    [globalXMax, globalYMax] = plotSpecificColonyHelperNoColony(handles, colony, colonyColors, globalXMax, globalYMax);
else
    %Specific colony data
    [globalXMax, globalYMax] = plotSpecificColonyHelperColony(handles, colony, colonyColors, globalXMax, globalYMax);
end

%Plots the given colonies
function colonyColors = helpPlot(hObject, handles, colN)
if isempty(handles.colonyColors)
    %init colonyColors only once
    colonyColors = rand(3, length(handles.colonies))';
else
    %do not want to keep changing colonyColors once initialized
    colonyColors = handles.colonyColors;
end

globalXMax = -Inf;
globalYMax = -Inf;
if ~isempty(colN)
    %plot colN, the specific colony desired by the user
    [globalXMax, globalYMax] = plotSpecificColony(handles, colN, colonyColors, globalXMax, globalYMax);
else
    for colony = 1 : length(handles.colonies)
        [currentXMax, currentYMax] = plotSpecificColony(handles, colony, colonyColors, globalXMax, globalYMax);
        if currentXMax > globalXMax
            globalXMax = currentXMax;
        end
        if currentYMax > globalYMax
            globalYMax = currentYMax;
        end
        if ~handles.withoutColonyInfo
            text(2, (- colony / 15 + 0.25), ['Colony: ' num2str(colony)], ...
                'Color', colonyColors(colony, :), 'Parent', handles.GraphAxes);
            hold on;
        end
    end
end

xlim(handles.GraphAxes, [0 globalXMax]);
ylim(handles.GraphAxes, [0 (globalYMax + 0.2)]);
xlabel(handles.GraphAxes, 'Frames');
set(handles.GraphAxes, 'Visible', 'on');
ylabel(handles.GraphAxes, 'Nuclear to Cytosolic Fluorescence Ratio');
set(handles.GraphAxes, 'FontSize', 16);
box on;
guidata(hObject, handles);

%Initialize graph view to visualize all colonies
function colonyColors = initGraph(hObject, handles)
if isfield(handles,'colonies')
    %originally set graphvis to view all colonies
    handles.coloniesToPlot = 1 : length(handles.colonies);
    colonyColors = helpPlot(hObject, handles, []);
end
if (length(handles.colonies) > 1 && ~handles.withoutColonyInfo)
    %now that you have plotted all colonies you can select specific colonies if
    %more than one colony exists
    set(handles.selectColoniesButton, 'Visible', 'on');
end
guidata(hObject,handles);

function maxIntensityProjections = getMaxIntensityProjectionsNoColonyHelper(images, hObject, handles, totalTimepts, maxIntensityProjections, numberOfChannels)
for channel = 1 : numberOfChannels
    fileNameCurrentChannel = getAndorFileName(images, handles.TIFMicroscopePosition, 0, 0, images.w(channel));
    startIdx = regexp(fileNameCurrentChannel, 'f{1}\d{4}');
    fileNameCurrentChannel(startIdx) = 'p'; %TODO: sloppy mutation - should fix
    structImgs = bfopen(fileNameCurrentChannel);
    for timept = 1 : size(structImgs{1}, 1)
        if channel == 1
            maxIntensityProjections{timept} = ...
                zeros(1024, 1024, numberOfChannels);
        end
        maxIntensityProjections{timept}(:, :, channel) = structImgs{1}{timept};
    end
end

function maxIntensityProjections = getMaxIntensityProjectionsColonyHelper(images, hObject, handles, totalTimepts, maxIntensityProjections, numberOfChannels)
for timept = 1 : totalTimepts
    %TODO: currently assumes all channels and timepoints are the
    %same dimentions 1024 X 1024 though could use metadata
    maxIntensityProjections{timept} = ...
        zeros(1024, 1024, numberOfChannels);
    for channel = 1 : numberOfChannels
        %get max intensity z slice
        currentChannel = ...
            andorMaxIntensityBF(images, handles.TIFMicroscopePosition, timept - 1, images.w(channel));
        maxIntensityProjections{timept}(:, :, channel) = currentChannel;
    end
end

%Generates the max intensity projection for each channel once and stores if
%necessary (depends on MAT file input); currently will only do this for
%colony segregated data
function maxIntensityProjections = getMaxIntensityProjections(images, hObject, handles, totalTimepts)
%TODO: parallelize using parfor()
%TODO: incorporate all time slices
%Init max intensity projection storage for all timepoints for all channels
maxIntensityProjections = cell(1, totalTimepts);
numberOfChannels = length(images.w);
if ~handles.withoutColonyInfo
    maxIntensityProjections = ...
        getMaxIntensityProjectionsColonyHelper(images, hObject, handles, totalTimepts, maxIntensityProjections, numberOfChannels);
else
    %No Colony Data means max intensity calculations are not necessary
    maxIntensityProjections = ...
        getMaxIntensityProjectionsNoColonyHelper(images, hObject, handles, totalTimepts, maxIntensityProjections, numberOfChannels);
end

%Initializes Image Viewer Panel based on user input
function [maxIntensityProjections, contrastLims, totalTimepts]  = initImageViewer(hObject, eventdata, handles)
axes(handles.ImageAxes);
images = readAndorDirectory(handles.imageDirectory);
fileNameCurrentChannel = getAndorFileName(images, handles.TIFMicroscopePosition, 0, 0, images.w(1)); %assume that there is at least one channel
structImgs = bfopen(fileNameCurrentChannel);
totalTimepts = size(structImgs{1}, 1);
maxIntensityProjections = getMaxIntensityProjections(images, hObject, handles, totalTimepts);
handles.maxIntensityProjections = maxIntensityProjections;
contrastLims = [];
handles.contrastLims = contrastLims;
set(handles.scrollImageSet, 'Min', 1);
set(handles.scrollImageSet, 'Max', length(maxIntensityProjections));
set(handles.scrollImageSet, 'Value', 1); %init at timept 1
set(handles.scrollImageSet, 'SliderStep', [1 / (length(maxIntensityProjections)) ,...
    10 / (length(maxIntensityProjections)) ]);
Nuclear_Callback(hObject, eventdata, handles);

% --- Executes on button press in loadTIFFILE. Occurs after loading .mat
% file
function loadTIFFILE_Callback(hObject, eventdata, handles)
% hObject    handle to loadTIFFILE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image.Directory = { {'uigetdir(''.'')'} };
image.MicroscopePosition = 0;
image.NoColonyInfo = 0;
image = StructDlg(image);
handles.TIFMicroscopePosition = image.MicroscopePosition;
handles.imgTimept = 1; %arbitrarily picked to start time from 1
handles.imageDirectory = image.Directory;
%TODO: NOT SURE HOW TO GET IF NO .MAT FILE PROVIDED
handles.withoutColonyInfo = image.NoColonyInfo;
[maxIntensityProjections, contrastLims, totalTimepts] = initImageViewer(hObject, eventdata, handles);
handles.totalTimepts = totalTimepts;
handles.contrastLims = contrastLims; %default value
handles.maxIntensityProjections = maxIntensityProjections;
handles.coloniesToPlot = [];
set(handles.Contrast, 'Visible', 'on');
set(handles.ContrastAdjSlider, 'Visible', 'on');
set(handles.scrollImageSet, 'Visible', 'on');
set(handles.ChannelButtonGroup, 'Visible', 'on');
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

%Update Image Display based on channel selection
handles.imgTimept = int32(get(handles.scrollImageSet, 'Value'));
hObjectContrastAdjSlider = handles.ContrastAdjSlider;
ContrastAdjSlider_Callback(hObjectContrastAdjSlider, eventdata, handles)
delete(allchild(handles.GraphAxes));
if handles.loadedMAT
    if isempty(handles.coloniesToPlot)
        for colony = 1 : length(handles.colonies)
            cellsToEnlarge = getXYColonyTime(handles, colony);
            helpPlot(hObject, handles, colony);
            if handles.withoutColonyInfo
                markerSz = 20;
            else
                markerSz = 50;
            end
            
            plot(handles.GraphAxes, cellsToEnlarge(:, 1), cellsToEnlarge(:, end), 'LineStyle', 'none', ...
                'Marker', '.', 'Markersize', markerSz, 'Color', handles.colonyColors(colony, :));
            hold on;
        end
    else
        cellsToEnlarge = getXYColonyTime(handles, handles.coloniesToPlot);
        helpPlot(hObject, handles, handles.coloniesToPlot);
        plot(handles.GraphAxes, cellsToEnlarge(:, 1), cellsToEnlarge(:, end), 'LineStyle', 'none', ...
            'Marker', '.', 'Markersize', 50, 'Color', handles.colonyColors(handles.coloniesToPlot, :));
        hold on;
    end
    %enlarge points at this time
    handles.cellsToEnlarge = cellsToEnlarge;
end
updateImageDisplay(handles)
guidata(hObject, handles);

% --- Executes on button press in Nuclear.
function Nuclear_Callback(hObject, eventdata, handles)
% hObject    handle to Nuclear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(allchild(handles.ImageAxes));
set(handles.Nuclear, 'Value', get(handles.Nuclear, 'Value'));
updateImageDisplay(handles);
guidata(hObject, handles);

% --- Executes on button press in Cytosolic.
function Cytosolic_Callback(hObject, eventdata, handles)
% hObject    handle to Cytosolic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(allchild(handles.ImageAxes));
set(handles.Cytosolic, 'Value', get(handles.Cytosolic, 'Value'));
updateImageDisplay(handles);
guidata(hObject, handles);

% --- Executes on button press in CellLabelling. Updates ImageAxes and
% Calls Fxn to handle GraphAxes Update.
function CellLabelling_Callback(hObject, eventdata, handles)
% hObject    handle to CellLabelling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(allchild(handles.ImageAxes));
set(handles.CellLabelling, 'Value', get(handles.CellLabelling, 'Value'));
updateImageDisplay(handles);
guidata(hObject, handles);

%Computes x, y coordinates of the cells based on .mat file input and is
%called when cell labelling is selected in slider
function cellInfo = getXYColonyTime(handles, colony)
%Inputs - handles is a structure that is globally visable to all functions
%in the GUI; colony is a positive integer corresponding to the selected colony
%Output - returns time, x, y coordinates to plot, and flurorescence data for a given
%colony if colonies exist otherwise everything is one giant colony

%TODO: redundant calculations here!!! cache calculated values
if handles.withoutColonyInfo
    timesAndCellData = cat(1, handles.colonies.cellData);
    coordinates = cat(1, timesAndCellData(handles.imgTimept).XY);
    someValNAN = (timesAndCellData(handles.imgTimept).nucLevel -...
        timesAndCellData(handles.imgTimept).background) ./ (timesAndCellData(handles.imgTimept).cytLevel - timesAndCellData(handles.imgTimept).background);
    nanidx = isnan(someValNAN);
    allSignal = someValNAN(nanidx ~= 1);
    dimCoord = coordinates(nanidx == 0, :);
    times = ones(1, length(allSignal)) .* double(handles.imgTimept);
    cellInfo = [times' dimCoord allSignal];% AN
else
    times = cat(1, handles.colonies(colony).cells.onframes);
    coordinates = cat(1, handles.colonies(colony).cells.position);
    allSignal = cat(1, handles.colonies(colony).cells.fluorData);
    signalData = allSignal(:, 2) ./ allSignal(:, 3);
    search = [times, coordinates, signalData];% AN
    cellsToPlot = find(search(:, 1) == (handles.imgTimept));
    cellInfo = search(cellsToPlot, :);
end

%Function updates cell labelling
function updateImageDisplayHelper(handles)
if handles.withoutColonyInfo
    cellInfo = getXYColonyTime(handles, 1);
    cellLabels = plot(handles.ImageAxes, cellInfo(:, 2), cellInfo(:, 3), 'LineStyle', 'none',...
        'Marker', '.', 'Markersize', 15, 'Color', handles.colonyColors(1, :));
    text(2, (- 1 / 15 + 0.25), ['Total Cells Found: ' int2str(length(cellInfo(:, 2)))],...
        'Color', handles.colonyColors(1, :), 'Parent', handles.GraphAxes);
    set(cellLabels, 'ButtonDownFcn', @(src, evnt)clickdisplay(src, evnt, handles));
else
    clear colony;
    for colony = 1 : length(handles.colonies)
        cellInfo = getXYColonyTime(handles, colony);
        cellLabels = plot(handles.ImageAxes, cellInfo(:, 2), cellInfo(:, 3), 'LineStyle', 'none',...
            'Marker', '.', 'Markersize', 30, 'Color', handles.colonyColors(colony, :));
        set(cellLabels, 'ButtonDownFcn', @(src, evnt)clickdisplay(src, evnt, handles));
    end
end

%Function Updates Image Display by checking UI
function updateImageDisplay(handles)
nuc = handles.maxIntensityProjections{handles.imgTimept}(:, :, 1);
cyto = handles.maxIntensityProjections{handles.imgTimept}(:, :, 2);
if get(handles.Nuclear, 'Value') && get(handles.Cytosolic, 'Value')
    displayedImg = imshow(imadjust(imfuse(nuc, cyto, 'ColorChannels', ...
        [1 2 0]), handles.contrastLims), 'Parent', handles.ImageAxes);
    
elseif get(handles.Nuclear, 'Value') && ~get(handles.Cytosolic, 'Value')
    displayedImg = imshow(imadjust(imfuse(nuc, nuc, 'ColorChannels', [1 2 0]),...
        handles.contrastLims), 'Parent', handles.ImageAxes);
    
elseif ~get(handles.Nuclear, 'Value') && get(handles.Cytosolic, 'Value')
    displayedImg = imshow(imadjust(imfuse(cyto, cyto, 'ColorChannels', [1 2 0]), ...
        handles.contrastLims), 'Parent', handles.ImageAxes);
else
    %EDGE CASE: ARBITRARILY SET NUCLEAR TO TRUE
    set(handles.Nuclear, 'Value', 1);
    displayedImg = imshow(imadjust(imfuse(nuc, nuc, 'ColorChannels', [1 2 0]),...
        handles.contrastLims), 'Parent', handles.ImageAxes);
end
hold on;
set(displayedImg, 'ButtonDownFcn', @(src, evnt)clickdisplay(src, evnt, handles));

%displays cell labels on ImageAxes
if get(handles.CellLabelling, 'Value')
    if handles.loadedMAT  %check if there is a relevant "lookup table" a.k.a. the .mat file
        updateImageDisplayHelper(handles);
    else
        errordlg('Please Load MATFILE to View Image Labels', 'Missing File Error');
    end
end

% --- Executes on button press in selectColoniesButton.
%Creates StructDlg for selection of specific colonies to plot and modifies
%the GRAPH AXES depending on desired colony visualization
function selectColoniesButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectColoniesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(allchild(handles.GraphAxes));
col = size(handles.colonies, 2);
userSelectedColonies.SelectedColonies = 0; % have the box empty before selection of colony to plot
userSelectedColonies = StructDlg(userSelectedColonies, 'Select Colonies');
colN = userSelectedColonies.SelectedColonies;
if colN > length(handles.colonies)
    errordlg('Invalid Colony', 'Select Valid Colony');
else
    helpPlot(hObject, handles, colN);
    handles.coloniesToPlot = userSelectedColonies.SelectedColonies;
    guidata(hObject, handles);
end

%Function updates click for Colony Data
function clickdisplayHelperColony(src, evnt, handles, mouseCurrentLocation)
allData = getXYColonyTime(handles, handles.coloniesToPlot); %call on all colonies
xyCoord = [allData(:, 2) allData(:, 3)];
xyMousePtDist = ipdm(xyCoord, mouseCurrentLocation, 'Subset', 'NearestNeighbor', 'Result', 'Structure');
[~, Idx] = min(xyMousePtDist.distance);
nearestCell = xyCoord(xyMousePtDist.rowindex(Idx), :);
mouseData = allData(ind2sub(size(allData), find(allData(:, 2) == nearestCell(1) & allData(:, 3) == nearestCell(2))), :);
plot(handles.GraphAxes, handles.imgTimept, mouseData(1, 4), ...
    'LineStyle', 'none', 'Marker', 'p', 'MarkerSize', 12, 'LineWidth', 3, ...
    'Color', [0 0 0], 'MarkerSize', 30);
currentColony = handles.colonies(handles.coloniesToPlot);
totalCells = length(currentColony.cells);
for eachCell = 2 : totalCells
    currentCell = currentColony.cells(eachCell);
    cytoToNuclearFluor = currentCell.fluorData(:, 2) ./ ...
        currentCell.fluorData(:, 3);
    if double(ismember(mouseData(4), cytoToNuclearFluor)) ~= 0
        plot(handles.GraphAxes, currentCell.onframes, cytoToNuclearFluor,...
            '-o', 'Color', 'k', 'LineWidth', 3, 'MarkerSize', 5,...
            'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'k');
        %display the "unique" tag (a number) of that cell that may or may
        %not be accurate
        displayedCellNo = text(nearestCell(1) + 12, nearestCell(2) - 5, int2str(eachCell), ...
            'Color', handles.colonyColors(handles.coloniesToPlot, :), 'FontSize', 14);
        set(displayedCellNo, 'ButtonDownFcn', @(src, evnt)clickdisplay(src, evnt, handles));
        hold on;
    else
        continue;
    end
end

function clickdisplayHelperNoColony(src, evnt, handles, mouseCurrentLocation)
allData = getXYColonyTime(handles, 1); %call on all colonies
xyCoord = [allData(:, 2) allData(:, 3)];
xyMousePtDist = ipdm(xyCoord, mouseCurrentLocation, 'Subset', 'NearestNeighbor', 'Result', 'Structure');
[~, Idx] = min(xyMousePtDist.distance);
nearestCell = xyCoord(xyMousePtDist.rowindex(Idx), :);
mouseData = allData(ind2sub(size(allData), find(allData(:, 2) == nearestCell(1) & allData(:, 3) == nearestCell(2))), :);
plot(handles.GraphAxes, handles.imgTimept, mouseData(1, 4), ...
    'LineStyle', 'none', 'Marker', 'p', 'MarkerSize', 12, 'LineWidth', 3, ...
    'Color', [0 0 0], 'MarkerSize', 30);

%Function that updates the GRAPH AXES an IMAGE depending on the selected
%cell (based on user input)
function clickdisplay(src, evnt, handles)
%Delete old cell marker(s) if it(they) exist(s)
delete(findobj(handles.GraphAxes.Children, 'Marker', 'pentagram'));
delete(findobj(handles.ImageAxes.Children, 'Type', 'Text'));
delete(findobj(handles.GraphAxes.Children, 'LineWidth', 3));
handles = guidata(src);
cursorPoint = get(handles.ImageAxes, 'CurrentPoint');
mouseCurrentLocation = [cursorPoint(1, 1) cursorPoint(1, 2)];
if get(handles.CellLabelling, 'Value') && ~isempty(handles.coloniesToPlot) && ~handles.withoutColonyInfo
    clickdisplayHelperColony(src, evnt, handles, mouseCurrentLocation)
elseif ~isempty(handles.coloniesToPlot) && ~handles.withoutColonyInfo
    errordlg('Please Select A Specific Colony', 'Missing Information');
else %handles.withoutColonyInfo is TRUE
    clickdisplayHelperNoColony(src, evnt, handles, mouseCurrentLocation)
end

function GraphAxes_CreateFcn(hObject, eventdata, handles)

% --- Executes on slider movement.
function ContrastAdjSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ContrastAdjSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
delete(allchild(handles.ImageAxes));
updateImageDisplay(handles);
if (get(hObject, 'Value') / 1000) >  (1 - (get(hObject, 'Value') / 10))
    handles.contrastLims = [(1 - (get(hObject, 'Value') / 10))...
        get(hObject, 'Value') / 1000];
else
    handles.contrastLims = [(get(hObject, 'Value') / 1000)...
        (1 - (get(hObject, 'Value') / 10))];
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ContrastAdjSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ContrastAdjSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor', [.9 .9 .9]);
end
set(hObject, 'Min', 0.5);
set(hObject, 'Max', 10);
set(hObject, 'Value', 1); %init at no contrast
set(hObject, 'SliderStep', [1, 1] ./ (hObject.Max - hObject.Min));
