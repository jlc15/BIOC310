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

% Last Modified by GUIDE v2.5 07-Feb-2017 11:20:03

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


% Executes on button press in loadMATFILEButton and constructs Pop-Up window and 
% populates handles structure with .mat analysis information
function loadMATFILEButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadMATFILEButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Output - handle structure with updated information from .mat file for
% graph rendering

%TODO maybe syntax to make empty string
inputWindow.Mat_File = { {'uigetfile(''.'')'} }; 
%TODO: throw exception if user tries to not put in units/etc.
inputWindow.Units = '';

inputWindow = StructDlg(inputWindow, 'Input Plotting Data');
handles.matFileName = inputWindow.Mat_File;
%reg ex matching
handles.units = inputWindow.Units; 
load(handles.matFileName);

%new data has colony information analysis integrated in .mat file
if exist('colonies','var')
    %maybe have different init graph calls here depends on if
    %colonies/cells provided
    handles.colonies = colonies;
elseif exist('cells', 'var')
    handles.cells = cells;
end

%possibly use isLoaded to tell user that the .mat file and the images
%don't match up
set(handles.isLoaded, 'BackgroundColor', [0, 1, 0]);
handles.highlightedCell = 1;
guidata(hObject,handles);
initGraph(hObject, handles);

%Plots the given colonies
function helpPlot(handles)
axes(handles.axes1);
%xlim MAX
globalMaxT = -Inf;
%ylim MAX
globalMaxCytoToNuclearFluor = -Inf;
for colony = 1 : length(handles.coloniesToPlot)
    %check if the current colony is going to be plotted
    if handles.coloniesToPlot(colony)  
        markerColor = rand(3, 1)';
        currentColony = handles.colonies(colony);
        totalCells = length(currentColony.cells);
        for cell = 1 : totalCells
            %get each cell
            currentCell = currentColony.cells(cell);
            cytoToNuclearFluor = currentCell.fluorData(:, 2) ./ ...
                currentCell.fluorData(:, 3);
            currentMaxT = max(currentCell.onframes);
            
            %update global xlim and ylim
            if currentMaxT > globalMaxT
                globalMaxT = currentMaxT;
            end
            
            if max(cytoToNuclearFluor) > globalMaxCytoToNuclearFluor
                globalMaxCytoToNuclearFluor = max(cytoToNuclearFluor);
            end
            plot(currentCell.onframes, cytoToNuclearFluor, '-o', 'MarkerSize', 2, ...
                'MarkerEdgeColor', markerColor, 'MarkerFaceColor', markerColor);
            hold on;
        end
    end
end

%legend('Colony 1', 'Colony 2');
xlim([0 globalMaxT]);
ylim([0 globalMaxCytoToNuclearFluor]);
xlabel(['Time ' handles.units]);
ylabel('Cytosolic to Nuclear Fluorescence Ratio');
set(gca, 'fontsize', 16);
box on;
hold off;

%Initialize graph view to visualize all colonies
function initGraph(hObject, handles)
if isfield(handles,'colonies')
    %originally set graphvis to view all colonies
    handles.coloniesToPlot = ones(1 : length(handles.colonies));
    helpPlot(handles);
end

if (length(handles.colonies) > 1)
    %now that you have plotted all colonies you can select specific colonies if
    %more than one colony exists
    set(handles.selectColoniesButton, 'Visible', 'on');  
end
guidata(hObject,handles);

% --- Executes on button press in loadTIFFILE.
function loadTIFFILE_Callback(hObject, eventdata, handles)
% hObject    handle to loadTIFFILE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image.Directory = { {'uigetdir(''.'')'} };
image.MicroscopePosition = 0;
image = StructDlg(image);
handles.microscopePosition = image.MicroscopePosition;
handles.displayedImageTimept = 0;

handles.imageDirectory = image.Directory;
initImageViewer(handles);
guidata(hObject, handles);


%helper fxn for initImageViewer that generates the max intensity projection
%for each channel
function getMaxIntensityProjections(images, handles)
%assume all channels are the same length
reader = bfGetReader([handles.imageDirectory, '/', images.prefix]);
%handles.maxIntensityProjections = cell(1, reader.getSizeT);
for timept = 1 : reader.getSizeT
    for channel = 1 : length(images.w)
        %get max intensity z slice
        currentChannel = ...
             andorMaxIntensityBF(images, handles.microscopePosition, timept, images.w(channel));
         handles.maxIntensityProjections(channel) = currentChannel;
    end      
end
    
%    


%img0 = ; %diff channels
%     img1 = andorMaxIntensityBF(ff,handles.pos,handles.currtime,1); %diff channels
%     axes(handles.axes1)
%     zz = zeros(size(img0));
%     img2show = {zz,zz,zz};
%     if get(handles.redbox,'Value')
%         img2show{1} = imadjust(img0);
%     end
%     if get(handles.greenbox,'Value')
%         img2show{2} = imadjust(img1);
%     end
%     showImg(img2show); hold on;


%Initializes Image Viewer Panel based on user input
function initImageViewer(handles)
%info about user selected image set
images = readAndorDirectory(handles.imageDirectory); 
% if images.p ~= handles.microscopePosition
%     disp('The desired images are not in this directory');
% end

%get max int projection if necessary
getMaxIntensityProjections(images, handles)
% axes(handles.axes1)
%     zz = zeros(size(img0));
%     img2show = {zz,zz,zz};
%     if get(handles.redbox,'Value')
%         img2show{1} = imadjust(img0);
%     end
%     if get(handles.greenbox,'Value')
%         img2show{2} = imadjust(img1);
%     end
%     showImg(img2show); hold on;
% end

% --- Executes on slider movement.
function scrollImageSet_Callback(hObject, eventdata, handles)
% hObject    handle to scrollImageSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function scrollImageSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scrollImageSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


%otherwise the user had a .mat file of the form 'cells' without 'colonies'


% --- Executes on button press in selectColoniesButton.
%Creates StructDlg for selection of specific colonies to plot
function selectColoniesButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectColoniesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set visability of all colonies to false
colonies = length(handles.coloniesToPlot);
userSelectedColonies.SelectedColonies(colonies) = 0;

%TODO: throw exception if user tries to not put in units/etc.
userSelectedColonies = StructDlg(userSelectedColonies, 'Select Colonies');
handles.coloniesToPlot = userSelectedColonies.SelectedColonies;
guidata(hObject,handles);
updatePlot(handles);

function updatePlot(handles)
helpPlot(handles)


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
