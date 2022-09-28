function varargout = wapi_gui(varargin)
%wapi_gui               Graphical interface for WAPI
%
% wapi_gui opens the graphical user interface for setting up and starting
% wavelet aided parametric imaging (WAPI). It executes wapi_wapi which
% actually performs the calculations.
%
% See also: WAPI_WAPI

%
%     Copyright (C) 2022 by Zsolt Cselnyi
%
%     The WAPI toolbox (collection of functions listed under the heading
%     "Proper WAPI toolbox functions (toolbox manifest)" in wapi_help.m
%     file) is free software: you can redistribute it and/or modify it
%     under the terms of the GNU General Public License as published by the
%     Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%     Author: Zsolt Cselnyi
%     e-mail: zsolt.cselenyi@ki.se
%
%     Version 2022-08-28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wapi_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @wapi_gui_OutputFcn, ...
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

% --- Executes just before wapi_gui is made visible.
function wapi_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wapi_gui (see VARARGIN)

% Choose default command line output for wapi_gui
if isfield(handles,'initialized')
    return; % was already started
end

handles.output = hObject;

names=fieldnames(handles);
wapiHelp=help('wapi_wapi');

for n=1:length(names)
    name=names{n};
    found=0;
    if ~isempty(strfind(name,'Label'))
        argName=strrep(name,'Label','');
        found=1;
    end
    if ~found && ~isempty(strfind(name,'Checkbox'))
        argName=strrep(name,'Checkbox','');
        found=2;
    end
    if ~found
        continue;
    end
    if found==1 && any(strcmp(names,[argName 'Edit']))
        argControl=[argName 'Edit'];
    elseif found==2
        argControl=name;
    else
        continue;
    end
    toks=regexp(wapiHelp,sprintf('%s\\s.+?</a>:\\s+(\\S.*?)\\n \\n', argName),'tokens','once');
    if isempty(toks)
        waitfor(errordlg(sprintf('No help found for argument %s in wapi_wapi!',argName),'WAPI','modal'));
        continue;
    end
    helpStr=toks{1};
    helpStr=regexprep(helpStr,'\n\s+','\n');
    helpStr=regexprep(helpStr,'<a href="matlab: help \w+">','');
    helpStr=strrep(helpStr,'</a>','');
    set(handles.(argControl),'TooltipString',helpStr);
end 

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes wapi_gui wait for user response (see UIRESUME)
% uiwait(handles.wapi_gui_figure);


% --- Outputs from this function are returned to the command line.
function varargout = wapi_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function infnameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to infnameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function infnameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to infnameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of infnameEdit as text
%        str2double(get(hObject,'String')) returns contents of infnameEdit as a double
infname=get(hObject,'String');

handles=infnameProcess(handles,infname);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tmsnameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmsnameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tmsnameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to tmsnameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmsnameEdit as text
%        str2double(get(hObject,'String')) returns contents of tmsnameEdit as a double
handles=checkTMS(handles);

% Update handles structure
guidata(hObject, handles);


function handles=checkTMS(handles)

tmsname=get(handles.tmsnameEdit,'String');

if handles.tmsrequired 
    set(handles.tmsnameLabel,'Visible','on');
    set(handles.tmsnameEdit,'Visible','on');
    set(handles.tmsnameBrowse,'Visible','on');
    if ~exist(tmsname,'file')
        set(handles.tmsnameLabel,'ForegroundColor',[1 0 0]);
    else
        set(handles.tmsnameLabel,'ForegroundColor',[0 0 0]);
    end
else
    set(handles.tmsnameLabel,'Visible','off');
    set(handles.tmsnameEdit,'Visible','off');
    set(handles.tmsnameBrowse,'Visible','off');
end

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doProcess(hObject, eventdata, handles, false);

function doProcess(hObject, eventdata, handles, compileOnly)
if nargin<4
	compileOnly=false;
end
infname=get(handles.infnameEdit,'String');
[pathname,n,ext]=fileparts(infname);
if ~exist(infname,'file')
	missing=true;
	if strcmp(ext,'.nis') && ~isempty(dir(fullfile(pathname, [n '_*.nii'])))
		missing=false;
	end
	if missing
	    errordlg(sprintf('%s missing required argument',get(handles.infnameLabel,'String')),'WAPI','modal');
		return;
	end
end
if handles.tmsrequired
    tmsname=get(handles.tmsnameEdit,'String');
    if ~exist(tmsname,'file')
        errordlg(sprintf('%s missing required argument',get(handles.tmsnameLabel,'String')),'WAPI','modal');
        return;
    end
else
    tmsname='';
end
reffname=get(handles.reffnameEdit,'String');
if ~exist(reffname,'file')
    errordlg(sprintf('%s missing required argument',get(handles.reffnameLabel,'String')),'WAPI','modal');
    return;
end
ncoeff=str2double(get(handles.ncoeffEdit,'String'));
depth=str2double(get(handles.depthEdit,'String'));
stationary=get(handles.stationaryCheckbox,'Value');
numpoints=str2double(get(handles.numpointsEdit,'String'));
outfname=get(handles.outfnameEdit,'String');
if isempty(outfname)
    errordlg(sprintf('%s missing required argument',get(handles.outfnameLabel,'String')),'WAPI','modal');
    return;
end
refmaskfile=get(handles.refmaskfileEdit,'String');
if ~isempty(refmaskfile) && ~exist(refmaskfile,'file')
    errordlg(sprintf('%s optional argument points to missing file',get(handles.refmaskfileLabel,'String')),'WAPI','modal');
    return;
end
weights=get(handles.weightsEdit,'String');
if ~isempty(weights) && ~exist(weights,'file')
    errordlg(sprintf('%s optional argument points to missing file',get(handles.weightsLabel,'String')),'WAPI','modal');
    return;
end
k2ref=str2double(get(handles.k2refEdit,'String'));
if isnan(k2ref)
    k2ref=[];
end
limsStr=get(handles.limsEdit,'String');
lims=str2num(limsStr); %#ok<ST2NM>
targetmaskfile=get(handles.targetmaskfileEdit,'String');
if ~isempty(targetmaskfile) && ~exist(targetmaskfile,'file')
    errordlg(sprintf('%s optional argument points to missing file',get(handles.targetmaskfileLabel,'String')),'WAPI','modal');
    return;
end
targetmaskpadding=str2double(get(handles.targetmaskpaddingEdit,'String'));
if isnan(targetmaskpadding)
    targetmaskpadding=0;
end
deleteoutputwd3=get(handles.deleteoutputwd3Checkbox,'Value')>0;
opts=struct('refmaskfile',refmaskfile,'weights',weights,'k2ref',k2ref,'lims',lims,'targetmaskfile',targetmaskfile,'targetmaskpadding',targetmaskpadding,'deleteoutputwd3',deleteoutputwd3);

if compileOnly
	try
		cmd=sprintf('outfnames=wapi_wapi(''%s'',''%s'',''%s'',%d,%d,%d,''%s'',%d,%s);',infname,tmsname,reffname,ncoeff,depth,stationary,outfname,numpoints,wapi_struct2str(opts,0));
	catch E
		rethrow(E);
	end
	clipboard('copy',cmd);
	msgbox(sprintf('WAPI command was successfully compiled and copied to the clipboard'), 'WAPI','modal');
else
	try
		set(handles.calculate,'Enable','off');
		tic;
		outfnames=wapi_wapi(infname,tmsname,reffname,ncoeff,depth,stationary,outfname,numpoints,opts);
		t=toc;
	catch E
		set(handles.calculate,'Enable','on');
		rethrow(E);
	end
	set(handles.calculate,'Enable','on');

	hrs=floor(t/60/60);
	trem=t-hrs*60*60;
	mins=floor(trem/60);
	trem=trem-mins*60;
	secs=trem;

	msgbox(sprintf('WAPI successfully finished!\n\nResults saved to:\n%s\nCalculations took:\n%d:%d:%0.3f (%d sec)', ...
		sprintf('%s\n',outfnames{:}), hrs, mins, secs, round(t)), 'WAPI','modal');
end

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);

% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.

set(handles.infnameLabel, 'ForegroundColor', [1 0 0]);
set(handles.infnameEdit, 'String', '');
handles.lastpath='';
set(handles.tmsnameLabel, 'ForegroundColor', [1 0 0]);
set(handles.tmsnameEdit,  'String', '');
set(handles.tmsnameLabel,'Visible','on');
set(handles.tmsnameEdit,'Visible','on');
set(handles.tmsnameBrowse,'Visible','on');
handles.tmsrequired=1;
set(handles.reffnameLabel', 'ForegroundColor', [1 0 0]);
set(handles.reffnameEdit,  'String', '');
set(handles.ncoeffEdit,'String',16);
set(handles.depthEdit,'String',3);
set(handles.numpointsEdit,'String',5);
set(handles.stationaryCheckbox,'Value',1);
set(handles.outfnameLabel,  'ForegroundColor', [1 0 0]);
set(handles.outfnameEdit,  'String', '');
handles.limsStr='';

handles.initialized=true;


% Update handles structure
guidata(handles.wapi_gui_figure, handles);



function refmaskfileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to refmaskfileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refmaskfileEdit as text
%        str2double(get(hObject,'String')) returns contents of refmaskfileEdit as a double


% --- Executes during object creation, after setting all properties.
function refmaskfileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refmaskfileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refmaskfileBrowse.
function refmaskfileBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to refmaskfileBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refmaskfile=get(handles.refmaskfileEdit,'String');
if isempty(refmaskfile)
    lastpath=handles.lastpath;
else
    lastpath=refmaskfile;
end
[filename, pathname] = uigetfile( { ...
    '*.nii', 'NIfTI-1 Images (*.nii)'}, ...
    'Pick an image file', lastpath);
if isequal(filename, 0)
    return;
end
refmaskfile=fullfile(pathname, filename);
set(handles.refmaskfileEdit,'String', refmaskfile);
handles.lastpath = pathname;

% set(handles.reffnameLabel,'ForegroundColor',[0 0 0]);

% Update handles structure
guidata(hObject, handles);




function weightsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to weightsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of weightsEdit as text
%        str2double(get(hObject,'String')) returns contents of weightsEdit as a double


% --- Executes during object creation, after setting all properties.
function weightsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to weightsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in weightsBrowse.
function weightsBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to weightsBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
weights=get(handles.weightsEdit,'String');
if isempty(weights)
    lastpath=handles.lastpath;
else
    lastpath=weights;
end
[filename, pathname] = uigetfile( { ...
    '*.txt', 'Weight factor files (*.txt)'}, ...
    'Pick a weight factor file', lastpath);
if isequal(filename, 0)
    return;
end
weights=fullfile(pathname, filename);
set(handles.weightsEdit,'String', weights);
handles.lastpath = pathname;

% set(handles.weightsLabel,'ForegroundColor',[0 0 0]);

% Update handles structure
guidata(hObject, handles);



function k2refEdit_Callback(hObject, eventdata, handles)
% hObject    handle to k2refEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2refEdit as text
%        str2double(get(hObject,'String')) returns contents of k2refEdit as a double
k2refStr = get(hObject, 'String');
if isempty(k2refStr)
    return;
end

k2ref=str2double(k2refStr);

if isnan(k2ref) || k2ref<0 || k2ref>10 || ~isequal(k2ref,real(k2ref)) || isinf(k2ref)
    set(hObject, 'String', 0);
    errordlg('Input must be a value between 0 and 10','Error');
end


% --- Executes during object creation, after setting all properties.
function k2refEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2refEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function limsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to limsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of limsEdit as text
%        str2double(get(hObject,'String')) returns contents of limsEdit as a double
limsStr=get(hObject,'String');

if isempty(limsStr)
    handles.limsStr=limsStr;
    % Update handles structure
    guidata(hObject, handles);
    return
end

lims=str2num(limsStr); %#ok<ST2NM>
if isempty(lims)
    set(hObject,'String',handles.limsStr);
    return;
end
sz=size(lims);
if ~(isequal(sz,[2 1]) || isequal(sz,[2 3]) || isequal(sz,[2 4]) || isequal(sz,[1 2]))
    set(hObject,'String',handles.limsStr);
    errordlg('Input must be a 2 element vector or a 2x3 or a 2x4 matrix','Error');
else
    handles.limsStr=limsStr;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function limsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to limsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reffnameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to reffnameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reffnameEdit as text
%        str2double(get(hObject,'String')) returns contents of reffnameEdit as a double
reffname=get(hObject,'String');
if ~exist(reffname,'file')
    set(handles.reffnameLabel,'ForegroundColor',[1 0 0]);
else
    set(handles.reffnameLabel,'ForegroundColor',[0 0 0]);
end


% --- Executes during object creation, after setting all properties.
function reffnameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reffnameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reffnameBrowse.
function reffnameBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to reffnameBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
reffname=get(handles.reffnameEdit,'String');
if isempty(reffname)
    lastpath=handles.lastpath;
else
    lastpath=reffname;
end
[filename, pathname] = uigetfile( { ...
    '*.txt;*.blo', 'Input function files (*.txt,*.blo)'; ...
    '*.txt', 'ROI TAC files (*.txt)'; ...
    '*.blo', 'Blood input files (*.blo)'}, ...
    'Pick an input function file', lastpath);
if isequal(filename, 0)
    return;
end
reffname=fullfile(pathname, filename);
set(handles.reffnameEdit,'String', reffname);
handles.lastpath = pathname;

set(handles.reffnameLabel,'ForegroundColor',[0 0 0]);

% Update handles structure
guidata(hObject, handles);




function ncoeffEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ncoeffEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncoeffEdit as text
%        str2double(get(hObject,'String')) returns contents of ncoeffEdit as a double
ncoeff = str2double(get(hObject, 'String'));
if isnan(ncoeff) || ncoeff<4 || ncoeff>128 || ~isequal(ncoeff,round(real(ncoeff))) || isinf(ncoeff)
    set(hObject, 'String', 16);
    errordlg('Input must be an integer between 4 and 128','Error');
end

% --- Executes during object creation, after setting all properties.
function ncoeffEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncoeffEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function depthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to depthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depthEdit as text
%        str2double(get(hObject,'String')) returns contents of depthEdit as a double
val = str2double(get(hObject, 'String'));
if isnan(val) || val<1 || val>8 || ~isequal(val,round(real(val))) || isinf(val)
    set(hObject, 'String', 3);
    errordlg('Input must be an integer between 1 and 8','Error');
end


% --- Executes during object creation, after setting all properties.
function depthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stationaryCheckbox.
function stationaryCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to stationaryCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stationaryCheckbox



function numpointsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to numpointsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numpointsEdit as text
%        str2double(get(hObject,'String')) returns contents of numpointsEdit as a double
val = str2double(get(hObject, 'String'));
if isnan(val) || val<2 || val>200 || ~isequal(val,round(real(val))) || isinf(val)
    set(hObject, 'String', 5);
    errordlg('Input must be an integer between 2 and 200','Error');
end


% --- Executes during object creation, after setting all properties.
function numpointsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numpointsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outfnameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to outfnameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outfnameEdit as text
%        str2double(get(hObject,'String')) returns contents of outfnameEdit as a double
outfname=get(hObject,'String');
handles=outfnameProcess(handles, outfname);

% Update handles structure
guidata(hObject, handles);


function handles=outfnameProcess(handles, outfname)

if nargin>1 % we have a pre-provided outfname, let's look at it
    [pathname,n, ext]=fileparts(outfname);
    handles.lastpath = pathname;
    if ~isempty(outfname) && ~strcmp(ext, '.nii')
        outfname=fullfile(pathname, [n '.nii']);
        set(handles.outfnameEdit, 'String', outfname);
    end
    if isempty(outfname)
        set(handles.outfnameLabel,'ForegroundColor',[1 0 0]);
    else
        set(handles.outfnameLabel,'ForegroundColor',[0 0 0]);
    end
else % no prep-provided infname, let's get one
    outfname=get(handles.outfnameEdit,'String');
    if isempty(outfname)
        lastpath=handles.lastpath;
    else
        lastpath=outfname;
    end
    [filename, pathname] = uiputfile( {'*.nii', 'NIfTI-1 Images (*.nii)'}, 'Specify basename for saving WAPI results', lastpath);
    if isequal(filename, 0)
        return;
    end
    outfname=fullfile(pathname, filename);
    [p,n, ext]=fileparts(outfname);
    if ~strcmp(ext, '.nii')
        outfname=fullfile(pathname, [n '.nii']);
    end
    handles.lastpath = pathname;
    
    set(handles.outfnameEdit,'String', outfname);
    set(handles.outfnameLabel,'ForegroundColor',[0 0 0]); % will not be empty
end


% --- Executes during object creation, after setting all properties.
function outfnameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outfnameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outfnameBrowse.
function outfnameBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to outfnameBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=outfnameProcess(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in infnameBrowse.
function infnameBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to infnameBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=infnameProcess(handles);

% Update handles structure
guidata(hObject, handles);



function handles=infnameProcess(handles,infname)

%        'name', 'nifti', ...
%        'desc', 'NIfTI-1 images', ...
%        'tms', 1);

formats=wapi_readim;

exts=fieldnames(formats);
filterspec={'', 'PET Images ('};
for i=1:length(exts)
    ext=exts{i};
	if strcmp(ext,'nis') % trick is needed, let user select first frame..
		extStr='_001.nii';
	else
		extStr=['.' ext];
	end
	desc=sprintf('%s', formats.(ext).desc);
    if i==1
        filterspec{1,1}=sprintf('*%s', extStr);
        filterspec{1,2}=sprintf('%s*%s', filterspec{1,2}, extStr);
    else
        filterspec{1,1}=sprintf('%s;*%s', filterspec{1,1}, extStr);
        filterspec{1,2}=sprintf('%s,*%s', filterspec{1,2}, extStr);
    end
    filterspec{1+i, 1}=sprintf('*%s', extStr);
    filterspec{1+i, 2}=desc;
end
filterspec{1,2}=[filterspec{1,2} ')'];
    
if nargin>1 % we have a pre-provided infname, let's look at it
    [pathname,n, ext]=fileparts(infname);
    handles.lastpath = pathname;
    if isempty(ext) || ~any(strcmp(exts, ext(2:end)))
        beep;
        errordlg(sprintf('Cannot handle PET images of type %s', ext),'error','modal');
        set(handles.infnameLabel,'ForegroundColor',[1 0 0]);
        handles.tmsrequired=1;
    else
        ext(1)=[];
        handles.tmsrequired=formats.(ext).tms;
    end
    if ~exist(infname,'file')
		if strcmp(ext,'nis')
			if ~exist(fullfile(pathname, [n '_001.nii']),'file')
		        set(handles.infnameLabel,'ForegroundColor',[1 1 0]); % may be ok for .nis virtual filename
			else
		        set(handles.infnameLabel,'ForegroundColor',[0 0 0]); % may be ok for .nis virtual filename
			end
		else
	        set(handles.infnameLabel,'ForegroundColor',[1 0 0]);
		end
    else
        set(handles.infnameLabel,'ForegroundColor',[0 0 0]);
    end
    handles=checkTMS(handles);
else % no prep-provided infname, let's get one
    infname=get(handles.infnameEdit,'String');
    if isempty(infname)
        lastpath=handles.lastpath;
	else
		[p,n, ext]=fileparts(infname);
		if strcmp(ext,'.nis')
			infname=fullfile(p, [n '_001.nii']);
		end
        lastpath=infname;
    end
    [filename, pathname] = uigetfile( filterspec, 'Pick a PET image', lastpath);
    if isequal(filename, 0)
        return;
    end
    infname=fullfile(pathname, filename);
    [p,n, ext]=fileparts(infname);
	if isequal(strfind([n ext],'_001.nii'),length([n ext])-7) % trick into .nis, a 4D image should not be named using _001.nii!!
		infname=fullfile(p, [n(1:end-4) '.nis']);
		ext='.nis';
	end
    if isempty(ext) || ~any(strcmp(exts, ext(2:end)))
        beep;
        errordlg(sprintf('Cannot handle PET images of type %s', ext),'error','modal');
        return;
    end
    ext(1)=[];
    handles.lastpath = pathname;
    
    set(handles.infnameEdit,'String', infname);
    set(handles.infnameLabel,'ForegroundColor',[0 0 0]); % will always exist since uigefile enforces
    
    handles.tmsrequired=formats.(ext).tms;
    handles=checkTMS(handles);
end


% --- Executes on button press in tmsnameBrowse.
function tmsnameBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to tmsnameBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmsname=get(handles.tmsnameEdit,'String');
if isempty(tmsname)
    lastpath=handles.lastpath;
else
    lastpath=tmsname;
end
[filename, pathname] = uigetfile( { ...
    '*.txt;*.tms', 'Frame info files (*.txt,*.tms)'; ...
    '*.txt', 'Frame timing info files (*.txt)'; ...
    '*.tms', 'Frame timing and scaling info files (*.tms)'}, ...
    'Pick a frame info file', lastpath);
if isequal(filename, 0)
    return;
end
tmsname=fullfile(pathname, filename);
set(handles.tmsnameEdit,'String', tmsname);
handles.lastpath = pathname;

handles=checkTMS(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doc wapi_wapi


function targetmaskfileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to targetmaskfileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetmaskfileEdit as text
%        str2double(get(hObject,'String')) returns contents of targetmaskfileEdit as a double


% --- Executes during object creation, after setting all properties.
function targetmaskfileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetmaskfileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in targetmaskfileBrowse.
function targetmaskfileBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to targetmaskfileBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
targetmaskfile=get(handles.targetmaskfileEdit,'String');
if isempty(targetmaskfile)
    lastpath=handles.lastpath;
else
    lastpath=targetmaskfile;
end
[filename, pathname] = uigetfile( { ...
    '*.nii', 'NIfTI-1 Images (*.nii)'}, ...
    'Pick an image file', lastpath);
if isequal(filename, 0)
    return;
end
targetmaskfile=fullfile(pathname, filename);
set(handles.targetmaskfileEdit,'String', targetmaskfile);
handles.lastpath = pathname;

% set(handles.reffnameLabel,'ForegroundColor',[0 0 0]);

% Update handles structure
guidata(hObject, handles);


function targetmaskpaddingEdit_Callback(hObject, eventdata, handles)
% hObject    handle to targetmaskpaddingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetmaskpaddingEdit as text
%        str2double(get(hObject,'String')) returns contents of targetmaskpaddingEdit as a double
targetmaskpaddingStr = get(hObject, 'String');
if isempty(targetmaskpaddingStr)
    return;
end

targetmaskpadding=str2double(targetmaskpaddingStr);

if isnan(targetmaskpadding) || targetmaskpadding<0 || targetmaskpadding>50 || ~isequal(targetmaskpadding,real(targetmaskpadding)) || isinf(targetmaskpadding) || targetmaskpadding~=round(targetmaskpadding)
    set(hObject, 'String', 1);
    errordlg('Input must be an integer value between 0 and 50','Error');
end


% --- Executes during object creation, after setting all properties.
function targetmaskpaddingEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetmaskpaddingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deleteoutputwd3Checkbox.
function deleteoutputwd3Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to deleteoutputwd3Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deleteoutputwd3Checkbox


% --- Executes during object creation, after setting all properties.
function deleteoutputwd3Checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deleteoutputwd3Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in compile.
function compile_Callback(hObject, eventdata, handles)
% hObject    handle to compile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doProcess(hObject, eventdata, handles, true);
