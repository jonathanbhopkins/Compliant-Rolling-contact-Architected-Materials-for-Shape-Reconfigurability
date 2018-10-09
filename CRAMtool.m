function varargout = CRAMtool(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CRAMtool_OpeningFcn, ...
    'gui_OutputFcn',  @CRAMtool_OutputFcn, ...
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

function CRAMtool_OpeningFcn(hObject, eventdata, handles, varargin)
% Front panel image
try
    axes(handles.lab)
    
    data = load('images.mat');
    handles.img_parameter = data.img_parameter;
    logo_Lab = data.img_logo;
    image(logo_Lab)
    
    axis off
    grid off
    set(handles.lab, 'Box', 'off')
    handles.lab.XRuler.Axle.LineStyle = 'none';
    handles.lab.YRuler.Axle.LineStyle = 'none';
catch
    errordlg('images.mat not found')
    return
end

axes(handles.axes)

% Simulation parameters
handles.num = 1e15; % starting force scale
handles.derivative_factor = 1e7; % Derivative factor, i.e. incdir = fabsize/derivative_factor
handles.targetnum = 10;
% Tensile Young's Modulus
handles.Et=0.55e9;
% Compressive Young's Modulus
handles.Ec=0.27e9;
% Tensile Yield Stress
handles.sigmat=27e6;
% Compressive Yield Stress
handles.sigmac=24e6;
% Density
handles.dens=2159;
% Poisson's Ratio
handles.v=0.46;
% Static Coefficient of Friction
handles.fric=0.07;

handles.defaults = [];
handles.defaultsSet = false;

% GUI operating parameters
handles.CRAMlocs = [];
handles.LocMat = [];
handles.Grounds = [];
handles.groundvec = [];
handles.Connections = [];
handles.PitchRadius = [];
handles.PitchCRAM = [];

handles.nf = [];
handles.fmag = [];
handles.dirMtx = [];

handles.nm = [];
handles.mmag = [];

handles.Frames = str2double(get(handles.step7_editFrames, 'String'));

% GUI display parameters
handles.step = 1;

handles.shadowOn = [0 0.447 0.741];
handles.shadowOff = [0.7 0.7 0.7];
handles.doneColor = [0 0.44 0.74];
handles.previousColor = [0.85 0.32 0.1];
handles.ForceColor = [0 0 1];
handles.MomentColor = [0 0 0];
handles.GravColor = [0 0.44 0.74];
handles.LineWidth = 2;
handles.LastGIF = [];

handles.text_force = [];
handles.text_moment = [];

handles.GIF.delay = 0.05;
handles.GIF.loopCount = inf;

set(handles.step1_done, 'ForegroundColor', handles.doneColor)
set(handles.step2_done, 'ForegroundColor', handles.doneColor)
set(handles.step3_done, 'ForegroundColor', handles.doneColor)
set(handles.step4_done, 'ForegroundColor', handles.doneColor)
set(handles.step5_done, 'ForegroundColor', handles.doneColor)
set(handles.step6_done, 'ForegroundColor', handles.doneColor)
set(handles.step7_done, 'ForegroundColor', handles.doneColor)
set(handles.step2_previous, 'ForegroundColor', handles.previousColor)
set(handles.step3_previous, 'ForegroundColor', handles.previousColor)
set(handles.step4_previous, 'ForegroundColor', handles.previousColor)
set(handles.step5_previous, 'ForegroundColor', handles.previousColor)
set(handles.step6_previous, 'ForegroundColor', handles.previousColor)
set(handles.step7_previous, 'ForegroundColor', handles.previousColor)
set(handles.step8_previous, 'ForegroundColor', handles.previousColor)

UpdateAxes(handles);
handles = UpdateDisplay(hObject, handles, handles.step);

% This code allows key presses (like 'l' or 'g' to be pressed even if the
% main GUI is not selected (such as when a button is in focus)
p = get(hObject, 'Children');
for k = 2:length(p)
    pp = get(p(k),'Children');
    for j = 1:length(pp)
        if ~isempty(get(pp(j),'Tag'))
            if isprop(pp(j),'KeyPressFcn')
                tag = get(pp(j),'Tag');
                str = ['handles. ' tag];
                fcnstr = '@(hObject,eventdata)CRAMtool(''figure1_KeyPressFcn'',hObject,eventdata,guidata(hObject))';
                set(eval(str),'KeyPressFcn',eval(fcnstr))
            end
            if isprop(pp(j), 'Units')
                tag = get(pp(j),'Tag');
                str = ['handles. ' tag];
                set(eval(str),'Units','pixels')
            end
        end
    end
end
set(handles.view_x,'KeyPressFcn',[])
set(handles.view_y,'KeyPressFcn',[])
set(handles.view_d,'KeyPressFcn',[])
set(handles.step7_editFrames,'KeyPressFcn',[])

% Choose default command line output for CRAMtool
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

function varargout = CRAMtool_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

%% GUI operation
function handles = simulation_settings_ClickedCallback(hObject, eventdata, handles)
defaults{1} = num2str(handles.derivative_factor);
defaults{2} = num2str(handles.targetnum);
str = {'Derivative factor: Make this number as large as possible without causing a singular matrix error.', ...
    'Step-size factor: The larger this number is, the more accurate the simulation will be but the longer it will take to run.'};
Answer = inputdlg(str, 'Simulation settings', 1, defaults);
if isempty(Answer); return; end

handles.derivative_factor = str2double(Answer{1});
handles.targetnum = str2double(Answer{2});

guidata(hObject, handles)
function menu_saveCRAM_ClickedCallback(hObject, eventdata, handles)
[filename, pathname] = uiputfile('*.mat', 'Pick .mat file');
if isequal(filename, 0) || isequal(pathname, 0) % Cancelled
    return
end
savefile = fullfile(pathname, filename);

if ~isempty(handles.CRAMlocs) && ~isempty(handles.Connections)
    InputLoc = handles.LocMat(:,2:3);
    Grounds = handles.groundvec;
elseif ~isempty(handles.CRAMlocs)
    InputLoc = handles.CRAMlocs;
    Grounds = handles.Grounds;
else
    return
end

InputConnect = handles.Connections;
[FOV(1), FOV(2), FOV(3)] = GetAxes(handles);

save(savefile, 'InputLoc', 'Grounds', 'InputConnect', 'FOV')

function figure1_KeyPressFcn(hObject, eventdata, handles)
key = eventdata.Key;
switch key
    case 'n'
        if handles.step == 1
            handles = step1_newCRAMclick_Callback(hObject, eventdata, handles);
        elseif handles.step == 2
            handles = step2_newGround_Callback(hObject, eventdata, handles);
        elseif handles.step == 3
            handles = step3_newConnection_Callback(hObject, eventdata, handles);
        end
    case 'r'
        if handles.step == 1
            handles = step1_removeCRAM_Callback(hObject, eventdata, handles);
        elseif handles.step == 2
            handles = step2_removeGround_Callback(hObject, eventdata, handles);
        elseif handles.step == 3
            handles = step3_removeConnection_Callback(hObject, eventdata, handles);
        elseif handles.step == 6
            handles = step6_remove_Callback(hObject, eventdata, handles);
        end
    case 'equal'
        value = str2double(get(handles.view_d, 'String'));
        set(handles.view_d, 'String', num2str(value / 1.1))
        UpdateAxes(handles);
    case 'hyphen'
        value = str2double(get(handles.view_d, 'String'));
        set(handles.view_d, 'String', num2str(value * 1.1))
        UpdateAxes(handles);
    case 'uparrow'
        pos = str2double(get(handles.view_y, 'String'));
        range = str2double(get(handles.view_d, 'String'));
        set(handles.view_y, 'String', num2str(pos + 0.1*range))
        UpdateAxes(handles);
    case 'downarrow'
        pos = str2double(get(handles.view_y, 'String'));
        range = str2double(get(handles.view_d, 'String'));
        set(handles.view_y, 'String', num2str(pos - 0.1*range))
        UpdateAxes(handles);
    case 'leftarrow'
        pos = str2double(get(handles.view_x, 'String'));
        range = str2double(get(handles.view_d, 'String'));
        set(handles.view_x, 'String', num2str(pos - 0.1*range))
        UpdateAxes(handles);
    case 'rightarrow'
        pos = str2double(get(handles.view_x, 'String'));
        range = str2double(get(handles.view_d, 'String'));
        set(handles.view_x, 'String', num2str(pos + 0.1*range))
        UpdateAxes(handles);
    case 'return'
        switch handles.step
            case 1; handles = step1_done_Callback(hObject, eventdata, handles);
            case 2; handles = step2_done_Callback(hObject, eventdata, handles);
            case 3; handles = step3_done_Callback(hObject, eventdata, handles);
            case 4; handles = step4_done_Callback(hObject, eventdata, handles);
            case 5; handles = step5_done_Callback(hObject, eventdata, handles);
            case 6; handles = step6_done_Callback(hObject, eventdata, handles);
            case 7; handles = step7_done_Callback(hObject, eventdata, handles);
        end
    case 'backspace'
        switch handles.step
            case 2; handles = step2_previous_Callback(hObject, eventdata, handles);
            case 3; handles = step3_previous_Callback(hObject, eventdata, handles);
            case 4; handles = step4_previous_Callback(hObject, eventdata, handles);
            case 5; handles = step5_previous_Callback(hObject, eventdata, handles);
            case 6; handles = step6_previous_Callback(hObject, eventdata, handles);
            case 7; handles = step7_previous_Callback(hObject, eventdata, handles);
            case 8; handles = step8_previous_Callback(hObject, eventdata, handles);
        end
end
guidata(hObject, handles)

function handles = UpdateInstructions(handles, step)
switch step
    case 1
        str = {'Define design window parameters using "Set axis limits".';
            ' - "Center:" defines the central point in the design window.';
            ' - "Field-of-view (zoom):" defines the range of the axes.';
            ' ';
            'Adjust settings in "1. Add circular cams" to place new cams.';
            ' - Select cam locations using "Click to add". (alt.: n key)';
            ' - Enter exact cam locations using "By value".';
            ' - Upload .mat file of desired locations using "Load .mat".';
            ' - Remove undesired cam locations using "Delete cam."            (alt.: r key) '
            ' - Clear the entire design window using "Clear all".';
            ' - Click "Next" when you have finished placing the cams.'};
    case 2
        str = {'Define which of the previously placed cams are held fixed';
            'using "2. Set grounds". Note that grounded cams are shown';
            'red while all other cams are shown blue.';
            ' - Select "New ground (click)" to define grounded cams.          (alt.: n key)';
            ' - Select "Remove ground" to un-ground a grounded cam.        (alt.: r key)';
            ' - Click "Previous" to return to the last step.';
            ' - Click "Next" when you have finished grounding cams.'};
    case 3
        str = {'Define which cams should be joined together by flexure';
            'straps using "3. Add connections".';
            ' - Click "New connection" and then click on the two cams';
            '   you''d like to join. A dotted black line will designate';
            '   connected cams. (alt.: n key)';
            ' - Click "Remove connection" to delete a joint. (alt.: r key)';
            ' - Click "Previous" to return to the last step.';
            ' - Click "Next" when you have finished joining the cams.'};
    case 4
        str = {'Define the pitch radius of a single circular cam using';
            '"4. Specify pitch radius" and the remaining radii will be';
            'calculated automatically.';
            ' - Click "Specify pitch radius" and then select the cam for';
            '   which you''d like to enter its pitch radius.';
            ' - Enter the pitch radius for that cam in the pop-up window.'
            ' - Click "Previous" to return to the last step.';
            ' - Click "Next" when you have finished defining cam radii.'};
    case 5
        str = {'Define geometric parameters and material properties using';
            '"5. Geometry and materials".';
            ' - Click "Geometry" to specify the geometric parameters.';
            '   They are initially set to reasonable default values that';
            '   can be changed in the pop-up window as desired.';
            '   An image will fill the design window with the labeled';
            '   geometric parameters to clarify what they are.';
            ' - Click "Materials" to specify the material properties';
            '   of the consituent material. The default properties are';
            '   those of Teflon.';
            ' - Click "Previous" to return to the last step.';
            ' - Click "Next" when you are finished entering the values.'};
    case 6
        str = {'Define loads on the ungrounded cams using "6. Add loads".';
            'Note that forces will always act at the cam centers and';
            'will follow the cams as they move.';
            ' - Click "Force" or "Moment" to add a load. Be sure to click';
            '   near the center of the cam that you''d like to load.';
            '   Pop-up windows will guide you in selecting the '
            '   magnitudes and directions of these loads.'
            ' - Click "Remove" to delete every load on the cam that you';
            '   select except for gravity.';
            ' - Click "Gravity?" to load all the ungrounded cams in the';
            '   downward direction using Earth''s gravitational pull.';
            ' - Click "Previous" to return to the last step.';
            ' - Click "Next" when you are finished loading the cams.'};
    case 7
        str = {'You are now ready to simulated your CRAM design for the';
            'loads specified. The simulation will generate an animation of';
            'the CRAM moving. Select the number of frames you''d like ';
            'this animation to contain using "Frames:" in '
            '"7. Run simulation".';
            ' - Click "Simulate" to run the simulation.';
            ' - Click "Stop" to terminate the simulation.';
            ' - Click "Previous" to return to the last step.';
            ' - Click "Next" when you are finished loading the cams.'};
    case 8
        str = {'You can now see how each layer should be CADed with ';
            'straight straps to enable the fabrication and assembly of ';
            'your design using "8. CAD to fabricate".';
            ' - Click "View" to see the undefomed layers.'
            ' - Click "Previous" to return to the last step.'};
end
set(handles.text_instructions, 'String', str)
function handles = UpdateDisplay(hObject, handles, step, clear_all)
% Step 1
set(handles.panel_axes, 'ShadowColor', handles.shadowOff)
set(handles.panel_addCRAMs, 'ShadowColor', handles.shadowOff)
set(handles.step1_loadLocs, 'Enable', 'Off')
set(handles.step1_newCRAMclick, 'Enable', 'Off')
set(handles.step1_newCRAMvalue, 'Enable', 'Off')
set(handles.step1_done, 'Enable', 'Off')
set(handles.step1_removeCRAM, 'Enable', 'Off')
set(handles.step1_clearAll, 'Enable', 'Off')
% Step 2
set(handles.panel_grounds, 'ShadowColor', handles.shadowOff)
set(handles.step2_newGround, 'Enable', 'Off')
set(handles.step2_removeGround, 'Enable', 'Off')
set(handles.step2_done, 'Enable', 'Off')
set(handles.step2_previous, 'Enable', 'Off')
% Step 3
set(handles.panel_connections, 'ShadowColor', handles.shadowOff)
set(handles.step3_newConnection, 'Enable', 'Off')
set(handles.step3_removeConnection, 'Enable', 'Off')
set(handles.step3_done, 'Enable', 'Off')
set(handles.step3_previous, 'Enable', 'Off')
% Step 4
set(handles.panel_matPropSettings, 'ShadowColor', handles.shadowOff)
set(handles.step4_geometrySettings, 'Enable', 'Off')
set(handles.step4_materialSettings, 'Enable', 'Off')
set(handles.step4_done, 'Enable', 'Off')
set(handles.step4_previous, 'Enable', 'Off')
% Step 5
set(handles.panel_settings, 'ShadowColor', handles.shadowOff)
set(handles.step5_pitchRadius, 'Enable', 'Off')
set(handles.step5_done, 'Enable', 'Off')
set(handles.step5_previous, 'Enable', 'Off')
% Step 6
set(handles.panel_forces, 'ShadowColor', handles.shadowOff)
set(handles.step6_done, 'Enable', 'Off')
set(handles.step6_previous, 'Enable', 'Off')
set(handles.step6_addForce, 'Enable', 'Off')
set(handles.step6_addMoment, 'Enable', 'Off')
set(handles.step6_remove, 'Enable', 'Off')
set(handles.step6_gravity, 'Enable', 'Off')
% Step 7
set(handles.panel_simulate, 'ShadowColor', handles.shadowOff)
set(handles.step7_simulate, 'Enable', 'Off')
set(handles.step7_editFrames, 'Enable', 'Off')
set(handles.step7_previous, 'Enable', 'Off')
set(handles.step7_stop, 'Enable', 'Off')
set(handles.step7_done, 'Enable', 'Off')
% Step 8
set(handles.step8_undeformed, 'Enable', 'Off')
set(handles.panel_fabricate, 'ShadowColor', handles.shadowOff)
set(handles.step8_previous, 'Enable', 'Off')

if nargin < 4
    switch step
        case 1
            set(handles.panel_axes, 'ShadowColor', handles.shadowOn)
            set(handles.panel_addCRAMs, 'ShadowColor', handles.shadowOn)
            set(handles.step1_loadLocs, 'Enable', 'On')
            set(handles.step1_newCRAMclick, 'Enable', 'On')
            set(handles.step1_newCRAMvalue, 'Enable', 'On')
            set(handles.step1_done, 'Enable', 'On')
            set(handles.step1_removeCRAM, 'Enable', 'On')
            set(handles.text_status, 'String', 'Step 1: Add circular cams by specifying their centroids.')
            set(handles.step1_clearAll, 'Enable', 'On')
        case 2
            set(handles.panel_grounds, 'ShadowColor', handles.shadowOn)
            set(handles.step2_newGround, 'Enable', 'On')
            set(handles.step2_removeGround, 'Enable', 'On')
            set(handles.step2_done, 'Enable', 'On')
            set(handles.step2_previous, 'Enable', 'On')
            set(handles.text_status, 'String', 'Step 2: Select which cams should be grounded or held fixed.')
        case 3
            set(handles.panel_connections, 'ShadowColor', handles.shadowOn)
            set(handles.step3_newConnection, 'Enable', 'On')
            set(handles.step3_removeConnection, 'Enable', 'On')
            set(handles.step3_done, 'Enable', 'On')
            set(handles.step3_previous, 'Enable', 'On')
            set(handles.text_status, 'String', 'Step 3: Define which cams should be connected together by flexure straps.')
        case 4
            set(handles.panel_settings, 'ShadowColor', handles.shadowOn)
            set(handles.step5_pitchRadius, 'Enable', 'On')
            set(handles.step5_done, 'Enable', 'On')
            set(handles.step5_previous, 'Enable', 'On')
            set(handles.text_status, 'String', 'Step 4: Define the pitch radii for every circular cam.')
        case 5
            set(handles.panel_matPropSettings, 'ShadowColor', handles.shadowOn)
            set(handles.step4_geometrySettings, 'Enable', 'On')
            set(handles.step4_materialSettings, 'Enable', 'On')
            set(handles.step4_done, 'Enable', 'On')
            set(handles.step4_previous, 'Enable', 'On')
            set(handles.text_status, 'String', 'Step 5: Input geometric parameters and material properties.')
        case 6
            set(handles.panel_forces, 'ShadowColor', handles.shadowOn)
            set(handles.step6_addForce, 'Enable', 'On')
            set(handles.step6_addMoment, 'Enable', 'On')
            set(handles.step6_remove, 'Enable', 'On')
            set(handles.step6_done, 'Enable', 'On')
            set(handles.step6_previous, 'Enable', 'On')
            set(handles.step6_gravity, 'Enable', 'On')
            set(handles.text_status, 'String', 'Step 6: Specify forces and moments on the desired cams.')
        case 7
            set(handles.panel_simulate, 'ShadowColor', handles.shadowOn)
            set(handles.step7_simulate, 'Enable', 'On')
            set(handles.step7_editFrames, 'Enable', 'On')
            set(handles.step7_previous, 'Enable', 'On')
            set(handles.panel_axes, 'ShadowColor', handles.shadowOn)
            set(handles.text_status, 'String', 'Step 7: Simulate your design for the given loads to generate an animation.')
            set(handles.step7_done, 'Enable', 'On')
        case 8
            set(handles.step8_undeformed, 'Enable', 'On')
            set(handles.panel_fabricate, 'ShadowColor', handles.shadowOn)
            set(handles.text_status, 'String', 'Step 8: View a CAD model of each undeformed layer so your design can be fabricated and assembled.')
            set(handles.step8_previous, 'Enable', 'On')
    end
end

handles = UpdateInstructions(handles, step);
guidata(hObject, handles)
function handles = drawSimple(hObject, handles)
try delete(handles.blue); end
try delete(handles.red); end
try
    for k = 1:length(handles.black)
        delete(handles.black{k})
    end
end
try
    for k = 1:length(handles.circle)
        delete(handles.circle{k})
    end
end

try
    delete(handles.axes.Children)
end

if handles.step >= 3
    x = handles.LocMat(:,2);
    y = handles.LocMat(:,3);
    grounds = handles.groundvec;
    Connections = handles.Connections;
else
    x = handles.CRAMlocs(:,1);
    y = handles.CRAMlocs(:,2);
    grounds = handles.Grounds;
    Connections = handles.Connections;
end

hold off;
handles.blue = plot(handles.axes, x, y, 'bo');
hold on;

if ~isempty(grounds)
    handles.red = plot(handles.axes, x(grounds), y(grounds), 'ro');
end

if ~isempty(Connections)
    for k = 1:size(Connections,1)
        I1 = Connections(k,1);
        I2 = Connections(k,2);
        handles.black{k} = plot(handles.axes, [x(I1), x(I2)], [y(I1), y(I2)], 'k--');
    end
end

UpdateAxes(handles);
guidata(hObject, handles)
function handles = DrawCircle(handles, xc, yc, r)
theta = linspace(0,2*pi, 200)';
x = r*cos(theta) + xc;
y = r*sin(theta) + yc;
handles.circle{1} = plot(handles.axes, x, y, 'k');
function handles = DrawCircles(handles)
for k = 1:size(handles.LocMat,1)
    handles.circle{k} = plaincircle(handles.LocMat(k,2), handles.LocMat(k,3), handles.LocMat(k,4));
end
UpdateAxes(handles);
function handles = drawPatchCRAM(handles)
FullJointMat = handles.FullJointMat;
LocMat = handles.LocMat;
jointnum = handles.jointnum;
C = handles.C;
boltr = handles.boltr;
t = handles.t;
W = handles.W;
del = handles.del;
dcut1 = handles.dcut1;
dcut2 = handles.dcut2;
O = handles.O;
Et = handles.Et;
Ec = handles.Ec;
v = handles.v;
fric = handles.fric;

if ~isempty(handles.text_force)
    for k = 1:length(handles.text_force)
        delete(handles.text_force{k})
    end
    handles.text_force = [];
end
if ~isempty(handles.text_moment)
    for k = 1:length(handles.text_moment)
        delete(handles.text_moment{k})
    end
    handles.text_moment = [];
end

drawCRAM(FullJointMat,LocMat,jointnum, C, boltr, t, W, dcut1, dcut2, Ec, v)
UpdateAxes(handles);
function handles = drawForces(handles)
axes(handles.axes)
hold on
xl = xlim;
yl = ylim;
xr = range(xl);
yr = range(yl);
drawscale = max([xr yr])/20;

if get(handles.step6_gravity, 'Value')
    [~, ~, gravforce] = UpdateGravity(handles);
else
    gravforce = [];
end
maxF = max([max(handles.fmag), max(gravforce)]);

if ~isempty(handles.fmag)
    for k = 1:size(handles.fmag,1)
        if handles.fmag(k) == 0; continue; end % Do not draw mag = 0
        forcedir = handles.dirMtx(k,:);
        x = handles.LocMat(handles.nf(k),2);
        y = handles.LocMat(handles.nf(k),3);
        str = ['F = ' num2str(handles.fmag(k)) ' N'];
        
        scale = handles.fmag(k) / maxF;
        scale = max([scale 0.1]);
        
        xc = forcedir(1)*drawscale*3.5*scale;
        yc = forcedir(2)*drawscale*3.5*scale;
        
        t = text(x+xc, y+yc, str);
        set(t, 'FontSize', 14)
        set(t, 'Color', handles.ForceColor)
        
        ext = get(t, 'Extent');
        h = ext(4);
        if yc < h/2
            yc = yc + sign(yc)*h/2;
            set(t, 'Position', [x+xc y+yc])
        end
        
        theta = atan2(forcedir(2), forcedir(1));
        if cos(theta) > 1/2
            set(t, 'HorizontalAlignment', 'Left')
        elseif cos(theta) < -1/2
            set(t, 'HorizontalAlignment', 'Right')
        else
            set(t, 'HorizontalAlignment', 'Center')
        end
        handles.text_force{k} = t;
        
        q = quiver(x, y, forcedir(1)*drawscale*3*scale, forcedir(2)*drawscale*3*scale);
        len = hypot(forcedir(1)*drawscale*3*scale, forcedir(2)*drawscale*3*scale);
        set(q, 'Color', handles.ForceColor);
        set(q, 'MaxHeadSize', drawscale/len*2)
        set(q, 'MaxHeadSize', 1)
        set(q, 'LineWidth', handles.LineWidth)
    end
end

if ~isempty(handles.mmag)
    for k = 1:size(handles.mmag,1)
        if handles.mmag(k) == 0; continue; end % Do not draw mag = 0
        x = handles.LocMat(handles.nm(k),2);
        y = handles.LocMat(handles.nm(k),3);
        str = ['M = ' num2str(abs(handles.mmag(k))) ' Nm'];
        t = text(x, y-drawscale*1.3, str);
        set(t, 'FontSize', 14)
        set(t, 'Color', handles.MomentColor)
        handles.text_moment{k} = t;
        
        drawCircularArrow(handles, x, y, drawscale/1.5, sign(handles.mmag(k)))
    end
end

% Check to see if gravity toggle is enabled
if get(handles.step6_gravity, 'Value')
    for k = 1:length(gravforce)
        F = gravforce(k);
        
        forcedir = [0 -1];
        x = handles.LocMat(k,2);
        y = handles.LocMat(k,3);
        str = ['F = ' num2str(F,3) ' N'];
        
        scale = F/maxF;
        scale = max([scale 0.1]);
        
        xc = forcedir(1)*drawscale*3.5*scale;
        yc = forcedir(2)*drawscale*3.5*scale;
        
        t = text(x+xc, y+yc, str);
        set(t, 'FontSize', 14)
        set(t, 'Color', handles.GravColor)
        
        ext = get(t, 'Extent');
        h = ext(4);
        if yc < h/2
            yc = yc + sign(yc)*h/2;
            set(t, 'Position', [x+xc y+yc])
        end
        
        set(t, 'HorizontalAlignment', 'Center')
        handles.text_force{k} = t;
        
        q = quiver(x, y, forcedir(1)*drawscale*3*scale, forcedir(2)*drawscale*3*scale);
        len = hypot(forcedir(1)*drawscale*3*scale, forcedir(2)*drawscale*3*scale);
        set(q, 'Color', handles.GravColor);
        set(q, 'MaxHeadSize', drawscale/len*2)
        set(q, 'LineWidth', handles.LineWidth)
    end
end
function drawGrounds(handles)
x = handles.LocMat(:,2);
y = handles.LocMat(:,3);
grounds = handles.groundvec;
if ~isempty(grounds)
    plot(handles.axes, x(grounds), y(grounds), 'ro');
end

function drawCircularArrow(handles, xc, yc, radius, dir)
% Circle
theta = linspace(pi/4,7*pi/4,200);
x = radius*cos(theta) + xc;
y = radius*sin(theta) + yc;
plot(x, y, 'color', handles.MomentColor, 'LineWidth', handles.LineWidth)

% Arrowhead
L = radius/2;
dir = sign(dir);
if dir == 1
    xe = x(end);
    ye = y(end);
    x1 = [xe xe-L];
    y1 = [ye ye];
    x2 = [xe xe];
    y2 = [ye ye-L];
else
    xe = x(1);
    ye = y(1);
    x1 = [xe xe-L];
    y1 = [ye ye];
    x2 = [xe xe];
    y2 = [ye ye+L];
end
plot(x1, y1, 'color', handles.MomentColor, 'LineWidth', handles.LineWidth)
plot(x2, y2, 'color', handles.MomentColor, 'LineWidth', handles.LineWidth)
function radius = SelectionRadius(handles)
dmax = str2double(get(handles.view_d,'String'));
radius = dmax/50;
function [D, Ddiag, Dtri] = DistMatrix(P1, P2)
n1 = size(P1,1);
n2 = size(P2,1);

D = zeros(n1, n2);
for r = 1:n1
    for c = 1:n2
        p1 = P1(r,:);
        p2 = P2(c,:);
        D(r,c) = sqrt( (p1(1)-p2(1))^2 + (p1(2)-p2(2))^2 );
    end
end

if nargout >= 2 && n1 == n2
    Ddiag = D;
    for r = 1:n1
        Ddiag(r,r) = NaN;
    end
end

if nargout > 2 && n1 == n2
    Dtri = triu(D);
    Dtri(Dtri == 0) = NaN;
end
function toolbar_hotkeys_ClickedCallback(hObject, eventdata, handles)
str = {'CRAM Tool hotkeys: ';
    ' ';
    'N: add new cam, ground, or connection (steps 1-3)';
    'R: remove cam, ground, or connection (steps 1-3)';
    ' ';
    'Arrow keys: pan the display area';
    '+/- keys: zoom the display area';
    ' ';
    'Return: proceed to next step';
    'Backspace: return to previous step'};
msgbox(str,'CRAM Tool hotkeys')
function r = range(vec)
r = max(vec) - min(vec);

%% Step 1: Set axes and input CRAM locations
function UpdateAxes(handles)
xc = str2double(get(handles.view_x, 'String'));
yc = str2double(get(handles.view_y, 'String'));
d  = str2double(get(handles.view_d, 'String'));

xMin = xc - d/2;
xMax = xc + d/2;
yMin = yc - d/2;
yMax = yc + d/2;

axis equal
xlim(handles.axes, [xMin xMax])
ylim(handles.axes, [yMin yMax])

xlabel('Meters')
function [xc, yc, d] = GetAxes(handles)
xc = str2double(get(handles.view_x, 'String'));
yc = str2double(get(handles.view_y, 'String'));
d  = str2double(get(handles.view_d, 'String'));
function SetAxes(handles, xc, yc, d)
set(handles.view_x, 'String', num2str(xc))
set(handles.view_y, 'String', num2str(yc))
set(handles.view_d, 'String', num2str(d))
UpdateAxes(handles);
guidata(hObject, handles)

function bool = BoundCheck(handles, x, y)
% Make sure the user clicked in the box
xl = xlim(handles.axes);
yl = ylim(handles.axes);
bx = and(x>xl(1), x<xl(2));
by = and(y>yl(1), y<yl(2));
bool = and(bx, by);
function view_x_Callback(hObject, eventdata, handles)
UpdateAxes(handles);
function view_x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function view_y_Callback(hObject, eventdata, handles)
UpdateAxes(handles);
function view_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function view_d_Callback(hObject, eventdata, handles)
UpdateAxes(handles);
function view_d_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function step1_loadLocs_Callback(hObject, eventdata, handles)
str = {'Use the Save Current CRAM Configuration icon in the GUI toolbar to automatically save a .mat file'};
uiwait(msgbox(str,'Load CRAM','modal'));

[filename, pathname] = uigetfile('*.mat', 'Pick .mat file');
if isequal(filename, 0) || isequal(pathname, 0) % Cancelled
    return
end
data = load(fullfile(pathname, filename));
handles.CRAMlocs = data.InputLoc;
handles.CRAMlocs = bsxfun(@minus, handles.CRAMlocs, mean(handles.CRAMlocs));

xr = range(handles.CRAMlocs(:,1));
yr = range(handles.CRAMlocs(:,2));
r = max([xr yr]);
set(handles.view_d, 'String', num2str(3*r))

set(handles.text_status, 'String', ['Loaded: ' filename ', press Next to continue'])
try
    handles.Grounds = data.Grounds;
catch
    handles.Grounds = [];
end
try
    handles.Connections = data.InputConnect;
catch
    handles.Connections = [];
end
try
    FOV = data.FOV;
    SetAxes(handles, FOV(1), FOV(2), FOV(3));
end

handles = drawSimple(hObject, handles);
guidata(hObject, handles)
function handles = step1_newCRAMclick_Callback(hObject, eventdata, handles)
if ~isempty(handles.Connections)
    Button = questdlg('Adding a new cam will delete existing connections','Connection warning','Proceed','Cancel','Proceed');
    switch Button
        case 'Proceed'
            handles.Connections = [];
        case 'Cancel'
            return
    end
end
handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x, y] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x, y); if ~bool; return; end

handles.CRAMlocs = [x, y; handles.CRAMlocs];
handles.Grounds = handles.Grounds + 1;

handles = drawSimple(hObject, handles);
guidata(hObject, handles)
function step1_newCRAMvalue_Callback(hObject, eventdata, handles)
if ~isempty(handles.Connections)
    Button = questdlg('Adding a new cam will delete existing connections','Connection warning','Proceed','Cancel','Proceed');
    switch Button
        case 'Proceed'
            handles.Connections = [];
        case 'Cancel'
            return
    end
end

Answer = inputdlg({'X', 'Y'}, 'Enter new cam centroid location', 1, {'0', '0'});
if isempty(Answer); return; end
x = str2double(Answer{1});
y = str2double(Answer{2});
% handles.CRAMlocs(end+1,:) = [x, y];
handles.CRAMlocs = [x, y; handles.CRAMlocs];
handles.Grounds = handles.Grounds + 1;
if ~isempty(handles.Connections); handles.Connections = handles.Connections + 1; end

handles = drawSimple(hObject, handles);
guidata(hObject, handles)
function handles = step1_removeCRAM_Callback(hObject, eventdata, handles)
if ~isempty(handles.Connections)
    Button = questdlg('Removing a cam will delete existing connections','Connection warning','Proceed','Cancel','Proceed');
    switch Button
        case 'Proceed'
            handles.Connections = [];
        case 'Cancel'
            return
    end
end

handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x, y] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x, y); if ~bool; return; end
D = DistMatrix([x y], handles.CRAMlocs);
[~,I] = min(D);
if D(I) < SelectionRadius(handles)
    handles.CRAMlocs(I,:) = [];
    
    handles.Grounds = [];
    handles.groundvec = [];
    handles.Connections = [];
    handles.PitchRadius = [];
    handles.PitchCRAM = [];
    
    guidata(hObject, handles)
end
handles = drawSimple(hObject, handles);
function handles = step1_done_Callback(hObject, eventdata, handles)
if isempty(handles.CRAMlocs)
    errordlg('No cams placed - define cam locations before proceeding');
    return
elseif size(handles.CRAMlocs,1) < 2
    errordlg('You must place at least two cams');
    return
end
handles.step = 2;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function step1_clearAll_Callback(hObject, eventdata, handles)
Button = questdlg('Are you sure you would like to clear all cam locations, connections, and grounds?', 'Clear all?', 'Confirm', 'Cancel', 'Cancel');
switch Button
    case 'Confirm'
        handles.CRAMlocs = zeros(0,2);
        handles.LocMat = zeros(0,4);
        handles.groundvec = [];
        handles.Grounds = [];
        handles.Connections = [];
        handles.PitchRadius = [];
        handles.PitchCRAM = [];
        
        handles.nf = [];
        handles.fmag = [];
        handles.dirMtx = [];
        
        handles.nm = [];
        handles.mmag = [];
        hold off;
        handles = drawSimple(hObject, handles);
end
guidata(hObject, handles)

%% Step 2: Specify ground locations
function handles = step2_newGround_Callback(hObject, eventdata, handles)
if ~isempty(handles.Connections)
    Button = questdlg('Adding a ground cam will delete existing connections','Connection warning','Proceed','Cancel','Proceed');
    switch Button
        case 'Proceed'
            handles.Connections = [];
        case 'Cancel'
            return
    end
end

handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x, y] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x, y); if ~bool; return; end
D = DistMatrix([x y], handles.CRAMlocs);
[~,I] = min(D);
if D(I) < SelectionRadius(handles)
    handles.Grounds = [handles.Grounds, I];
    handles.Grounds = unique(handles.Grounds);
    handles = drawSimple(hObject, handles);
    guidata(hObject, handles)
end
function handles = step2_removeGround_Callback(hObject, eventdata, handles)
if ~isempty(handles.Connections)
    Button = questdlg('Removing a ground will delete existing connections','Connection warning','Proceed','Cancel','Proceed');
    switch Button
        case 'Proceed'
            handles.Connections = [];
        case 'Cancel'
            return
    end
end
handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x, y] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x, y); if ~bool; return; end
D = DistMatrix([x y], handles.CRAMlocs);
[~,I] = min(D);
if D(I) < SelectionRadius(handles)
    handles.Grounds(handles.Grounds == I) = [];
    handles.Connections = [];
    
    guidata(hObject, handles)
end
handles = drawSimple(hObject, handles);
function handles = step2_done_Callback(hObject, eventdata, handles)
if isempty(handles.Grounds)
    errordlg('No ground specified');
    return
elseif length(handles.Grounds) == size(handles.CRAMlocs,1)
    errordlg('Not all cams can be grounded');
    return
end
handles.step = 3;
handles = UpdateDisplay(hObject, handles, handles.step);
handles = ChangeOrder(handles);
guidata(hObject, handles)
function handles = step2_previous_Callback(hObject, eventdata, handles)
handles.step = 1;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function handles = ChangeOrder(handles)
[rg, cg]=size(handles.Grounds);
tempg=zeros(cg,2);
for i=1:cg
    tempg(i,:)=handles.CRAMlocs(handles.Grounds(1,i),:);
end
[rm, cm]=size(handles.CRAMlocs);
tempm=zeros((rm-cg),cm);
count=0;
for i=1:rm
    if(i~=handles.Grounds)
        count=count+1;
        tempm(count,:)=handles.CRAMlocs(i,:);
    end
end
handles.LocMat=zeros(rm,cm+2);
handles.LocMat(1:(rm-cg),2:3)=tempm;
handles.LocMat((rm-cg+1):rm,2:3)=tempg;

handles.groundvec = (rm-cg+1):rm;

%% Step 3: Specify connections
function handles = step3_newConnection_Callback(hObject, eventdata, handles)
handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x1, y1] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x1, y1); if ~bool; return; end
D = DistMatrix([x1, y1], handles.LocMat(:,2:3));
[~,I1] = min(D);
if D(I1) >= SelectionRadius(handles)
    return;
end

handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x2, y2] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x2, y2); if ~bool; return; end
D = DistMatrix([x2, y2], handles.LocMat(:,2:3));
[~,I2] = min(D);
if D(I2) >= SelectionRadius(handles)
    return;
end

if I1 == I2
    return;
end

Is = [I1, I2];
handles.Connections = [handles.Connections; min(Is), max(Is)];
handles.Connections = unique(handles.Connections, 'rows');
handles = drawSimple(hObject, handles);
UpdateAxes(handles);
guidata(hObject, handles)
function handles = step3_autoConnection_Callback(hObject, eventdata, handles)
default = size(handles.LocMat,1)-1;
Answer = inputdlg({'Connections to auto-assign'},'Auto-assign connections',1, {num2str(default)});
if isempty(Answer); return; end
n = str2double(Answer{1});
[~,~,D] = DistMatrix(handles.LocMat(:,2:3), handles.LocMat(:,2:3));
d = sort(D(:));
for k = 1:n
    [r, c] = find(D == d(k));
    D(r, c) = NaN;
    i1(k,1) = min([r, c]);
    i2(k,1) = max([r, c]);
end

Is = [i1, i2];
handles.Connections = [handles.Connections; Is];
handles.Connections = unique(handles.Connections, 'rows');
handles = drawSimple(hObject, handles);
UpdateAxes(handles);
guidata(hObject, handles)
function handles = step3_removeConnection_Callback(hObject, eventdata, handles)
handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x1, y1] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x1, y1); if ~bool; return; end
D = DistMatrix([x1, y1], handles.LocMat(:,2:3));
[~,I1] = min(D);
if D(I1) >= SelectionRadius(handles)
    handles = drawSimple(hObject, handles);
    return;
end

handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x2, y2] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x2, y2); if ~bool; return; end
D = DistMatrix([x2, y2], handles.LocMat(:,2:3));
[~,I2] = min(D);
if D(I2) >= SelectionRadius(handles)
    handles = drawSimple(hObject, handles);
    return;
end

perm1 = bsxfun(@eq, handles.Connections, [I1, I2]); bool1 = sum(perm1, 2) == 2;
perm2 = bsxfun(@eq, handles.Connections, [I2, I1]); bool2 = sum(perm2, 2) == 2;
bool = or(bool1, bool2);
if any(bool)
    handles.Connections(bool,:) = [];
    handles.Connections = unique(handles.Connections, 'rows');
    handles = drawSimple(hObject, handles);
    guidata(hObject, handles)
end
function handles = step3_previous_Callback(hObject, eventdata, handles)
handles.step = 2;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function handles = step3_done_Callback(hObject, eventdata, handles)
if isempty(handles.Connections)
    errordlg('No connections placed - define connections before proceeding')
    return
end
handles.step = 4;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)

%% Step 4
function step4_geometrySettings_Callback(hObject, eventdata, handles)
parameters = handles.img_parameter;
axes(handles.axes)

try
    xtv = handles.axes.XAxis.TickValues;
    ytv = handles.axes.YAxis.TickValues;
catch
    xtv = [];
    ytv = [];
end

hold off
image(parameters)
if ~isempty(xtv)
    handles.axes.XAxis.TickValues = [];
    handles.axes.YAxis.TickValues = [];
end

prompts = {'Strap thickness, t [m]', ...
    'Smallest fabricatable feature size, \delta [m]', ...
    'Amount the initial strap lengths must stretch when they are deformed and assembled within the CRAM, \Delta [m]', ...
    'Layer thickness, W [m]', ...
    'Number of layers, O (must be a multiple of 4 since CRAMS consist of sets of four layers of alternating straps)', ...
    'Number of strap thickness lengths over which the straps are connected to their cams, C', ...
    'Diameter of the bolt or pin holes used to hold the layers together, B [m]'};

title = 'Specify geometric variables';

options.Interpreter = 'tex';
if isempty(handles.defaults) || ~handles.defaultsSet
    %     Rp = min(handles.LocMat(:,4));
    Rp = handles.PitchRadius;
    defaults{1} = num2str(Rp/20); % t
    defaults{2} = num2str(Rp/20); % fabsize, delta
    defaults{3} = num2str(Rp/20); % del, Delta
    defaults{4} = num2str(Rp/2); % W
    defaults{5} = num2str(4); % O
    defaults{6} = num2str(2); % C
    defaults{7} = num2str(Rp/10); % B
    handles.defaults = defaults;
else
    defaults = {num2str(handles.t), num2str(handles.fabsize), num2str(handles.del), num2str(handles.W), num2str(handles.O), num2str(handles.C), num2str(handles.boltr*2)};
end

Answer = inputdlg(prompts, title, 1, defaults, options);

if ~isempty(xtv)
    handles.axes.XAxis.TickValues = xtv;
    handles.axes.YAxis.TickValues = ytv;
end

if isempty(Answer);
    handles = drawSimple(hObject, handles);
    hold on
    handles = DrawCircles(handles);
    UpdateAxes(handles);
    return;
end
handles.t = str2double(Answer{1});
handles.fabsize = str2double(Answer{2}); % will later become handles.t
handles.del = str2double(Answer{3});
handles.W = str2double(Answer{4});
handles.O = str2double(Answer{5});
handles.C = str2double(Answer{6});
handles.boltr = str2double(Answer{7})/2;

handles.defaultsSet = true;

if mod(handles.O, 4) ~= 0
    handles.O = 4*ceil(handles.O/4);
    warning(['Number of layers is not divisible by 4! Rounding to ' num2str(handles.O) ' layers'])
end

if ~isempty(handles.PitchRadius)
    handles = drawSimple(hObject, handles);
    hold off
    handles = QuickRadii(hObject, handles);
    
    UpdateAxes(handles);
else
    handles = drawSimple(hObject, handles);
    UpdateAxes(handles);
end

guidata(hObject, handles)
function step4_materialSettings_Callback(hObject, eventdata, handles)
defaults = {num2str(handles.Et), num2str(handles.Ec), num2str(handles.sigmat),num2str(handles.sigmac), ...
    num2str(handles.dens), num2str(handles.v), num2str(handles.fric)};
prompts = {'Tensile Young''s modulus [N/m^2]', ...
    'Compressive Young''s modulus [N/m^2]', ...
    'Tensile yield stress [N/m^2]', ...
    'Compressive yield stress [N/m^2]', ...
    'Density [kg/m^3]', ...
    'Poisson''s ratio [m/m]', ...
    'Static coefficient of friction [N/N]'};

Answer = inputdlg(prompts, 'Specify material properties', 1, defaults);

if isempty(Answer); return; end
handles.Et = str2double(Answer{1});
handles.Ec = str2double(Answer{2});
handles.sigmat = str2double(Answer{3});
handles.sigmac = str2double(Answer{4});
handles.dens = str2double(Answer{5});
handles.v = str2double(Answer{6});
handles.fric = str2double(Answer{7});

if ~isempty(handles.PitchRadius)
    handles = drawSimple(hObject, handles);
    hold off
    handles = QuickRadii(hObject, handles);
    
    UpdateAxes(handles);
else
    handles = drawSimple(hObject, handles);
    UpdateAxes(handles);
end

guidata(hObject, handles)
function handles = step4_previous_Callback(hObject, eventdata, handles)
handles = drawSimple(hObject, handles);
handles = DrawCircles(handles);
handles.step = 4;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function handles = step4_done_Callback(hObject, eventdata, handles)
if isempty(handles.defaults) || ~handles.defaultsSet
    errordlg('Specify geometric parameters before proceeding')
    return
end
handles = CalculateRadii(hObject, handles);

try
    handles = drawPatchCRAM(handles);
    try drawForces(handles); end
    drawGrounds(handles);
    UpdateAxes(handles);
end

handles.step = 6;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function handles = CalculateRadii(hObject, handles)

LocMat = handles.LocMat;
radius = handles.PitchRadius;
cam = handles.PitchCRAM;

InputConnect = handles.Connections;
[rm, ~]=size(handles.LocMat);
[~, cg]=size(handles.groundvec);

[jointnum, cj]= size(InputConnect);
JointMat=zeros(jointnum,cj+1);
JointMat(1:jointnum,1:2)=InputConnect;
for i=1:jointnum
    pos1=[LocMat(InputConnect(i,1), 2), LocMat(InputConnect(i,1), 3)];
    pos2=[LocMat(InputConnect(i,2), 2), LocMat(InputConnect(i,2), 3)];
    JointMat(i,3)=sqrt(dot((pos1-pos2),(pos1-pos2))); % this puts the distance between the centers of each joint in the third column of JointMat
end

LocMat(cam,4)=radius;
count=0;
camhistory=zeros(1,rm);
camhistory(1,1)=cam;
repeat=prod(camhistory);
while (repeat==0)
    for i=1:jointnum
        if(JointMat(i,1)==cam || JointMat(i,2)==cam)
            radius=JointMat(i,3)-LocMat(cam,4);
            if(JointMat(i,1)==cam)
                LocMat(JointMat(i,2),4)=radius;
                if(JointMat(i,2)~=camhistory)
                    count=count+1;
                    cam=JointMat(i,2);
                    camhistory(1,count+1)=cam;
                else
                    for j=1:rm
                        if(j~=camhistory)
                            for k=1:jointnum
                                if(JointMat(k,1)==j && ismember(JointMat(k,2),camhistory)==1)
                                    radius=JointMat(k,3)-LocMat(JointMat(k,2),4);
                                    LocMat(JointMat(k,1),4)=radius;
                                    cam=JointMat(k,1);
                                    count=count+1;
                                    camhistory(1,count+1)=cam;
                                elseif(JointMat(k,2)==j && ismember(JointMat(k,1),camhistory)==1)
                                    radius=JointMat(k,3)-LocMat(JointMat(k,1),4);
                                    LocMat(JointMat(k,2),4)=radius;
                                    cam=JointMat(k,2);
                                    count=count+1;
                                    camhistory(1,count+1)=cam;
                                end
                            end
                        end
                    end
                end
            else
                LocMat(JointMat(i,1),4)=radius;
                if(JointMat(i,1)~=camhistory)
                    count=count+1;
                    cam=JointMat(i,1);
                    camhistory(1,count+1)=cam;
                else
                    for j=1:rm
                        if(j~=camhistory)
                            for k=1:jointnum
                                if(JointMat(k,1)==j && ismember(JointMat(k,2),camhistory)==1)
                                    radius=JointMat(k,3)-LocMat(JointMat(k,2),4);
                                    LocMat(JointMat(k,1),4)=radius;
                                    cam=JointMat(k,1);
                                    count=count+1;
                                    camhistory(1,count+1)=cam;
                                elseif(JointMat(k,2)==j && ismember(JointMat(k,1),camhistory)==1)
                                    radius=JointMat(k,3)-LocMat(JointMat(k,1),4);
                                    LocMat(JointMat(k,2),4)=radius;
                                    cam=JointMat(k,2);
                                    count=count+1;
                                    camhistory(1,count+1)=cam;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    repeat=prod(camhistory);
end
bad=0;
for i=1:rm
    if(LocMat(i,4)<=0)
        bad=1;
    end
end
[YoN, culpritvec]=overlap(LocMat);
if(bad==1)
    h = warndlg('The specified design is not geometrically feasible. Try entering a different pitch radius or rearrange the locations of the circular cams.');
    set(h, 'WindowStyle', 'modal')
elseif(YoN==1)
    h = warndlg('Warning: Some of the circular cams chosen are overlapping. Please correct the design by reducing their radii accordingly.');
    set(h, 'WindowStyle', 'modal')
    %     hold on
    %     for i=1:rm
    %         plaincircle(LocMat(i,2),LocMat(i,3),LocMat(i,4))
    %     end
else
    %%%%%% the code now draws the geometrically compatible pitch circles
    %     for i=1:rm
    %         plaincircle(LocMat(i,2),LocMat(i,3),LocMat(i,4))
    %     end
    
    %%%%%% the code now determines the optimal strap attachment angle for all the cams
    %%%%%%
    % Add twice as many rows to FullJointMat by s swaping the order and direction of the joint
    FullJointMat=zeros(2*jointnum,4);
    FullJointMat(1:jointnum,1:2)=JointMat(:,1:2);
    FullJointMat((jointnum+1):(2*jointnum),1)=JointMat(:,2);
    FullJointMat((jointnum+1):(2*jointnum),2)=JointMat(:,1);
    % First make the strap angle halfway between neighboring joint angles
    % unless the strap angle would be more than pi/2 then leave it that way
    for i=1:(2*jointnum)
        jointvec=LocMat(FullJointMat(i,2),2:3)-LocMat(FullJointMat(i,1),2:3);
        unjovec=jointvec/sqrt(dot(jointvec,jointvec));
        if(unjovec(1,2)>=0)
            anjointvec=acos(dot([1, 0], unjovec));
        else
            anjointvec=(2*pi)-acos(dot([1, 0], unjovec));
        end
        FullJointMat(i,3)=anjointvec;
    end
    temp=zeros(1,2);
    for i=1:rm
        count=1;
        for j=1:(2*jointnum)
            if(FullJointMat(j,1)==i)
                temp(count,1)=j;
                temp(count,2)=FullJointMat(j,3);
                count=count+1;
            end
        end
        count=count-1;
        orderangle=sort(temp(:,2));
        for k=1:count
            for m=1:count
                if(temp(m,2)==orderangle(k,1))
                    orderangle(k,2)=temp(m,1);
                end
            end
        end
        temp=zeros(1,2);
        if(count==1)
            FullJointMat(orderangle(1,2),4)=(pi/2);
        else
            for n=1:count-1
                if(((orderangle(n+1,1)-orderangle(n,1))/2) > (pi/2))
                    FullJointMat(orderangle(n,2),4)=(pi/2);
                else
                    FullJointMat(orderangle(n,2),4)=(orderangle(n+1,1)-orderangle(n,1))/2;
                end
            end
            if(((orderangle(1,1)-(orderangle(count,1)-(2*pi)))/2) > (pi/2))
                FullJointMat(orderangle(count,2),4)=(pi/2);
            else
                FullJointMat(orderangle(count,2),4)=(orderangle(1,1)-(orderangle(count,1)-(2*pi)))/2;
            end
        end
        orderangle=zeros(1,2);
    end
    % Now shorten the longer straps so that everything is geometrically compatible on the front layer
    temp=zeros(1,2);
    for i=1:rm
        count=1;
        for j=1:(2*jointnum)
            if(FullJointMat(j,1)==i)
                temp(count,1)=j;
                temp(count,2)=FullJointMat(j,3);
                count=count+1;
            end
        end
        count=count-1;
        orderangle=sort(temp(:,2));
        for k=1:count
            for m=1:count
                if(temp(m,2)==orderangle(k,1))
                    orderangle(k,2)=temp(m,1);
                end
            end
        end
        temp=zeros(1,2);
        if(count>1)
            for n=1:count-1
                dist1=FullJointMat(orderangle(n,2),4)*LocMat(FullJointMat(orderangle(n,2),1),4);
                for o=1:(2*jointnum)
                    if(FullJointMat(orderangle(n+1,2),1)==FullJointMat(o,2) && FullJointMat(orderangle(n+1,2),2)==FullJointMat(o,1))
                        dist2=FullJointMat(o,4)*LocMat(FullJointMat(o,1),4);
                        numo=o;
                    end
                end
                disttot=(FullJointMat(orderangle(n+1,2),3)-FullJointMat(orderangle(n,2),3))*LocMat(FullJointMat(orderangle(n,2),1),4);
                if((dist1+dist2)>disttot)
                    if(dist1>=dist2)
                        FullJointMat(orderangle(n,2),4)=FullJointMat(orderangle(n,2),4)+((disttot-(dist1+dist2))/LocMat(FullJointMat(orderangle(n,2),1),4));
                    else
                        FullJointMat(numo,4)=FullJointMat(numo,4)+((disttot-(dist1+dist2))/LocMat(FullJointMat(numo,1),4));
                    end
                end
            end
            dist1=FullJointMat(orderangle(count,2),4)*LocMat(FullJointMat(orderangle(count,2),1),4);
            for o=1:(2*jointnum)
                if(FullJointMat(orderangle(1,2),1)==FullJointMat(o,2) && FullJointMat(orderangle(1,2),2)==FullJointMat(o,1))
                    dist2=FullJointMat(o,4)*LocMat(FullJointMat(o,1),4);
                    numo=o;
                end
            end
            disttot=((2*pi)+FullJointMat(orderangle(1,2),3)-FullJointMat(orderangle(count,2),3))*LocMat(FullJointMat(orderangle(count,2),1),4);
            if((dist1+dist2)>disttot)
                if(dist1>=dist2)
                    FullJointMat(orderangle(count,2),4)=FullJointMat(orderangle(count,2),4)+((disttot-(dist1+dist2))/LocMat(FullJointMat(orderangle(count,2),1),4));
                else
                    FullJointMat(numo,4)=FullJointMat(numo,4)+((disttot-(dist1+dist2))/LocMat(FullJointMat(numo,1),4));
                end
            end
        end
    end
    % Now shorten more of the straps so that the back layer is also geometrically compatible
    for pp=1:jointnum
        dist1a=FullJointMat(pp,4)*LocMat(FullJointMat(pp,1),4);
        for kk=jointnum+1:(2*jointnum)
            if(FullJointMat(pp,1)==FullJointMat(kk,2) && FullJointMat(pp,2)==FullJointMat(kk,1))
                dist2a=FullJointMat(kk,4)*LocMat(FullJointMat(kk,1),4);
                numkk=kk;
            end
        end
        if(cutdec(dist2a,10)>cutdec(dist1a,10))
            FullJointMat(numkk,4)=dist1a/(LocMat(FullJointMat(numkk,1),4));
        elseif(cutdec(dist2a,10)<cutdec(dist1a,10))
            FullJointMat(pp,4)=dist2a/(LocMat(FullJointMat(pp,1),4));
        end
    end
    
    Et = handles.Et;
    Ec = handles.Ec;
    sigmat = handles.sigmat;
    sigmac = handles.sigmac;
    dens = handles.dens;
    v = handles.v;
    fric = handles.fric;
    fabsize = handles.fabsize;
    del = handles.del;
    W = handles.W;
    O = handles.O;
    C = handles.C;
    boltr = handles.boltr;
    
    Rmin=min(LocMat(:,4));
    if(sigmat/Et<sigmac/Ec)
        tmax=2*Rmin*(sigmat/Et);
    else
        tmax=2*Rmin*(sigmac/Et);
    end
    t=handles.t; % We pick this as an example
    
    dcut1=fabsize;
    dcut2=fabsize;
    
    for j=1:(2*jointnum)
        Rp1=LocMat(FullJointMat(j,1),4);
        Rp2=LocMat(FullJointMat(j,2),4);
        Rb1=Rp1-(t/2);
        Rb2=Rp2-(t/2);
        FullJointMat(j,5)=(pi/2)-(FullJointMat(j,4)-(C*t/Rb1));
        for k=1:(2*jointnum)
            if(FullJointMat(j,1)==FullJointMat(k,2) && FullJointMat(j,2)==FullJointMat(k,1))
                FullJointMat(j,6)=(pi/2)-(FullJointMat(k,4)-(C*t/Rb2));
            end
        end
        FullJointMat(j,7)=Rp1;
        FullJointMat(j,8)=Rp2;
    end
    FullJointMat(:,4:7)=FullJointMat(:,5:8);
    FullJointMattemp=FullJointMat(:,1:7);
    FullJointMat=FullJointMattemp;
    
    %%% now the code determines if anything will have yielded when the lattice
    %%% is assembled. In the process, it will calculate the assembled strap
    %%% lengths of each joint and add those to an eighth column within
    %%% FullJointMat
    
    for i=1:(2*jointnum)
        FullJointMat(i,8)=(((pi/2)-FullJointMat(i,4))*FullJointMat(i,6))+(((pi/2)-FullJointMat(i,5))*FullJointMat(i,7));
    end
    yield=0;
    for i=1:(2*jointnum)
        stress1=((Et*del)/(FullJointMat(i,8)-del))+((Et*t)/(2*FullJointMat(i,6)));
        stress2=((Et*del)/(FullJointMat(i,8)-del))-((Ec*t)/(2*FullJointMat(i,6)));
        if(stress1>=sigmat || stress2>=sigmat || stress2>=sigmac)
            yield=1;
        end
    end
    
    [dofn]=DOF(FullJointMat,LocMat,jointnum,rm,cg);
    
    if(yield==1)
        h = warndlg({'The current design will yield when it is assembled. We recommend that you use a different material, make the straps thinner, and/or make the smallest cam radius larger.'; ...
            ' ';
            ['Note that the CRAM designed possesses ' num2str(dofn) ' degrees of freedom.']}, 'Yield warning and degree of freedom count');
    else
        h = msgbox({'Congratulations. Your current design will not yield when it is assembled.'; ...
            ' ';
            ['Note that the CRAM designed possesses ' num2str(dofn) ' degrees of freedom.']}, 'Degree of freedom count');
    end
    set(h, 'WindowStyle', 'modal')
    
    handles.dcut1 = dcut1;
    handles.dcut2 = dcut2;
    handles.t = t;
    handles.FullJointMat = FullJointMat;
    handles.LocMat = LocMat;
    handles.jointnum = jointnum;
    
    figure(handles.figure1)
    axes(handles.axes)
    
    guidata(hObject, handles)
end

%% Step 5
function step5_pitchRadius_Callback(hObject, eventdata, handles)
handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x, y] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x, y); if ~bool; return; end
D = DistMatrix([x y], handles.LocMat(:,2:3));
[~,I] = min(D);
if D(I) < SelectionRadius(handles)
    handles.PitchCRAM = I;
    guidata(hObject, handles)
else
    return
end

if handles.LocMat(I,4) > 0
    default_value = handles.LocMat(I,4);
else
    [~, Ddiag] = DistMatrix(handles.LocMat(:,2:3), handles.LocMat(:,2:3));
    default_value = min(Ddiag(:)/2);
end
Answer = inputdlg({'Pitch radius of selected cam'},'Set cam pitch radius',1,{num2str(default_value)});
if isempty(Answer); return; end
handles.PitchRadius = str2double(Answer{1});

handles = QuickRadii(hObject, handles);
% handles = drawSimple(hObject, handles);

guidata(hObject, handles)
function handles = step5_previous_Callback(hObject, eventdata, handles)
if ~isempty(handles.nf) || ~isempty(handles.nm) || ~isempty(handles.PitchRadius)
    Button = questdlg('Returning to previous step will clear pitch radii and all loads', 'Warning', 'Proceed', 'Cancel', 'Proceed');
    switch Button
        case 'Proceed'
            handles.nf = [];
            handles.fmag = [];
            handles.dirMtx = [];
            
            handles.nm = [];
            handles.mmag = [];
            
            handles.PitchRadius = [];
            handles.LocMat(:,4) = 0;
            handles.PitchCRAM = [];
        case 'Cancel'
            return
    end
end
handles = drawSimple(hObject, handles);

handles.step = 3;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function handles = step5_done_Callback(hObject, eventdata, handles)
if isempty(handles.PitchRadius) || isempty(handles.PitchCRAM)
    errordlg('No pitch specified!')
    return
end

handles.step = 5;

handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function handles = QuickRadii(hObject, handles)
LocMat = handles.LocMat;
radius = handles.PitchRadius;
cam = handles.PitchCRAM;

InputConnect = handles.Connections;
[rm, ~]=size(handles.LocMat);
[~, cg]=size(handles.groundvec);

[jointnum, cj]= size(InputConnect);
JointMat=zeros(jointnum,cj+1);
JointMat(1:jointnum,1:2)=InputConnect;
for i=1:jointnum
    pos1=[LocMat(InputConnect(i,1), 2), LocMat(InputConnect(i,1), 3)];
    pos2=[LocMat(InputConnect(i,2), 2), LocMat(InputConnect(i,2), 3)];
    JointMat(i,3)=sqrt(dot((pos1-pos2),(pos1-pos2))); % this puts the distance between the centers of each joint in the third column of JointMat
end

LocMat(cam,4)=radius;
count=0;
camhistory=zeros(1,rm);
camhistory(1,1)=cam;
repeat=prod(camhistory);
while (repeat==0)
    for i=1:jointnum
        if(JointMat(i,1)==cam || JointMat(i,2)==cam)
            radius=JointMat(i,3)-LocMat(cam,4);
            if(JointMat(i,1)==cam)
                LocMat(JointMat(i,2),4)=radius;
                if(JointMat(i,2)~=camhistory)
                    count=count+1;
                    cam=JointMat(i,2);
                    camhistory(1,count+1)=cam;
                else
                    for j=1:rm
                        if(j~=camhistory)
                            for k=1:jointnum
                                if(JointMat(k,1)==j && ismember(JointMat(k,2),camhistory)==1)
                                    radius=JointMat(k,3)-LocMat(JointMat(k,2),4);
                                    LocMat(JointMat(k,1),4)=radius;
                                    cam=JointMat(k,1);
                                    count=count+1;
                                    camhistory(1,count+1)=cam;
                                elseif(JointMat(k,2)==j && ismember(JointMat(k,1),camhistory)==1)
                                    radius=JointMat(k,3)-LocMat(JointMat(k,1),4);
                                    LocMat(JointMat(k,2),4)=radius;
                                    cam=JointMat(k,2);
                                    count=count+1;
                                    camhistory(1,count+1)=cam;
                                end
                            end
                        end
                    end
                end
            else
                LocMat(JointMat(i,1),4)=radius;
                if(JointMat(i,1)~=camhistory)
                    count=count+1;
                    cam=JointMat(i,1);
                    camhistory(1,count+1)=cam;
                else
                    for j=1:rm
                        if(j~=camhistory)
                            for k=1:jointnum
                                if(JointMat(k,1)==j && ismember(JointMat(k,2),camhistory)==1)
                                    radius=JointMat(k,3)-LocMat(JointMat(k,2),4);
                                    LocMat(JointMat(k,1),4)=radius;
                                    cam=JointMat(k,1);
                                    count=count+1;
                                    camhistory(1,count+1)=cam;
                                elseif(JointMat(k,2)==j && ismember(JointMat(k,1),camhistory)==1)
                                    radius=JointMat(k,3)-LocMat(JointMat(k,1),4);
                                    LocMat(JointMat(k,2),4)=radius;
                                    cam=JointMat(k,2);
                                    count=count+1;
                                    camhistory(1,count+1)=cam;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    repeat=prod(camhistory);
end
bad=0;
for i=1:rm
    if(LocMat(i,4)<=0)
        bad=1;
    end
end
[YoN, culpritvec]=overlap(LocMat);
if(bad==1)
    h = warndlg('The specified design is not geometrically feasible. Try entering a different pitch radius or rearrange the locations of the circular cams.');
    set(h, 'WindowStyle', 'modal')
elseif(YoN==1)
    h = warndlg('Warning: Some of the circular cams chosen are overlapping. Please correct the design by reducing their radii accordingly.');
    set(h, 'WindowStyle', 'modal')
    figure(handles.figure1)
    axes(handles.axes)
    %     hold on
    %     for i=1:rm
    %         plaincircle(LocMat(i,2),LocMat(i,3),LocMat(i,4))
    %     end
    axis equal
else
    %%%%%% the code now draws the geometrically compatible pitch circles
    %     for i=1:rm
    %         plaincircle(LocMat(i,2),LocMat(i,3),LocMat(i,4))
    %     end
    axis equal
    
    %%%%%% the code now determines the optimal strap attachment angle for all the cams
    %%%%%%
    % Add twice as many rows to FullJointMat by s swaping the order and direction of the joint
    FullJointMat=zeros(2*jointnum,4);
    FullJointMat(1:jointnum,1:2)=JointMat(:,1:2);
    FullJointMat((jointnum+1):(2*jointnum),1)=JointMat(:,2);
    FullJointMat((jointnum+1):(2*jointnum),2)=JointMat(:,1);
    % First make the strap angle halfway between neighboring joint angles
    % unless the strap angle would be more than pi/2 then leave it that way
    for i=1:(2*jointnum)
        jointvec=LocMat(FullJointMat(i,2),2:3)-LocMat(FullJointMat(i,1),2:3);
        unjovec=jointvec/sqrt(dot(jointvec,jointvec));
        if(unjovec(1,2)>=0)
            anjointvec=acos(dot([1, 0], unjovec));
        else
            anjointvec=(2*pi)-acos(dot([1, 0], unjovec));
        end
        FullJointMat(i,3)=anjointvec;
    end
    temp=zeros(1,2);
    for i=1:rm
        count=1;
        for j=1:(2*jointnum)
            if(FullJointMat(j,1)==i)
                temp(count,1)=j;
                temp(count,2)=FullJointMat(j,3);
                count=count+1;
            end
        end
        count=count-1;
        orderangle=sort(temp(:,2));
        for k=1:count
            for m=1:count
                if(temp(m,2)==orderangle(k,1))
                    orderangle(k,2)=temp(m,1);
                end
            end
        end
        temp=zeros(1,2);
        if(count==1)
            FullJointMat(orderangle(1,2),4)=(pi/2);
        else
            for n=1:count-1
                if(((orderangle(n+1,1)-orderangle(n,1))/2) > (pi/2))
                    FullJointMat(orderangle(n,2),4)=(pi/2);
                else
                    FullJointMat(orderangle(n,2),4)=(orderangle(n+1,1)-orderangle(n,1))/2;
                end
            end
            if(((orderangle(1,1)-(orderangle(count,1)-(2*pi)))/2) > (pi/2))
                FullJointMat(orderangle(count,2),4)=(pi/2);
            else
                FullJointMat(orderangle(count,2),4)=(orderangle(1,1)-(orderangle(count,1)-(2*pi)))/2;
            end
        end
        orderangle=zeros(1,2);
    end
    % Now shorten the longer straps so that everything is geometrically compatible on the front layer
    temp=zeros(1,2);
    for i=1:rm
        count=1;
        for j=1:(2*jointnum)
            if(FullJointMat(j,1)==i)
                temp(count,1)=j;
                temp(count,2)=FullJointMat(j,3);
                count=count+1;
            end
        end
        count=count-1;
        orderangle=sort(temp(:,2));
        for k=1:count
            for m=1:count
                if(temp(m,2)==orderangle(k,1))
                    orderangle(k,2)=temp(m,1);
                end
            end
        end
        temp=zeros(1,2);
        if(count>1)
            for n=1:count-1
                dist1=FullJointMat(orderangle(n,2),4)*LocMat(FullJointMat(orderangle(n,2),1),4);
                for o=1:(2*jointnum)
                    if(FullJointMat(orderangle(n+1,2),1)==FullJointMat(o,2) && FullJointMat(orderangle(n+1,2),2)==FullJointMat(o,1))
                        dist2=FullJointMat(o,4)*LocMat(FullJointMat(o,1),4);
                        numo=o;
                    end
                end
                disttot=(FullJointMat(orderangle(n+1,2),3)-FullJointMat(orderangle(n,2),3))*LocMat(FullJointMat(orderangle(n,2),1),4);
                if((dist1+dist2)>disttot)
                    if(dist1>=dist2)
                        FullJointMat(orderangle(n,2),4)=FullJointMat(orderangle(n,2),4)+((disttot-(dist1+dist2))/LocMat(FullJointMat(orderangle(n,2),1),4));
                    else
                        FullJointMat(numo,4)=FullJointMat(numo,4)+((disttot-(dist1+dist2))/LocMat(FullJointMat(numo,1),4));
                    end
                end
            end
            dist1=FullJointMat(orderangle(count,2),4)*LocMat(FullJointMat(orderangle(count,2),1),4);
            for o=1:(2*jointnum)
                if(FullJointMat(orderangle(1,2),1)==FullJointMat(o,2) && FullJointMat(orderangle(1,2),2)==FullJointMat(o,1))
                    dist2=FullJointMat(o,4)*LocMat(FullJointMat(o,1),4);
                    numo=o;
                end
            end
            disttot=((2*pi)+FullJointMat(orderangle(1,2),3)-FullJointMat(orderangle(count,2),3))*LocMat(FullJointMat(orderangle(count,2),1),4);
            if((dist1+dist2)>disttot)
                if(dist1>=dist2)
                    FullJointMat(orderangle(count,2),4)=FullJointMat(orderangle(count,2),4)+((disttot-(dist1+dist2))/LocMat(FullJointMat(orderangle(count,2),1),4));
                else
                    FullJointMat(numo,4)=FullJointMat(numo,4)+((disttot-(dist1+dist2))/LocMat(FullJointMat(numo,1),4));
                end
            end
        end
    end
    % Now shorten more of the straps so that the back layer is also geometrically compatible
    for pp=1:jointnum
        dist1a=FullJointMat(pp,4)*LocMat(FullJointMat(pp,1),4);
        for kk=jointnum+1:(2*jointnum)
            if(FullJointMat(pp,1)==FullJointMat(kk,2) && FullJointMat(pp,2)==FullJointMat(kk,1))
                dist2a=FullJointMat(kk,4)*LocMat(FullJointMat(kk,1),4);
                numkk=kk;
            end
        end
        if(cutdec(dist2a,10)>cutdec(dist1a,10))
            FullJointMat(numkk,4)=dist1a/(LocMat(FullJointMat(numkk,1),4));
        elseif(cutdec(dist2a,10)<cutdec(dist1a,10))
            FullJointMat(pp,4)=dist2a/(LocMat(FullJointMat(pp,1),4));
        end
    end
    
    %     Rmin=min(LocMat(:,4));
    
    handles.LocMat = LocMat;
    
    handles = drawSimple(hObject, handles);
    hold on
    handles = DrawCircles(handles);
    
    guidata(hObject, handles)
end

%% Step 6: Add loads
function step6_addForce_Callback(hObject, eventdata, handles)
handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x, y] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x, y); if ~bool; return; end
D = DistMatrix([x y], handles.LocMat(:,2:3));
[~,I] = min(D);
if any(I == handles.groundvec)
    h = warndlg('Grounded cams are not allowed to be loaded.', 'Ground load warning');
    set(h, 'WindowStyle', 'modal')
    return
end
if D(I) >= SelectionRadius(handles)
    return
end

Answer = inputdlg({'Force magnitude [N]', 'X-component of force', 'Y-component of force'},'Enter force',1,{'10', '1', '1'});
if isempty(Answer); return; end

dirforce = [str2double(Answer{2}), str2double(Answer{3})];
if size(dirforce,1) ~= 1 || size(dirforce,2) ~= 2
    errordlg('Force direction must be a 1x2 vector');
    return
end
if (dot(dirforce,dirforce)==0)
    unitdirforce = [0, 0];
else
    unitdirforce = dirforce/sqrt(dot(dirforce,dirforce)) * sign(str2double(Answer{1}));
end

if any(handles.nf == I)
    handles.fmag(handles.nf == I,1) = abs(str2double(Answer{1}));
    handles.nf(handles.nf == I,1) = I; % Cam to load
    handles.dirMtx(handles.nf == I,1:2) = unitdirforce;
else
    handles.fmag(end+1,1) = str2double(Answer{1});
    handles.nf(end+1,1) = I; % Cam to load
    handles.dirMtx(end+1,1:2) = unitdirforce;
end

handles = drawPatchCRAM(handles);
drawForces(handles);
drawGrounds(handles);
UpdateAxes(handles);

guidata(hObject, handles)
function step6_addMoment_Callback(hObject, eventdata, handles)
handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x, y] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x, y); if ~bool; return; end
D = DistMatrix([x y], handles.LocMat(:,2:3));
[~,I] = min(D);
if any(I == handles.groundvec)
    h = warndlg('Grounded cams are not allowed to be loaded.', 'Ground load warning');
    set(h, 'WindowStyle', 'modal')
    return
end
if D(I) >= SelectionRadius(handles)
    return
end

Answer = inputdlg({'Moment magnitude [Nm] and direction (negative is clockwise)'},'Enter Moment',1,{'1'});
if isempty(Answer); return; end

if any(handles.nm == I)
    handles.mmag(handles.nm == I,1) = str2double(Answer{1});
    handles.nm(handles.nm == I,1) = I; % Cam to load
else
    handles.mmag(end+1,1) = str2double(Answer{1});
    handles.nm(end+1,1) = I; % Cam to load
end

handles = drawPatchCRAM(handles);
drawForces(handles);
drawGrounds(handles);
UpdateAxes(handles);

guidata(hObject, handles)
function handles = step6_remove_Callback(hObject, eventdata, handles)
handles = UpdateDisplay(hObject, handles, handles.step, 'clear');
[x, y] = cinput(1);
handles = UpdateDisplay(hObject, handles, handles.step);

bool = BoundCheck(handles, x, y); if ~bool; return; end
D = DistMatrix([x y], handles.LocMat(:,2:3));
[~,I] = min(D);
if D(I) >= SelectionRadius(handles)
    return
end

fI = handles.nf == I;
handles.fmag(fI,:) = [];
handles.dirMtx(fI,:) = [];
handles.nf(fI,:) = [];

mI = handles.nm == I;
handles.mmag(mI,:) = [];
handles.nm(mI,:) = [];

handles = drawPatchCRAM(handles);
drawForces(handles);
UpdateAxes(handles);

guidata(hObject, handles)
function handles = step6_previous_Callback(hObject, eventdata, handles)
handles = drawSimple(hObject, handles);
handles = DrawCircles(handles);
UpdateAxes(handles);
handles.step = 5;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function handles = step6_done_Callback(hObject, eventdata, handles)
if isempty(handles.nf) && isempty(handles.nm) && get(handles.step6_gravity,'Value') == 0
    errordlg('No loads specified');
    return
end
if isempty(handles.mmag)
    handles.mmag = 0;
end
try
    handles = drawPatchCRAM(handles);
    drawForces(handles);
    drawGrounds(handles);
    UpdateAxes(handles);
end
handles.step = 7;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function step6_gravity_Callback(hObject, eventdata, handles)
value = get(handles.step6_gravity, 'Value');
handles = drawPatchCRAM(handles);
drawForces(handles);
drawGrounds(handles);
UpdateAxes(handles);
guidata(hObject, handles)
function [handles, Wrench, gravforce] = UpdateGravity(handles)
value = get(handles.step6_gravity, 'Value');

[rm, ~]=size(handles.CRAMlocs);
[~, cg]=size(handles.groundvec);

Wrench=zeros(3*(rm-cg),1);
LocMat = handles.LocMat;
W = handles.W;
O = handles.O;
t = handles.t;
dens = handles.dens;

gravforce=zeros(rm-cg,1);

if value == 1
    unitgravdir = [0, -1];
    gravforce=zeros(rm-cg,1);
    gravsumload=zeros(rm-cg,1);
    fincg=zeros(rm-cg,1);
    for i=1:rm-cg
        gravforce(i,1)=pi*((LocMat(i,4)-(t/2))^2)*W*O*dens*9.81;
        fincg(i,1)=gravforce(i,1)/handles.num;
        camforce=fincg(i,1)*unitgravdir;
        camLoc=[LocMat(i,2), LocMat(i,3)];
        Wrench((3*(i-1))+1:(3*i),1)=Wrench((3*(i-1))+1:(3*i),1)+[camforce(1,1); camforce(1,2); ((camLoc(1,1)*camforce(1,2))-(camLoc(1,2)*camforce(1,1)))];
    end
    handles.fincg = fincg;
else
    handles.fincg = zeros(rm-cg,1);
end
handles.gravforce = gravforce;
handles.Wrench = Wrench;
function [handles, Wrench] = CreateWrench(handles)
[handles, Wrench] = UpdateGravity(handles);

if length(handles.nf) > length(handles.nm)
    len = length(handles.nf);
    handles.nm(end+1:len,1) = 1;
    handles.mmag(end+1:len,1) = 0;
elseif length(handles.nf) < length(handles.nm)
    len = length(handles.nm);
    handles.nf(end+1:len,1) = 1;
    handles.fmag(end+1:len,1) = 0;
    handles.dirMtx(len,1:2) = [0 0];
end

nf = handles.nf;
nm = handles.nm;

for cct = 1:length(handles.nf)
    fmag = handles.fmag(cct, 1);
    dirMtx = handles.dirMtx(cct,1:2);
    finc = fmag/handles.num;
    camforce = finc * dirMtx;
    camLoc=[handles.LocMat(nf(cct,1),2), handles.LocMat(nf(cct,1),3)];
    Wrench((3*(nf(cct,1)-1))+1:(3*nf(cct,1)),1) = Wrench((3*(nf(cct,1)-1))+1:(3*nf(cct,1)),1)+[camforce(1,1); camforce(1,2); ((camLoc(1,1)*camforce(1,2))-(camLoc(1,2)*camforce(1,1)))];
end

for cct = 1:length(handles.nm)
    minc = handles.mmag/handles.num;
    Wrench((3*(nm(cct,1)-1))+1:(3*nm(cct,1)),1) = Wrench((3*(nm(cct,1)-1))+1:(3*nm(cct,1)),1)+[0; 0; (minc(cct,1))];
end
handles.Wrench = Wrench;

%% Step 7: Simulate
function step7_simulate_Callback(hObject, eventdata, handles)
set(handles.step7_stop, 'Enable', 'On')
set(handles.step7_stop, 'Value', 0)
set(handles.text_status, 'String', 'Running simulation...')
UpdateAxes(handles);

handles.LastGIF = [];
axes(handles.axes)
set(handles.axes, 'Box', 'off')
set(handles.axes, 'XTick', [])
set(handles.axes, 'YTick', [])
xl = xlim;
yl = ylim;

FullJointMat = handles.FullJointMat;
LocMat = handles.LocMat;
jointnum = handles.jointnum;
C = handles.C;
boltr = handles.boltr;
t = handles.t;
W = handles.W;
del = handles.del;
dcut1 = handles.dcut1;
dcut2 = handles.dcut2;
O = handles.O;
Et = handles.Et;
Ec = handles.Ec;
v = handles.v;
fric = handles.fric;
fabsize = handles.fabsize;
frames = str2double(get(handles.step7_editFrames, 'String'));

rm = size(handles.LocMat,1);
cg = size(handles.groundvec,2);

[handles, Wrench] = CreateWrench(handles);
gravforce = handles.gravforce;
incder=fabsize/handles.derivative_factor; % derivative increment
targetnum=handles.targetnum;

dirMtx = handles.dirMtx;
gravsumload=zeros(rm-cg,1);
finc = handles.fmag/handles.num;
minc = handles.mmag/handles.num;
fincg = handles.fincg;
fyes = sum(handles.fmag) > 0;
myes = sum(handles.mmag) > 0;
cct = length(handles.nf);

fmag = handles.fmag;
mmag = handles.mmag;
fmagsum = zeros(size(handles.fmag));
mmagsum = zeros(size(handles.mmag));
nf = handles.nf;
nm = handles.nm;

Rmax=max(LocMat(:,4));
targetphi=fabsize/(targetnum*Rmax);
targetxy=fabsize/targetnum;

%%% Now the code starts to run the FEA
grav = get(handles.step6_gravity, 'Value');
if(grav==1)
    loadfull=gravforce(1,1);
    load=gravsumload(1,1);
elseif(fyes==1)
    temp=1;
    while(temp<=cct)
        if(abs(fmag(temp,1))>0)
            loadfull=fmag(temp,1);
            load=fmagsum(temp,1);
            fnumber=temp;
            break
        end
        temp=temp+1;
    end
elseif(myes==1)
    temp=1;
    while(temp<=cct)
        if(abs(mmag(temp,1))>0)
            loadfull=mmag(temp,1);
            load=mmagsum(temp,1);
            mnumber=temp;
            break
        end
        temp=temp+1;
    end
else
    loadfull=0;
    load=0;
end

count=1;
while (abs(load)<abs(loadfull))
    nope=0;
    while (nope<1)
        oldLocMat=LocMat;
        [stiffmat]=BuildStiffnessMatrix(FullJointMat, LocMat, jointnum, rm, cg, t, W, del, dcut1, dcut2, O, Et, Ec, v, fric, incder);
        Twistglobal=pinv(stiffmat)*Wrench;
        for m=1:(rm-cg)
            transM=[1, 0, 0; LocMat(m,3), 1, 0; -LocMat(m,2), 0, 1];
            Tinc=real(transM\(Twistglobal((3*(m-1))+1:(3*m),1)));
            LocMat(m,1)=LocMat(m,1)+Tinc(1,1);
            LocMat(m,2)=LocMat(m,2)+Tinc(2,1);
            LocMat(m,3)=LocMat(m,3)+Tinc(3,1);
        end
        DifferenceMat=abs(LocMat-oldLocMat);
        phiangles=DifferenceMat(:,1);
        xytrans=DifferenceMat(:,2:3);
        maxphidiff=max(phiangles);
        maxxydiff=max(max(xytrans));
        if(maxphidiff==0 && maxxydiff==0)
            if(grav==1)
                temp=1;
                while(temp<=rm-cg)
                    fincg(temp,1)=fincg(temp,1)*10;
                    temp=temp+1;
                end
            end
            if(fyes==1)
                temp=1;
                while(temp<=cct)
                    finc(temp,1)=finc(temp,1)*10;
                    temp=temp+1;
                end
            end
            if(myes==1)
                temp=1;
                while(temp<=cct)
                    minc(temp,1)=minc(temp,1)*10;
                    temp=temp+1;
                end
            end
            LocMat=oldLocMat;
        elseif((maxxydiff<(targetxy/2)) && (maxphidiff<(targetphi/2)))
            if(maxphidiff>maxxydiff)
                if(grav==1)
                    temp=1;
                    while(temp<=rm-cg)
                        fincg(temp,1)=fincg(temp,1)*(targetphi/maxphidiff);
                        temp=temp+1;
                    end
                end
                if(fyes==1)
                    temp=1;
                    while(temp<=cct)
                        finc(temp,1)=finc(temp,1)*(targetphi/maxphidiff);
                        temp=temp+1;
                    end
                end
                if(myes==1)
                    temp=1;
                    while(temp<=cct)
                        minc(temp,1)=minc(temp,1)*(targetphi/maxphidiff);
                        temp=temp+1;
                    end
                end
                LocMat=oldLocMat;
            else
                if(grav==1)
                    temp=1;
                    while(temp<=rm-cg)
                        fincg(temp,1)=fincg(temp,1)*(targetxy/maxxydiff);
                        temp=temp+1;
                    end
                end
                if(fyes==1)
                    temp=1;
                    while(temp<=cct)
                        finc(temp,1)=finc(temp,1)*(targetxy/maxxydiff);
                        temp=temp+1;
                    end
                end
                if(myes==1)
                    temp=1;
                    while(temp<=cct)
                        minc(temp,1)=minc(temp,1)*(targetxy/maxxydiff);
                        temp=temp+1;
                    end
                end
                LocMat=oldLocMat;
            end
        elseif((maxxydiff>targetxy) || (maxphidiff>targetphi))
            if(grav==1)
                temp=1;
                while(temp<=rm-cg)
                    fincg(temp,1)=fincg(temp,1)/2;
                    temp=temp+1;
                end
            end
            if(fyes==1)
                temp=1;
                while(temp<=cct)
                    finc(temp,1)=finc(temp,1)/2;
                    temp=temp+1;
                end
            end
            if(myes==1)
                temp=1;
                while(temp<=cct)
                    minc(temp,1)=minc(temp,1)/2;
                    temp=temp+1;
                end
            end
            LocMat=oldLocMat;
        else
            nope=1;
            if(grav==1)
                temp=1;
                while(temp<=rm-cg)
                    gravsumload(temp,1)=gravsumload(temp,1)+fincg(temp,1);
                    temp=temp+1;
                end
            end
            if(fyes==1)
                temp=1;
                while(temp<=cct)
                    fmagsum(temp,1)=fmagsum(temp,1)+finc(temp,1);
                    temp=temp+1;
                end
            end
            if(myes==1)
                temp=1;
                while(temp<=cct)
                    mmagsum(temp,1)=mmagsum(temp,1)+minc(temp,1);
                    temp=temp+1;
                end
            end
            
            if(grav==1)
                loadfull=gravforce(1,1);
                load=gravsumload(1,1);
            elseif(fyes==1)
                temp=1;
                while(temp<=cct)
                    if(abs(fmag(temp,1))>0)
                        loadfull=fmag(temp,1);
                        load=fmagsum(temp,1);
                        fnumber=temp;
                        break
                    end
                    temp=temp+1;
                end
            elseif(myes==1)
                temp=1;
                while(temp<=cct)
                    if(abs(mmag(temp,1))>0)
                        loadfull=mmag(temp,1);
                        load=mmagsum(temp,1);
                        mnumber=temp;
                        break
                    end
                    temp=temp+1;
                end
            end
        end
        %%%%
        Wrench=zeros(3*(rm-cg),1);
        unitgravdir = [0 -1];
        if(grav==1)
            for i=1:rm-cg
                camforce=fincg(i,1)*unitgravdir;
                camLoc=[LocMat(i,2), LocMat(i,3)];
                Wrench((3*(i-1))+1:(3*i),1)=Wrench((3*(i-1))+1:(3*i),1)+[camforce(1,1); camforce(1,2); ((camLoc(1,1)*camforce(1,2))-(camLoc(1,2)*camforce(1,1)))];
            end
        end
        if(fyes==1)
            temp=1;
            while(temp<=cct)
                camforce=finc(temp,1)*dirMtx(temp,1:2);
                camLoc=[LocMat(nf(temp,1),2), LocMat(nf(temp,1),3)];
                Wrench((3*(nf(temp,1)-1))+1:(3*nf(temp,1)),1)=Wrench((3*(nf(temp,1)-1))+1:(3*nf(temp,1)),1)+[camforce(1,1); camforce(1,2); ((camLoc(1,1)*camforce(1,2))-(camLoc(1,2)*camforce(1,1)))];
                temp=temp+1;
            end
        end
        if(myes==1)
            temp=1;
            while(temp<=cct)
                Wrench((3*(nm(temp,1)-1))+1:(3*nm(temp,1)),1)=Wrench((3*(nm(temp,1)-1))+1:(3*nm(temp,1)),1)+[0; 0; minc(temp,1)];
                temp=temp+1;
            end
        end
    end
    %%% Draw an updated image
    if((100*abs(load)/abs(loadfull))>=(100*count/frames))
        drawCRAM(FullJointMat,LocMat,jointnum, C, boltr, t, W, dcut1, dcut2, Ec, v)
        axis equal
        xlim(xl); ylim(yl);
        
        cap = getframe(handles.axes);
        image = cap.cdata;
        [handles.LastGIF(:,:,1,count)] = grayscale(image);
        
        set(handles.text_frames, 'String', ['Frames: ' num2str(count) ' of'])
        count = count+1;
    end
    
    if get(handles.step7_stop, 'Value') == 1
        set(handles.step7_stop, 'Value', 0)
        set(handles.step7_stop, 'Enable', 'Off')
        
        set(handles.text_status, 'String', 'Simulation stopped')
        
        str = ['Simulation stopped with ' num2str(count-1) ' frames rendered. Click "Save GIF" to save the animation file. If the results seem strange, you can change the simulation settings by clicking "Change simulation settings". Click "Dismiss" to close the window.'];
        Choice = questdlg(str, 'Simulation complete', 'Save GIF', 'Change simulation settings', 'Dismiss', 'Dismiss');
        switch Choice
            case 'Save GIF'
                step7_saveGIF_Callback(hObject, eventdata, handles);
            case 'Change simulation settings'
                handles = simulation_settings_ClickedCallback(hObject, eventdata, handles);
        end
        
        guidata(hObject, handles)
        return
    end
end

set(handles.text_status, 'String', 'Simulation complete')

str = ['Simulation completed with ' num2str(count-1) ' frames rendered. Click "Save GIF" to save the animation file. If the results seem strange, you can change the simulation settings by clicking "Change simulation settings". Click "Dismiss" to close the window.'];
Choice = questdlg(str, 'Simulation complete', 'Save GIF', 'Change simulation settings', 'Dismiss', 'Dismiss');
switch Choice
    case 'Save GIF'
        step7_saveGIF_Callback(hObject, eventdata, handles);
    case 'Change simulation settings'
        handles = simulation_settings_ClickedCallback(hObject, eventdata, handles);
end

guidata(hObject, handles)
function step7_editFrames_Callback(hObject, eventdata, handles)
handles.Frames = str2double(get(handles.step7_editFrames, 'String'));
guidata(hObject, handles)
function step7_editFrames_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function handles = step7_previous_Callback(hObject, eventdata, handles)
handles = drawPatchCRAM(handles);
drawForces(handles);
drawGrounds(handles);
UpdateAxes(handles);
handles.step = 6;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)
function step7_saveGIF_Callback(hObject, eventdata, handles)
[filename, pathname] = uiputfile('*.gif', 'Save GIF animation as');
if isequal(filename,0) || isequal(pathname,0)
    return
end

prompts = {'Delay time between frames (sec)', 'Loop count', 'Animate in both directions? (bool)', 'Open once saved? (bool)'};
defaults = {num2str(handles.GIF.delay), num2str(handles.GIF.loopCount), '1', '1'};
Answer = inputdlg(prompts,'Write GIF animation',1,defaults);
if isempty(Answer); return; end
handles.GIF.delay = str2double(Answer{1});
handles.GIF.loopCount = str2double(Answer{2});
bidirectional = str2double(Answer{3});
open_bool = str2double(Answer{4});

if bidirectional
    forward = handles.LastGIF(:,:,:,1:end-1);
    backward = handles.LastGIF(:,:,:,end:-1:2);
    WriteGIF = cat(4, forward, backward);
else
    WriteGIF = handles.LastGIF;
end

imwrite(WriteGIF, fullfile(pathname, filename), 'DelayTime', handles.GIF.delay, 'LoopCount', handles.GIF.loopCount)

if open_bool
    winopen(fullfile(pathname, filename));
end

guidata(hObject, handles)
function step7_stop_Callback(hObject, eventdata, handles)
function handles = step7_done_Callback(hObject, eventdata, handles)
handles.step = 8;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)

%% Step 8: Fabricate
function step8_undeformed_Callback(hObject, eventdata, handles)
FullJointMat = handles.FullJointMat;
LocMat = handles.LocMat;
jointnum = handles.jointnum;
C = handles.C;
boltr = handles.boltr;
t = handles.t;
del = handles.del;
dcut1 = handles.dcut1;
dcut2 = handles.dcut2;

rm = size(LocMat,1);

figure;
drawlayerundeformed(FullJointMat,LocMat,jointnum, rm, C, boltr, t, del, dcut1, dcut2);
xlim(xlim*1.1)
ylim(ylim*1.1)
axis equal
set(gcf, 'color', 'w')
xlabel('Meters')
ylabel('Meters')
title('Unassembled CRAM layer')
function handles = step8_previous_Callback(hObject, eventdata, handles)
handles.step = 7;
handles = UpdateDisplay(hObject, handles, handles.step);
guidata(hObject, handles)

%% Other functions
function menu_saveScreenshot_ClickedCallback(hObject, eventdata, handles)
[filename, pathname] = uiputfile({'*.png'; '*.tiff'; '*.bmp'; '*.jpg'}, 'Save screenshot as');
if isequal(filename,0) || isequal(pathname,0)
    return;
end

cap = getframe(handles.axes);
image = cap.cdata;

imwrite(image, fullfile(pathname, filename));
function bstraplattice(Loc,Lsvec,t,extra,rotangle,veccam1)
% This funtion draws a rectangular strap of thickness t that begins at Loc
% and points in the direction and length of Lsvec. It also adds an extra
% length to both sides of length extra.

Loc(1,1)=-Loc(1,1);
Lsvec(1,1)=-Lsvec(1,1);
Ls=sqrt(dot(Lsvec,Lsvec));
unitLs=Lsvec/Ls;
perpLs=[unitLs(1,2), -unitLs(1,1)];
strapx=[(Loc(1,1)+((t/2)*perpLs(1,1))-(extra*unitLs(1,1))), (Loc(1,1)-((t/2)*perpLs(1,1))-(extra*unitLs(1,1))), (Loc(1,1)-((t/2)*perpLs(1,1))+((Ls+extra)*unitLs(1,1))), (Loc(1,1)+((t/2)*perpLs(1,1))+((Ls+extra)*unitLs(1,1)))];
strapy=[(Loc(1,2)+((t/2)*perpLs(1,2))-(extra*unitLs(1,2))), (Loc(1,2)-((t/2)*perpLs(1,2))-(extra*unitLs(1,2))), (Loc(1,2)-((t/2)*perpLs(1,2))+((Ls+extra)*unitLs(1,2))), (Loc(1,2)+((t/2)*perpLs(1,2))+((Ls+extra)*unitLs(1,2)))];
strapz = zeros(size(strapy));

%rotate them back to the global coordinate system
xunitg=(strapx*cos(rotangle))-(strapy*sin(rotangle))+veccam1(1,1);
yunitg=(strapx*sin(rotangle))+(strapy*cos(rotangle))+veccam1(2,1);
strapx=xunitg;
strapy=yunitg;

patch(strapx, strapy, strapz, [0.8 0.8 0.8], 'LineStyle', 'none');
function [stiffmat]=BuildStiffnessMatrix(FullJointMat, LocMat, jointnum, rm, cg, t, W, del, dcut1, dcut2, O, Et, Ec, v, fric, incder)
% This function builds the lattice's stiffness matrix

Rave=mean(LocMat(:,4));
incderth=incder/Rave;

stiffmat=zeros(3*(rm-cg),3*(rm-cg));
for j=1:(rm-cg)
    for k=1:(2*jointnum)
        if(FullJointMat(k,2)==j)
            phi1=LocMat(FullJointMat(k,1),1);
            phi2=LocMat(FullJointMat(k,2),1);
            veccam1=[LocMat(FullJointMat(k,1),2); LocMat(FullJointMat(k,1),3)];
            veccam2=[LocMat(FullJointMat(k,2),2); LocMat(FullJointMat(k,2),3)];
            difvec=veccam2-veccam1;
            rotangle=FullJointMat(k,3)-(pi/2)+phi1;
            RotMat=[cos(rotangle), -sin(rotangle); sin(rotangle), cos(rotangle)];
            rotvec=RotMat\difvec;
            Rp1=FullJointMat(k,6);
            Rp2=FullJointMat(k,7);
            Rb1=Rp1-(t/2);
            Rb2=Rp2-(t/2);
            beta1=FullJointMat(k,4);
            beta2=FullJointMat(k,5);
            phi=cutdec(phi2-phi1,10);
            xcom=cutdec(rotvec(1,1),10);
            ycom=cutdec(rotvec(2,1),10);
            
            [PitchLoc]=PitchLocation(phi, xcom, ycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v);
            PitchLocRot=(RotMat*[PitchLoc(1,1); PitchLoc(1,2)]);
            ActualPitchLoc=PitchLocRot+veccam1;
            D=sqrt(dot(PitchLocRot,PitchLocRot));
            unitPitchvec=PitchLocRot/D;
            if(unitPitchvec(2,1)>=0)
                anjoint=acos(dot([1; 0], unitPitchvec));
            else
                anjoint=(2*pi)-acos(dot([1; 0], unitPitchvec));
            end
            n1xa=cos(anjoint-(pi/2));
            n1ya=sin(anjoint-(pi/2));
            n2xa=cos(anjoint);
            n2ya=sin(anjoint);
            transtwist1=[1, 0, 0; veccam1(2,1), n1xa, n2xa; -veccam1(1,1), n1ya, n2ya];
            transtwist2=[1, 0, 0; ActualPitchLoc(2,1), n1xa, n2xa; -ActualPitchLoc(1,1), n1ya, n2ya];
            
            CorrectM=[1, 0, 0; -D, 1, 0; 0, 0, 1];
            
            n1xb=cos(rotangle);
            n1yb=sin(rotangle);
            n2xb=cos(rotangle+(pi/2));
            n2yb=sin(rotangle+(pi/2));
            dispan=acos(dot(unitPitchvec,[n1xb; n1yb]));
            
            difff=ActualPitchLoc-veccam2;
            L2toPitch=RotMat\difff;
            n1xc=cos(dispan-(pi/2));
            n1yc=sin(dispan-(pi/2));
            n2xc=cos(dispan);
            n2yc=sin(dispan);
            transtangent=[1, 0, 0; L2toPitch(2,1), n1xc, n2xc; -L2toPitch(1,1), n1yc, n2yc];
            Tphi=[incderth; 0; 0];
            Txcom=[0; incder; 0];
            Tycom=[0; 0; incder];
            Tphitran=transtangent*Tphi;
            Txcomtran=transtangent*Txcom;
            Tycomtran=transtangent*Tycom;
            
            transwrench=[n1xc, n2xc, 0; n1yc, n2yc, 0; ((PitchLoc(1,1)*n1yc)-(PitchLoc(1,2)*n1xc)), ((PitchLoc(1,1)*n2yc)-(PitchLoc(1,2)*n2xc)), 1];
            [Wtot]=jointloads(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, O, Et, Ec, v, fric);
            Wtot=transwrench\Wtot;
            
            [Wfinalphi1]=jointloads(phi+Tphitran(1,1), xcom+Tphitran(2,1), ycom+Tphitran(3,1), Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, O, Et, Ec, v, fric);
            Wfinalphi1=transwrench\Wfinalphi1;
            [Wfinalxcom1]=jointloads(phi+Txcomtran(1,1), xcom+Txcomtran(2,1), ycom+Txcomtran(3,1), Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, O, Et, Ec, v, fric);
            Wfinalxcom1=transwrench\Wfinalxcom1;
            [Wfinalycom1]=jointloads(phi+Tycomtran(1,1), xcom+Tycomtran(2,1), ycom+Tycomtran(3,1), Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, O, Et, Ec, v, fric);
            Wfinalycom1=transwrench\Wfinalycom1;
            
            [Wfinalphi2]=jointloads(phi-Tphitran(1,1), xcom-Tphitran(2,1), ycom-Tphitran(3,1), Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, O, Et, Ec, v, fric);
            Wfinalphi2=transwrench\Wfinalphi2;
            [Wfinalxcom2]=jointloads(phi-Txcomtran(1,1), xcom-Txcomtran(2,1), ycom-Txcomtran(3,1), Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, O, Et, Ec, v, fric);
            Wfinalxcom2=transwrench\Wfinalxcom2;
            [Wfinalycom2]=jointloads(phi-Tycomtran(1,1), xcom-Tycomtran(2,1), ycom-Tycomtran(3,1), Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, O, Et, Ec, v, fric);
            Wfinalycom2=transwrench\Wfinalycom2;
            
            dermat=zeros(3,3);
            dermat(:,1)=(((Wfinalphi1-Wtot)/incderth)+((Wtot-Wfinalphi2)/incderth))/2;
            dermat(:,2)=(((Wfinalxcom1-Wtot)/incder)+((Wtot-Wfinalxcom2)/incder))/2;
            dermat(:,3)=(((Wfinalycom1-Wtot)/incder)+((Wtot-Wfinalycom2)/incder))/2;
            stiff=abs(dermat);
            
            globalwrench=[n1xb, n2xb, 0; n1yb, n2yb, 0; ((veccam1(1,1)*n1yb)-(veccam1(2,1)*n1xb)), ((veccam1(1,1)*n2yb)-(veccam1(2,1)*n2xb)), 1];
            
            stiffmat((3*(j-1))+1:(3*j),(3*(j-1))+1:(3*j))=stiffmat((3*(j-1))+1:(3*j),(3*(j-1))+1:(3*j))+(globalwrench*transwrench*stiff*pinv(transtwist2));
            if(FullJointMat(k,1)<=(rm-cg))
                a=FullJointMat(k,1);
                stiffmat((3*(j-1))+1:(3*j),(3*(a-1))+1:(3*a))=stiffmat((3*(j-1))+1:(3*j),(3*(a-1))+1:(3*a))-(globalwrench*transwrench*stiff*CorrectM*pinv(transtwist1));
            end
        end
    end
    
    %%% This portion sums up additional collision stiffnesses on cams that bump into
    %%% other cams to which they are not joined by straps. Such stiffness prevents
    %%% the cams from passing through each other as their lattice is deformed
    
    [YoN, culpritvec]=overlap(LocMat);
    if(YoN==1)
        for i=1:rm
            gonogo=1;
            for n=1:(2*jointnum)
                if(FullJointMat(n,1)==j && FullJointMat(n,2)==i)
                    gonogo=0;
                end
            end
            if(gonogo==1 && culpritvec(j,i)>0)
                veccam1=[LocMat(i,2); LocMat(i,3)];
                veccam2=[LocMat(j,2); LocMat(j,3)];
                difvec=veccam2-veccam1;
                Rp1=LocMat(i,4);
                Rp2=LocMat(j,4);
                Rp1=Rp1-(t/2);
                Rp2=Rp2-(t/2);
                xcom=difvec(1,1);
                ycom=difvec(2,1);
                D=sqrt(dot(difvec,difvec));
                unitdifvec=difvec/D;
                if(unitdifvec(2,1)>=0)
                    anjoint=acos(dot([1; 0], unitdifvec));
                else
                    anjoint=(2*pi)-acos(dot([1; 0], unitdifvec));
                end
                n1xa=cos(anjoint-(pi/2));
                n1ya=sin(anjoint-(pi/2));
                n2xa=cos(anjoint);
                n2ya=sin(anjoint);
                transtwist1=[1, 0, 0; veccam1(2,1), n1xa, n2xa; -veccam1(1,1), n1ya, n2ya];
                transtwist2=[1, 0, 0; veccam2(2,1), n1xa, n2xa; -veccam2(1,1), n1ya, n2ya];
                
                CorrectM=[1, 0, 0; -D, 1, 0; 0, 0, 1];
                
                n1xc=cos((pi/2)-anjoint);
                n1yc=sin((pi/2)-anjoint);
                n2xc=cos(pi-anjoint);
                n2yc=sin(pi-anjoint);
                transtangent=[1, 0, 0; 0, n1xc, n2xc; 0, n1yc, n2yc];
                Txcom=[0; incder; 0];
                Tycom=[0; 0; incder];
                Txcomtran=transtangent\Txcom;
                Tycomtran=transtangent\Tycom;
                
                transwrench=[n1xa, n2xa, 0; n1ya, n2ya, 0; ((xcom*n1ya)-(ycom*n1xa)), ((xcom*n2ya)-(ycom*n2xa)), 1];
                [Wtot]=contactloads(xcom, ycom, Rp1, Rp2, W, O, Ec, v, incder);
                Wtot=transwrench\Wtot;
                
                [Wfinalxcom1]=contactloads(xcom+Txcomtran(2,1), ycom+Txcomtran(3,1), Rp1, Rp2, W, O, Ec, v, incder);
                Wfinalxcom1=transwrench\Wfinalxcom1;
                [Wfinalycom1]=contactloads(xcom+Tycomtran(2,1), ycom+Tycomtran(3,1), Rp1, Rp2, W, O, Ec, v, incder);
                Wfinalycom1=transwrench\Wfinalycom1;
                
                [Wfinalxcom2]=contactloads(xcom-Txcomtran(2,1), ycom-Txcomtran(3,1), Rp1, Rp2, W, O, Ec, v, incder);
                Wfinalxcom2=transwrench\Wfinalxcom2;
                [Wfinalycom2]=contactloads(xcom-Tycomtran(2,1), ycom-Tycomtran(3,1), Rp1, Rp2, W, O, Ec, v, incder);
                Wfinalycom2=transwrench\Wfinalycom2;
                
                dermat=zeros(3,3);
                dermat(:,1)=[0; 0; 0];
                dermat(:,2)=(((Wfinalxcom1-Wtot)/incder)+((Wtot-Wfinalxcom2)/incder))/2;
                dermat(:,3)=(((Wfinalycom1-Wtot)/incder)+((Wtot-Wfinalycom2)/incder))/2;
                
                stiff=abs(dermat);
                
                globalwrench=[1, 0, 0; 0, 1, 0; -veccam1(2,1), veccam1(1,1), 1];
                
                stiffmat((3*(j-1))+1:(3*j),(3*(j-1))+1:(3*j))=stiffmat((3*(j-1))+1:(3*j),(3*(j-1))+1:(3*j))+(globalwrench*transwrench*stiff*pinv(transtwist2));
                if(i<=(rm-cg))
                    stiffmat((3*(j-1))+1:(3*j),(3*(i-1))+1:(3*i))=stiffmat((3*(j-1))+1:(3*j),(3*(i-1))+1:(3*i))-(globalwrench*transwrench*stiff*CorrectM*pinv(transtwist1));
                end
            end
        end
    end
end
function bwrapstrapconnectlattice(x,y,r,t,Z,beta,C,phi,rotangle,veccam1)
% draws back layer strap sector where x and y are center coordinates and r is
% the base circle radius and t is the strap thickness. Z is the angle
% from the horizontal to the perpendicular line of the joint. beta is defined
% in the figure, C is how many strap thickness lengths the strap attaches to
% the cams, and phi is how much the cam rotates

inc=100;
connectprt=(C*t/r);

% connector sector points
an1=Z-beta+connectprt+phi;
an2=Z-beta+phi;
thc = an1:-abs(an1-an2)/inc:an2;
xunit1c = (r+t)*cos(thc);
xunitc = [xunit1c, 0]+x;
yunit1c = (r+t)*sin(thc);
yunitc = [yunit1c, 0]+y;

%reflect these about the correct axis
ang=(Z-(pi/2));
refxunitc=(xunitc*cos(2*ang))+(yunitc*sin(2*ang));
refyunitc=(xunitc*sin(2*ang))+(yunitc*(-cos(2*ang)));
xunitc=refxunitc;
yunitc=refyunitc;
zunitc=zeros(size(yunitc));

%rotate them back to the global coordinate system
xunitg=(xunitc*cos(rotangle))-(yunitc*sin(rotangle))+veccam1(1,1);
yunitg=(xunitc*sin(rotangle))+(yunitc*cos(rotangle))+veccam1(2,1);
xunitc=xunitg;
yunitc=yunitg;

patch(xunitc, yunitc, zunitc, [0.8 0.8 0.8], 'LineStyle', 'none');
function bwrapstrapthetalattice(x,y,r,t,Z,dcut,beta,theta,phi,rotangle,veccam1)
% draws back layer strap sector where x and y are center coordinates and r is
% the base circle radius and t is the strap thickness. Z is the angle
% from the horizontal to the perpendicular line of the joint. dcut is the
% fabricated cut by the strap that defines alpha. beta and theta are defined
% in the figure, and phi is how much the cam rotates

inc=100;
alpha=acos((r-dcut)/r);

% theta sector points
an3=Z-beta-alpha+phi;
an4=Z-beta-alpha-theta+phi;
ththe = an3:-abs(an3-an4)/inc:an4;
xunit1the = (r+t)*cos(ththe);
xunitthe = [xunit1the, 0]+x;
yunit1the = (r+t)*sin(ththe);
yunitthe = [yunit1the, 0]+y;

%reflect these about the correct axis
ang=(Z-(pi/2));
refxunitthe=(xunitthe*cos(2*ang))+(yunitthe*sin(2*ang));
refyunitthe=(xunitthe*sin(2*ang))+(yunitthe*(-cos(2*ang)));
xunitthe=refxunitthe;
yunitthe=refyunitthe;
zunitthe=zeros(size(yunitthe));

%rotate them back to the global coordinate system
xunitg=(xunitthe*cos(rotangle))-(yunitthe*sin(rotangle))+veccam1(1,1);
yunitg=(xunitthe*sin(rotangle))+(yunitthe*cos(rotangle))+veccam1(2,1);
xunitthe=xunitg;
yunitthe=yunitg;

patch(xunitthe, yunitthe, zunitthe, [0.8 0.8 0.8], 'LineStyle', 'none');
function circlelattice(x,y,r,boltr,phi,rotangle,veccam1)
% draws circle where x and y are center coordinates and r is the radius
% boltr is the radius of the bolt and phi is how much the circle rotates

inc=50;
th = 0:pi/inc:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
zunit = zeros(size(yunit));
%rotate them back to the global coordinate system
xunitg=(xunit*cos(rotangle))-(yunit*sin(rotangle))+veccam1(1,1);
yunitg=(xunit*sin(rotangle))+(yunit*cos(rotangle))+veccam1(2,1);
xunit=xunitg;
yunit=yunitg;
patch(xunit, yunit, zunit, [.6 .6 .6], 'LineStyle', 'none');

xunit1 = (boltr*cos(th)) + x + ((2*r/3)*cos(phi));
yunit1 = (boltr*sin(th)) + y + ((2*r/3)*sin(phi));
%rotate them back to the global coordinate system
xunitg=(xunit1*cos(rotangle))-(yunit1*sin(rotangle))+veccam1(1,1);
yunitg=(xunit1*sin(rotangle))+(yunit1*cos(rotangle))+veccam1(2,1);
xunit1=xunitg;
yunit1=yunitg;
patch(xunit1, yunit1, 'k', 'LineStyle', 'none');

xunit2 = (boltr*cos(th)) + x + ((2*r/3)*cos(phi-(pi/2)));
yunit2 = (boltr*sin(th)) + y + ((2*r/3)*sin(phi-(pi/2)));
%rotate them back to the global coordinate system
xunitg=(xunit2*cos(rotangle))-(yunit2*sin(rotangle))+veccam1(1,1);
yunitg=(xunit2*sin(rotangle))+(yunit2*cos(rotangle))+veccam1(2,1);
xunit2=xunitg;
yunit2=yunitg;
patch(xunit2, yunit2, 'k', 'LineStyle', 'none');

xunit3 = (boltr*cos(th)) + x - ((2*r/3)*cos(phi));
yunit3 = (boltr*sin(th)) + y - ((2*r/3)*sin(phi));
%rotate them back to the global coordinate system
xunitg=(xunit3*cos(rotangle))-(yunit3*sin(rotangle))+veccam1(1,1);
yunitg=(xunit3*sin(rotangle))+(yunit3*cos(rotangle))+veccam1(2,1);
xunit3=xunitg;
yunit3=yunitg;
patch(xunit3, yunit3, 'k', 'LineStyle', 'none');

xunit4 = (boltr*cos(th)) + x - ((2*r/3)*cos(phi-(pi/2)));
yunit4 = (boltr*sin(th)) + y - ((2*r/3)*sin(phi-(pi/2)));
%rotate them back to the global coordinate system
xunitg=(xunit4*cos(rotangle))-(yunit4*sin(rotangle))+veccam1(1,1);
yunitg=(xunit4*sin(rotangle))+(yunit4*cos(rotangle))+veccam1(2,1);
xunit4=xunitg;
yunit4=yunitg;
patch(xunit4, yunit4, 'k', 'LineStyle', 'none');
function circlelattice2(x,y,r,boltr,phi,rotangle,veccam1)
% draws circle where x and y are center coordinates and r is the radius
% boltr is the radius of the bolt and phi is how much the circle rotates

format long

hold on
inc=50;
th = 0:pi/inc:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
zunit = zeros(size(yunit));
%rotate them back to the global coordinate system
xunitg=(xunit*cos(rotangle))-(yunit*sin(rotangle))+veccam1(1,1);
yunitg=(xunit*sin(rotangle))+(yunit*cos(rotangle))+veccam1(2,1);
xunit=xunitg;
yunit=yunitg;
patch(xunit, yunit, zunit, [.6 .6 .6], 'LineStyle', 'none');

xunit1 = (boltr*cos(th)) + x + ((2*r/3)*cos(phi));
yunit1 = (boltr*sin(th)) + y + ((2*r/3)*sin(phi));
%rotate them back to the global coordinate system
xunitg=(xunit1*cos(rotangle))-(yunit1*sin(rotangle))+veccam1(1,1);
yunitg=(xunit1*sin(rotangle))+(yunit1*cos(rotangle))+veccam1(2,1);
xunit1=xunitg;
yunit1=yunitg;
patch(xunit1, yunit1, [1 1 1], 'LineStyle', 'none');

xunit2 = (boltr*cos(th)) + x + ((2*r/3)*cos(phi-(pi/2)));
yunit2 = (boltr*sin(th)) + y + ((2*r/3)*sin(phi-(pi/2)));
%rotate them back to the global coordinate system
xunitg=(xunit2*cos(rotangle))-(yunit2*sin(rotangle))+veccam1(1,1);
yunitg=(xunit2*sin(rotangle))+(yunit2*cos(rotangle))+veccam1(2,1);
xunit2=xunitg;
yunit2=yunitg;
patch(xunit2, yunit2, [1 1 1], 'LineStyle', 'none');

xunit3 = (boltr*cos(th)) + x - ((2*r/3)*cos(phi));
yunit3 = (boltr*sin(th)) + y - ((2*r/3)*sin(phi));
%rotate them back to the global coordinate system
xunitg=(xunit3*cos(rotangle))-(yunit3*sin(rotangle))+veccam1(1,1);
yunitg=(xunit3*sin(rotangle))+(yunit3*cos(rotangle))+veccam1(2,1);
xunit3=xunitg;
yunit3=yunitg;
patch(xunit3, yunit3, [1 1 1], 'LineStyle', 'none');

xunit4 = (boltr*cos(th)) + x - ((2*r/3)*cos(phi-(pi/2)));
yunit4 = (boltr*sin(th)) + y - ((2*r/3)*sin(phi-(pi/2)));
%rotate them back to the global coordinate system
xunitg=(xunit4*cos(rotangle))-(yunit4*sin(rotangle))+veccam1(1,1);
yunitg=(xunit4*sin(rotangle))+(yunit4*cos(rotangle))+veccam1(2,1);
xunit4=xunitg;
yunit4=yunitg;
patch(xunit4, yunit4, [1 1 1], 'LineStyle', 'none');

hold off
function [Wload, Ttan]=compressionregime1(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta2, vec1, vec2, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 1

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp2*sin(alpha2/2))-(Rp2*theta2);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=(2*Rp2*sin(alpha2/2))+(Rp2*((exp(fric*theta2)-1)/fric))+(Li*exp(fric*theta2));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
R2vec=vec2-[xcom, ycom];
perpvec=-[R2vec(1,2), -R2vec(1,1)];
unitperp=(perpvec)/sqrt(dotFast((perpvec),(perpvec)));
forcecorner2=Talpha2*unitperp;
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

% sum of reaction normal forces from the strap on the cam
dd=phi-beta2-alpha2;
sd=sin(dd);
cd=cos(dd);
sdt=sin(dd-real(theta2));
cdt=cos(dd-real(theta2));
Ncon=Talpha2*exp(fric*dd)*(2/(1+(fric^2)));
NFx=Ncon*(((sd/(2*exp(fric*dd)))-(fric*cd/(2*exp(fric*dd))))-((sdt/(2*exp(fric*(dd-real(theta2)))))-(fric*cdt/(2*exp(fric*(dd-real(theta2)))))));
NFy=Ncon*(((-fric*sd/(2*exp(fric*dd)))-(cd/(2*exp(fric*dd))))-((-fric*sdt/(2*exp(fric*(dd-real(theta2)))))-(cdt/(2*exp(fric*(dd-real(theta2)))))));
NMz=(xcom*NFy)-(ycom*NFx);
Wforcenormal=[NFx; NFy; NMz];

% sum of reaction tangential friction forces from the strap on the cam
fFx=-fric*NFy;
fFy=fric*NFx;
fMz=(fFy*xcom)-(fFx*ycom)+(Rb2*Talpha2*(exp(fric*real(theta2))-1));
Wforcetan=[fFx; fFy; fMz];

Wforce=Wcomp+Wforceattach+Wforcecorner+Wforcenormal+Wforcetan;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % moment to resist the strap bent around the cams
% Momenttheta2=(E*W*(t^3))/(12*Rp2);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     rotbeta1=(2*pi)-acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% else
%     rotbeta1=acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% end
% if(Talpha2>0)
%     Momentbeta1=LivingHingeK*rotbeta1;
%     Momentalpha2=LivingHingeK*(alpha2/2);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentbeta1=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta2+Momentbeta1+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime2(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, vec1, vec2, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 2

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls+(2*Rp2*sin(alpha2/2))-Lo,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=Lo;
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
if(Ls==0)
    perpLoc=[Locforce2(1,2), -Locforce2(1,1)];
    forcecorner2=Talpha2*(perpLoc/sqrt(dotFast(perpLoc,perpLoc)));
else
    forcecorner2=Talpha2*(Lsvec/Ls);
end
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

Wforce=Wcomp+Wforceattach+Wforcecorner;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     rotbeta1=(2*pi)-acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% else
%     rotbeta1=acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% end
% slopevec=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
% anbetween=acos(dotFast(slopevec,Lsvec/Ls));
% if(Talpha2>0)
%     Momentbeta1=LivingHingeK*(rotbeta1);
%     Momentalpha2=LivingHingeK*(anbetween);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentbeta1=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momentbeta1+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime3(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, Loc, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 3

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Lo,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=Lo;
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap
Locforce=Loc+Lsvec;
force=Talpha2*(Lsvec/Ls);
Wten=[force(1,1); force(1,2); ((Locforce(1,1)*force(1,2))-(Locforce(1,2)*force(1,1)))]; % tension force from strap

Wforce=Wcomp+Wten;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     rotbeta1=(2*pi)-acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% else
%     rotbeta1=acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% end
% vec1=Loc+Lsvec;
% R2vec=vec1-[xcom, ycom];
% perpvec=[R2vec(1,2), -R2vec(1,1)];
% unitperp=perpvec/sqrt(dotFast(perpvec,perpvec));
% if(unitperp(1,2)>=0)
%     anperp=acos(dotFast([1,0], unitperp));
% else
%     anperp=-acos(dotFast([1,0], unitperp));
% end
% unitLs=-Lsvec/Ls;
% if(unitLs(1,2)>=0)
%     anLs=acos(dotFast([1,0], unitLs));
% else
%     anLs=-acos(dotFast([1,0], unitLs));
% end
% anbetween=anLs-anperp;
%
% if(Talpha2>0)
%     Momentbeta1=LivingHingeK*(rotbeta1);
%     Momentbeta2=LivingHingeK*(-anbetween);
% else
%     Momentbeta1=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momentbeta1+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime4(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta2, vec1, vec2, Loc, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 4

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp2*sin(alpha2/2))-(Rp2*theta2);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=(2*Rp2*sin(alpha2/2))+(Rp2*((exp(fric*theta2)-1)/fric))+(Li*exp(fric*theta2));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
R2vec=vec2-[xcom, ycom];
perpvec=-[R2vec(1,2), -R2vec(1,1)];
unitperp=(perpvec)/sqrt(dotFast((perpvec),(perpvec)));
forcecorner2=Talpha2*unitperp;
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

% sum of reaction normal forces from the strap on the cam
dd=phi-beta2-alpha2;
sd=sin(dd);
cd=cos(dd);
sdt=sin(dd-real(theta2));
cdt=cos(dd-real(theta2));
Ncon=Talpha2*exp(fric*dd)*(2/(1+(fric^2)));
NFx=Ncon*(((sd/(2*exp(fric*dd)))-(fric*cd/(2*exp(fric*dd))))-((sdt/(2*exp(fric*(dd-real(theta2)))))-(fric*cdt/(2*exp(fric*(dd-real(theta2)))))));
NFy=Ncon*(((-fric*sd/(2*exp(fric*dd)))-(cd/(2*exp(fric*dd))))-((-fric*sdt/(2*exp(fric*(dd-real(theta2)))))-(cdt/(2*exp(fric*(dd-real(theta2)))))));
NMz=(xcom*NFy)-(ycom*NFx);
Wforcenormal=[NFx; NFy; NMz];

% sum of reaction tangential friction forces from the strap on the cam
fFx=-fric*NFy;
fFy=fric*NFx;
fMz=(fFy*xcom)-(fFx*ycom)+(Rb2*Talpha2*(exp(fric*real(theta2))-1));
Wforcetan=[fFx; fFy; fMz];

% reaction force at bend corner from the strap
Locbend=Loc+Lsvec;
R2bend=Locbend-[xcom, ycom];
dircurve=[R2bend(1,2), -R2bend(1,1)];
unitdircurve=dircurve/sqrt(dotFast(dircurve,dircurve));
forcebend1=Talpha2*exp(fric*theta2)*unitdircurve;
forcebend2=Talpha2*exp(fric*theta2)*(Lsvec/Ls);
forcebend=forcebend1+forcebend2;
Wforcebend=[forcebend(1,1); forcebend(1,2); ((Locbend(1,1)*forcebend(1,2))-(Locbend(1,2)*forcebend(1,1)))];

Wforce=Wcomp+Wforceattach+Wforcecorner+Wforcenormal+Wforcetan+Wforcebend;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % moment to resist the strap bent around the cams
% Momenttheta2=(E*W*(t^3))/(12*Rp2);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=0;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% anLs=acos(dotFast([1,0],Lsvec/Ls));
% rotbeta1=((pi/2)-beta1)-anLs;
% R2vec=(Loc+Lsvec)-[xcom, ycom];
% perpvec=[R2vec(1,2), -R2vec(1,1)];
% unitperp=perpvec/sqrt(dotFast(perpvec,perpvec));
% unitLs=-Lsvec/Ls;
% anbetween=acos(dotFast(unitperp,unitLs));
%
% if(Talpha2>0)
%     Momentbeta1=LivingHingeK*(-rotbeta1);
%     Momentbend=LivingHingeK*(anbetween);
%     Momentalpha2=LivingHingeK*(alpha2/2);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentbeta1=0;
%     Momentbend=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta2+Momentbeta1+Momentbend+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime5(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta2, attach1, vec1, vec2, Loc, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 5

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp2*sin(alpha2/2))-(Rp2*theta2);
Ls=sqrt(dotFast(Lsvec,Lsvec));
strapmag=sqrt(dotFast((Loc-attach1),(Loc-attach1)));
stretch=cutdec(Ls+strapmag-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=(2*Rp2*sin(alpha2/2))+(Rp2*((exp(fric*theta2)-1)/fric))+(Li*exp(fric*theta2));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
R2vec=vec2-[xcom, ycom];
perpvec=-[R2vec(1,2), -R2vec(1,1)];
unitperp=(perpvec)/sqrt(dotFast((perpvec),(perpvec)));
forcecorner2=Talpha2*unitperp;
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

% sum of reaction normal forces from the strap on the cam
dd=phi-beta2-alpha2;
sd=sin(dd);
cd=cos(dd);
sdt=sin(dd-real(theta2));
cdt=cos(dd-real(theta2));
Ncon=Talpha2*exp(fric*dd)*(2/(1+(fric^2)));
NFx=Ncon*(((sd/(2*exp(fric*dd)))-(fric*cd/(2*exp(fric*dd))))-((sdt/(2*exp(fric*(dd-real(theta2)))))-(fric*cdt/(2*exp(fric*(dd-real(theta2)))))));
NFy=Ncon*(((-fric*sd/(2*exp(fric*dd)))-(cd/(2*exp(fric*dd))))-((-fric*sdt/(2*exp(fric*(dd-real(theta2)))))-(cdt/(2*exp(fric*(dd-real(theta2)))))));
NMz=(xcom*NFy)-(ycom*NFx);
Wforcenormal=[NFx; NFy; NMz];

% sum of reaction tangential friction forces from the strap on the cam
fFx=-fric*NFy;
fFy=fric*NFx;
fMz=(fFy*xcom)-(fFx*ycom)+(Rb2*Talpha2*(exp(fric*real(theta2))-1));
Wforcetan=[fFx; fFy; fMz];

% reaction force at bend corner from the strap
Locbend=Loc+Lsvec;
R2bend=Locbend-[xcom, ycom];
dircurve=[R2bend(1,2), -R2bend(1,1)];
unitdircurve=dircurve/sqrt(dotFast(dircurve,dircurve));
forcebend1=Talpha2*exp(fric*theta2)*unitdircurve;
forcebend2=Talpha2*exp(fric*theta2)*(Lsvec/Ls);
forcebend=forcebend1+forcebend2;
Wforcebend=[forcebend(1,1); forcebend(1,2); ((Locbend(1,1)*forcebend(1,2))-(Locbend(1,2)*forcebend(1,1)))];

Wforce=Wcomp+Wforceattach+Wforcecorner+Wforcenormal+Wforcetan+Wforcebend;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % moment to resist the strap bent around the cams
% Momenttheta2=(E*W*(t^3))/(12*Rp2);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% anstrap=acos(dotFast([1,0],(Loc-attach1)/sqrt(dotFast((Loc-attach1),(Loc-attach1)))));
% rotbeta1=((pi/2)-beta1)-anstrap;
% ancomp=acos(dotFast([1,0],(Lsvec/Ls)));
% anbetween1=anstrap-ancomp;
% R2vec=(Loc+Lsvec)-[xcom, ycom];
% perpvec=[R2vec(1,2), -R2vec(1,1)];
% unitperp=perpvec/sqrt(dotFast(perpvec,perpvec));
% unitLs=-Lsvec/Ls;
% anbetween2=acos(dotFast(unitperp,unitLs));
%
% if(Talpha2>0)
%     Momentbeta1=LivingHingeK*(-rotbeta1);
%     Momentbend1=LivingHingeK*(-anbetween1);
%     Momentbend2=LivingHingeK*(anbetween2);
%     Momentalpha2=LivingHingeK*(alpha2/2);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentbeta1=0;
%     Momentbend1=0;
%     Momentbend2=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta2+Momentbeta1+Momentbend1+Momentbend2+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime6(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, theta2, vec1, vec2, Loc, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 6

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp1*sin(alpha1/2))-(2*Rp2*sin(alpha2/2))-(Rp1*theta1)-(Rp2*theta2);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=(2*Rp2*sin(alpha2/2))+(2*Rp1*sin(alpha1/2)*exp(fric*(theta2-theta1)))+(Rp2*((exp(fric*theta2)-1)/fric))+(Rp1*((exp(fric*theta1)-1)/fric)*exp(fric*(theta2-theta1)))+(Li*exp(fric*theta2));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
R2vec=vec2-[xcom, ycom];
perpvec=-[R2vec(1,2), -R2vec(1,1)];
unitperp=(perpvec)/sqrt(dotFast((perpvec),(perpvec)));
forcecorner2=Talpha2*unitperp;
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

% sum of reaction normal forces from the strap on the cam
dd=phi-beta2-alpha2;
sd=sin(dd);
cd=cos(dd);
sdt=sin(dd-real(theta2));
cdt=cos(dd-real(theta2));
Ncon=Talpha2*exp(fric*dd)*(2/(1+(fric^2)));
NFx=Ncon*(((sd/(2*exp(fric*dd)))-(fric*cd/(2*exp(fric*dd))))-((sdt/(2*exp(fric*(dd-real(theta2)))))-(fric*cdt/(2*exp(fric*(dd-real(theta2)))))));
NFy=Ncon*(((-fric*sd/(2*exp(fric*dd)))-(cd/(2*exp(fric*dd))))-((-fric*sdt/(2*exp(fric*(dd-real(theta2)))))-(cdt/(2*exp(fric*(dd-real(theta2)))))));
NMz=(xcom*NFy)-(ycom*NFx);
Wforcenormal=[NFx; NFy; NMz];

% sum of reaction tangential friction forces from the strap on the cam
fFx=-fric*NFy;
fFy=fric*NFx;
fMz=(fFy*xcom)-(fFx*ycom)+(Rb2*Talpha2*(exp(fric*real(theta2))-1));
Wforcetan=[fFx; fFy; fMz];

% reaction force at bend corner from the strap
Locbend=Loc+Lsvec;
R2bend=Locbend-[xcom, ycom];
dircurve=[R2bend(1,2), -R2bend(1,1)];
unitdircurve=dircurve/sqrt(dotFast(dircurve,dircurve));
forcebend1=Talpha2*exp(fric*theta2)*unitdircurve;
forcebend2=Talpha2*exp(fric*theta2)*(Lsvec/Ls);
forcebend=forcebend1+forcebend2;
Wforcebend=[forcebend(1,1); forcebend(1,2); ((Locbend(1,1)*forcebend(1,2))-(Locbend(1,2)*forcebend(1,1)))];

Wforce=Wcomp+Wforceattach+Wforcecorner+Wforcenormal+Wforcetan+Wforcebend;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % moment to resist the strap bent around the cams
% Momenttheta1=-(E*W*(t^3))/(12*Rp1);
% Momenttheta2=(E*W*(t^3))/(12*Rp2);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     ancomp=-acos(dotFast([1,0],(Lsvec/Ls)));
% else
%     ancomp=acos(dotFast([1,0],(Lsvec/Ls)));
% end
% perpvec1=[Loc(1,2),-Loc(1,1)];
% unitperp1=perpvec1/sqrt(dotFast(perpvec1,perpvec1));
% if(unitperp1(1,2)<0)
%     anperp1=-acos(dotFast([1,0],unitperp1));
% else
%     anperp1=acos(dotFast([1,0],unitperp1));
% end
% anbetween1=anperp1-ancomp;
% anR2=Loc+Lsvec-[xcom, ycom];
% perpvec2=[anR2(1,2),-anR2(1,1)];
% unitperp2=-perpvec2/sqrt(dotFast(perpvec2,perpvec2));
% if(unitperp2<0)
%     anperp2=-acos(dotFast([1,0],unitperp2));
% else
%     anperp2=acos(dotFast([1,0],unitperp2));
% end
% anbetween2=anperp2-ancomp;
%
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-alpha1/2);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentbend1=LivingHingeK*(-anbetween1);
%     Momentbend2=LivingHingeK*(anbetween2);
%     Momentalpha2=LivingHingeK*(alpha2/2);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentbend1=0;
%     Momentbend2=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta1+Momenttheta2+Momentalpha1+Momentbeta1+Momentbend1+Momentbend2+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime7(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, vec1, vec2, Loc, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 7

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp1*sin(alpha1/2))-(Rp1*theta1);
Ls=sqrt(dotFast(Lsvec,Lsvec));
strapmag=sqrt(dotFast((vec1-(Loc+Lsvec)),(vec1-(Loc+Lsvec))));
stretch=cutdec(Ls+strapmag-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch)*(exp(fric*theta1));
    denominat=(2*Rp1*sin(alpha1/2))+(Rp1*((exp(fric*theta1)-1)/fric))+(Li*exp(fric*theta1));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at bend from the strap
Locbend=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcebend1=Talpha2*dirforce;
forcebend2=Talpha2*(Lsvec/Ls);
forcebend=forcebend1+forcebend2;
Wforcebend=[forcebend(1,1); forcebend(1,2); ((Locbend(1,1)*forcebend(1,2))-(Locbend(1,2)*forcebend(1,1)))];

Wforce=Wcomp+Wforceattach+Wforcebend;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % moment to resist the strap bent around the cams
% Momenttheta1=-(E*W*(t^3))/(12*Rp1);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     ancomp=-acos(dotFast([1,0],(Lsvec/Ls)));
% else
%     ancomp=acos(dotFast([1,0],(Lsvec/Ls)));
% end
% perpvec1=[Loc(1,2),-Loc(1,1)];
% unitperp1=perpvec1/sqrt(dotFast(perpvec1,perpvec1));
% if(unitperp1(1,2)<0)
%     anperp1=-acos(dotFast([1,0],unitperp1));
% else
%     anperp1=acos(dotFast([1,0],unitperp1));
% end
% anbetween1=anperp1-ancomp;
% slopevec=vec1-(Loc+Lsvec);
% unitslopevec=slopevec/sqrt(dotFast(slopevec,slopevec));
% anbetween2=acos(dotFast(unitslopevec,Lsvec/Ls));
% R2vec=vec1-[xcom,ycom];
% perpvec=[R2vec(1,2),-R2vec(1,1)];
% unitperp=-perpvec/sqrt(dotFast(perpvec,perpvec));
% rotanbeta2=acos(dotFast(unitslopevec,unitperp));
%
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-alpha1/2);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentbend1=LivingHingeK*(-anbetween1);
%     Momentbend2=LivingHingeK*(anbetween2);
%     Momentbeta2=LivingHingeK*(rotanbeta2);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentbend1=0;
%     Momentbend2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta1+Momentalpha1+Momentbeta1+Momentbend1+Momentbend2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime8(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, vec1, Loc, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 8

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp1*sin(alpha1/2))-(Rp1*theta1);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch)*(exp(fric*theta1));
    denominat=(2*Rp1*sin(alpha1/2))+(Rp1*((exp(fric*theta1)-1)/fric))+(Li*exp(fric*theta1));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap
Locforce=Loc+Lsvec;
if(Ls==0)
    perpLoc=[Locforce(1,2), -Locforce(1,1)];
    force=Talpha2*(perpLoc/sqrt(dotFast(perpLoc,perpLoc)));
else
    force=Talpha2*(Lsvec/Ls);
end
Wforceattach=[force(1,1); force(1,2); ((Locforce(1,1)*force(1,2))-(Locforce(1,2)*force(1,1)))]; % tension force from strap

Wforce=Wcomp+Wforceattach;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % moment to resist the strap bent around the cams
% Momenttheta1=-(E*W*(t^3))/(12*Rp1);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     ancomp=-acos(dotFast([1,0],(Lsvec/Ls)));
% else
%     ancomp=acos(dotFast([1,0],(Lsvec/Ls)));
% end
% perpvec1=[Loc(1,2),-Loc(1,1)];
% unitperp1=perpvec1/sqrt(dotFast(perpvec1,perpvec1));
% if(unitperp1(1,2)<0)
%     anperp1=-acos(dotFast([1,0],unitperp1));
% else
%     anperp1=acos(dotFast([1,0],unitperp1));
% end
% anbetween1=anperp1-ancomp;
% R2vec=vec1-[xcom,ycom];
% perpvec=[R2vec(1,2),-R2vec(1,1)];
% unitperp=-perpvec/sqrt(dotFast(perpvec,perpvec));
% if(unitperp(1,2)<0)
%     anperp2=-acos(dotFast([1,0],unitperp));
% else
%     anperp2=acos(dotFast([1,0],unitperp));
% end
% rotanbeta2=anperp2-ancomp;
%
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-alpha1/2);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentbend1=LivingHingeK*(-anbetween1);
%     Momentbeta2=LivingHingeK*(rotanbeta2);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentbend1=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta1+Momentalpha1+Momentbeta1+Momentbend1+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime9(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, Loc, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 9

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp1*sin(alpha1/2))-(Rp1*theta1);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch)*(exp(fric*theta1));
    denominat=(2*Rp1*sin(alpha1/2))+(Rp1*((exp(fric*theta1)-1)/fric))+(Li*exp(fric*theta1));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap
Locforce=Loc+Lsvec;
if(Ls==0)
    perpLoc=[Locforce(1,2), -Locforce(1,1)];
    force=Talpha2*(perpLoc/sqrt(dotFast(perpLoc,perpLoc)));
else
    force=Talpha2*(Lsvec/Ls);
end
Wforceattach=[force(1,1); force(1,2); ((Locforce(1,1)*force(1,2))-(Locforce(1,2)*force(1,1)))]; % tension force from strap

Wforce=Wcomp+Wforceattach;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % moment to resist the strap bent around the cams
% Momenttheta1=-(E*W*(t^3))/(12*Rp1);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% vec1=Loc+Lsvec;
% R2vec=vec1-[xcom, ycom];
% perpvec=[R2vec(1,2), -R2vec(1,1)];
% unitperp=perpvec/sqrt(dotFast(perpvec,perpvec));
% if(unitperp(1,2)>0)
%     anperp=(2*pi)-acos(dotFast([1,0],unitperp));
% else
%     anperp=acos(dotFast([1,0],unitperp));
% end
% unitLs=-Lsvec/Ls;
% if(unitLs(1,2)>0)
%     anLs=(2*pi)-acos(dotFast([1,0],unitLs));
% else
%     anLs=acos(dotFast([1,0],unitLs));
% end
% anbetween=anLs-anperp;
%
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-alpha1/2);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentbeta2=LivingHingeK*(anbetween);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta1+Momentalpha1+Momentbeta1+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime10(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, attach1, attalpha1, Loc, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 10

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls+(2*Rp1*sin(alpha1/2))-Lo,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=Lo;
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap
Locforce=Loc+Lsvec;
force=Talpha2*(Lsvec/Ls);
Wforceattach=[force(1,1); force(1,2); ((Locforce(1,1)*force(1,2))-(Locforce(1,2)*force(1,1)))]; % tension force from strap

Wforce=Wcomp+Wforceattach;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% slopevec1=(attalpha1-attach1)/sqrt(dotFast((attalpha1-attach1),(attalpha1-attach1)));
% rotalpha1=acos(dotFast([1, 0],slopevec1))-acos(dotFast([1, 0],Lsvec/Ls));
% vec1=Loc+Lsvec;
% R2vec=vec1-[xcom, ycom];
% perpvec=[R2vec(1,2), -R2vec(1,1)];
% unitperp=perpvec/sqrt(dotFast(perpvec,perpvec));
% if(unitperp(1,2)>0)
%     anperp=(2*pi)-acos(dotFast([1,0],unitperp));
% else
%     anperp=acos(dotFast([1,0],unitperp));
% end
% unitLs=-Lsvec/Ls;
% if(unitLs(1,2)>0)
%     anLs=(2*pi)-acos(dotFast([1,0],unitLs));
% else
%     anLs=acos(dotFast([1,0],unitLs));
% end
% anbetween=anLs-anperp;
%
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-rotalpha1);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentbeta2=LivingHingeK*(anbetween);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momentalpha1+Momentbeta1+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload, Ttan]=compressionregime11(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, Loc, Lsvec, Fc)
% This function calculates the loads on a layer in the compression scenario
% within regime 11

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Lo,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=Lo;
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% compression force
Loccomp=[xcom, ycom];
comp=Fc*(-Loccomp/sqrt(dotFast(Loccomp,Loccomp)));
Wcomp=[comp(1,1); comp(1,2); ((Loccomp(1,1)*comp(1,2))-(Loccomp(1,2)*comp(1,1)))]; % compression force

% force to resist tension in the strap
Locforce=Loc+Lsvec;
force=Talpha2*(Lsvec/Ls);
Wforceattach=[force(1,1); force(1,2); ((Locforce(1,1)*force(1,2))-(Locforce(1,2)*force(1,1)))]; % tension force from strap

Wforce=Wcomp+Wforceattach;

fvec=[Wforce(1,1), Wforce(2,1), 0];
unitcenter=Loccomp/sqrt(dotFast(Loccomp,Loccomp));
unitcenter(1,3)=0;
Ttanv=-cross(unitcenter,fvec);
Ttan=Ttanv(1,3);

% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     rotbeta1=(2*pi)-acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% else
%     rotbeta1=acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% end
% vec1=Loc+Lsvec;
% R2vec=vec1-[xcom, ycom];
% perpvec=[R2vec(1,2), -R2vec(1,1)];
% unitperp=perpvec/sqrt(dotFast(perpvec,perpvec));
% if(unitperp(1,2)>0)
%     anperp=(2*pi)-acos(dotFast([1,0],unitperp));
% else
%     anperp=acos(dotFast([1,0],unitperp));
% end
% unitLs=-Lsvec/Ls;
% if(unitLs(1,2)>0)
%     anLs=(2*pi)-acos(dotFast([1,0],unitLs));
% else
%     anLs=acos(dotFast([1,0],unitLs));
% end
% anbetween=anLs-anperp;
%
% if(Talpha2>0)
%     Momentbeta1=LivingHingeK*(rotbeta1);
%     Momentbeta2=LivingHingeK*(anbetween);
% else
%     Momentbeta1=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momentbeta1+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wtot]=contactloads(xcom, ycom, Rp1, Rp2, W, O, Ec, v, incder)
%This function calculates the contact load when two cams collide that are
%not joined together by straps

centers=[xcom, ycom];
D=sqrt(dotFast(centers,centers));
if(cutdec(D-Rp1-Rp2,10)>=0)
    Fc=0;
else
    Compress=cutdec(abs(D-Rp1-Rp2),10);
    Fc=incder;
    error=1e10;
    ct=0;
    while (error > 1e-10)
        ct=ct+1;
        b=4*sqrt((2*Fc*((1-(v^2))/Ec))/(pi*W*((1/Rp1)+(1/Rp2))));
        func=(((2*Fc)/(pi*W))*((1-(v^2))/Ec)*((2/3)+(log(8*Rp1/b))+(log(8*Rp2/b))))-Compress;
        binc=4*sqrt((2*(Fc+incder)*((1-(v^2))/Ec))/(pi*W*((1/Rp1)+(1/Rp2))));
        funcinc=(((2*(Fc+incder))/(pi*W))*((1-(v^2))/Ec)*((2/3)+(log(8*Rp1/binc))+(log(8*Rp2/binc))))-Compress;
        deriv=(funcinc-func)/incder;
        newFc=Fc-(func/deriv);
        error=abs(Fc-newFc);
        Fc=newFc;
        if(ct>1e4)
            break
        end
    end
end
unitcenters=-centers/D;
force=O*Fc*unitcenters;
Wtot=[force(1,1); force(1,2); (xcom*force(1,2))-(ycom*force(1,1))];
function [num]=cutdec(X,N)

num=(round(X*(10^N)))/(10^N);
function [dofn]=DOF(FullJointMat,LocMat,jointnum,rm,cg)
% This function calculates how many degrees of freeom a CRAM achieves
% jointnum is the number of joints in the lattice
% rm is the number of circular cams
% cg is the number of grounded cams

cutoffcam=rm-(cg-1);
for a=1:jointnum
    if(FullJointMat(a,1)>cutoffcam)
        FullJointMat(a,1)=cutoffcam;
    elseif(FullJointMat(a,2)>cutoffcam)
        FullJointMat(a,2)=cutoffcam;
    end
end
for b=1:jointnum
    if(FullJointMat(b,1)==FullJointMat(b,2))
        for c=b:jointnum
            FullJointMat(b,:)=FullJointMat(b+1,:);
        end
        jointnum=jointnum-1;
    end
end
Cmat=zeros(jointnum,cutoffcam);
for i=1:jointnum
    for j=1:rm-(cg-1)
        if(FullJointMat(i,1)==j)
            Cmat(i,j)=-1;
        elseif(FullJointMat(i,2)==j)
            Cmat(i,j)=1;
        end
    end
end
Qmat=transpose(null(transpose(Cmat),'r'));
[row, col]=size(Qmat);
FT=zeros(6*row,col);
w=[0;0;1];
for k=1:jointnum
    Lc=[LocMat(FullJointMat(k,1),2); LocMat(FullJointMat(k,1),3); 0];
    Lv=Lc+[(FullJointMat(k,6)*cos(FullJointMat(k,3))); (FullJointMat(k,6)*sin(FullJointMat(k,3))); 0];
    Twist=[w; cross(Lv,w)];
    for m=0:row-1
        FT((6*m)+1:(6*m)+6,k)=Qmat(m+1,k)*Twist;
    end
end
Soln=null(FT);
[rsoln, csoln]=size(Soln);
dofn=csoln;
function product = dotFast(x, y)
product = x * y';
function drawCRAM(FullJointMat,LocMat,jointnum, C, boltr, t, W, dcut1, dcut2, Ec, v)

hold off;
plot(LocMat(1,2), LocMat(1,3), 'k.');
hold on;

%%%% draw back layer
for i=1:jointnum
    phi1=LocMat(FullJointMat(i,1),1);
    phi2=LocMat(FullJointMat(i,2),1);
    veccam2=[LocMat(FullJointMat(i,2),2); LocMat(FullJointMat(i,2),3)];
    veccam1=[LocMat(FullJointMat(i,1),2); LocMat(FullJointMat(i,1),3)];
    difvec=veccam2-veccam1;
    rotangle=FullJointMat(i,3)-(pi/2)+phi1;
    rotvec=[cos(rotangle), -sin(rotangle); sin(rotangle), cos(rotangle)]\difvec;
    Rp1=FullJointMat(i,6);
    Rp2=FullJointMat(i,7);
    Rb1=Rp1-(t/2);
    Rb2=Rp2-(t/2);
    beta1=FullJointMat(i,4);
    beta2=FullJointMat(i,5);
    phi=phi2-phi1;
    xcom=rotvec(1,1);
    ycom=rotvec(2,1);
    
    bphi=-phi;
    bxcom=-xcom;
    bycom=ycom;
    
    % back layer
    [btension, bregime, battach1, battalpha1, bvec1, bvec2, btheta1, btheta2, bLoc, bLsvec, bb, bFc]=kinematics(bphi, bxcom, bycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v);
    
    if (btension==1) % tension
        if (bregime==1)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(bxcom,bycom,Rb2,t,2*pi,dcut2,beta2,btheta2,bphi,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==2)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==3)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==4)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(bxcom,bycom,Rb2,t,2*pi,dcut2,beta2,btheta2,bphi,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==5)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==6)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==7)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,btheta1,0,rotangle,veccam1)
            bwrapstrapthetalattice(bxcom,bycom,Rb2,t,2*pi,dcut2,beta2,btheta2,bphi,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==8)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,btheta1,0,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==9)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,btheta1,0,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        end
    else % compression
        if (bregime==1)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(bxcom,bycom,Rb2,t,2*pi,dcut2,beta2,btheta2,bphi,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==2)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==3)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==4)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(bxcom,bycom,Rb2,t,2*pi,dcut2,beta2,btheta2,bphi,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==5)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
            bwrapstrapthetalattice(bxcom,bycom,Rb2,t,2*pi,dcut2,beta2,btheta2,bphi,rotangle,veccam1)
        elseif (bregime==6)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,btheta1,0,rotangle,veccam1)
            bwrapstrapthetalattice(bxcom,bycom,Rb2,t,2*pi,dcut2,beta2,btheta2,bphi,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==7)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bstraplattice(bvec1,bvec2-bvec1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,btheta1,0,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==8)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,btheta1,0,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==9)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bwrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,btheta1,0,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==10)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(battach1,battalpha1-battach1,t,t/4,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        elseif (bregime==11)
            bwrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            bwrapstrapconnectlattice(bxcom,bycom,Rb2,t,2*pi,beta2,C,bphi,rotangle,veccam1)
            bstraplattice(bLoc,bLsvec,t,t/2,rotangle,veccam1)
        end
    end
end

%%%% draw front layer wrapped straps
for i=1:jointnum
    phi1=LocMat(FullJointMat(i,1),1);
    phi2=LocMat(FullJointMat(i,2),1);
    veccam2=[LocMat(FullJointMat(i,2),2); LocMat(FullJointMat(i,2),3)];
    veccam1=[LocMat(FullJointMat(i,1),2); LocMat(FullJointMat(i,1),3)];
    difvec=veccam2-veccam1;
    rotangle=FullJointMat(i,3)-(pi/2)+phi1;
    rotvec=[cos(rotangle), -sin(rotangle); sin(rotangle), cos(rotangle)]\difvec;
    Rp1=FullJointMat(i,6);
    Rp2=FullJointMat(i,7);
    Rb1=Rp1-(t/2);
    Rb2=Rp2-(t/2);
    beta1=FullJointMat(i,4);
    beta2=FullJointMat(i,5);
    phi=phi2-phi1;
    xcom=rotvec(1,1);
    ycom=rotvec(2,1);
    
    % front layer
    [tension, regime, attach1, attalpha1, vec1, vec2, theta1, theta2, Loc, Lsvec, b, Fc]=kinematics(phi, xcom, ycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v);
    
    if (tension==1) % tension
        if (regime==1)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(xcom,ycom,Rb2,t,2*pi,dcut2,beta2,theta2,phi,rotangle,veccam1)
        elseif (regime==2)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
        elseif (regime==3)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
        elseif (regime==4)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(xcom,ycom,Rb2,t,2*pi,dcut2,beta2,theta2,phi,rotangle,veccam1)
        elseif (regime==5)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
        elseif (regime==6)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
        elseif (regime==7)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,theta1,0,rotangle,veccam1)
            wrapstrapthetalattice(xcom,ycom,Rb2,t,2*pi,dcut2,beta2,theta2,phi,rotangle,veccam1)
        elseif (regime==8)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,theta1,0,rotangle,veccam1)
        elseif (regime==9)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,theta1,0,rotangle,veccam1)
        end
    else % compression
        if (regime==1)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(xcom,ycom,Rb2,t,2*pi,dcut2,beta2,theta2,phi,rotangle,veccam1)
        elseif (regime==2)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
        elseif (regime==3)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
        elseif (regime==4)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(xcom,ycom,Rb2,t,2*pi,dcut2,beta2,theta2,phi,rotangle,veccam1)
        elseif (regime==5)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(xcom,ycom,Rb2,t,2*pi,dcut2,beta2,theta2,phi,rotangle,veccam1)
        elseif (regime==6)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,theta1,0,rotangle,veccam1)
            wrapstrapthetalattice(xcom,ycom,Rb2,t,2*pi,dcut2,beta2,theta2,phi,rotangle,veccam1)
        elseif (regime==7)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,theta1,0,rotangle,veccam1)
        elseif (regime==8)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,theta1,0,rotangle,veccam1)
        elseif (regime==9)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
            wrapstrapthetalattice(0,0,Rb1,t,pi,dcut1,beta1,theta1,0,rotangle,veccam1)
        elseif (regime==10)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
        elseif (regime==11)
            wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
            wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
        end
    end
end

%%%% draw circles
for i=1:jointnum
    phi1=LocMat(FullJointMat(i,1),1);
    phi2=LocMat(FullJointMat(i,2),1);
    veccam2=[LocMat(FullJointMat(i,2),2); LocMat(FullJointMat(i,2),3)];
    veccam1=[LocMat(FullJointMat(i,1),2); LocMat(FullJointMat(i,1),3)];
    difvec=veccam2-veccam1;
    rotangle=FullJointMat(i,3)-(pi/2)+phi1;
    rotvec=[cos(rotangle), -sin(rotangle); sin(rotangle), cos(rotangle)]\difvec;
    Rp1=FullJointMat(i,6);
    Rp2=FullJointMat(i,7);
    Rb1=Rp1-(t/2);
    Rb2=Rp2-(t/2);
    beta1=FullJointMat(i,4);
    beta2=FullJointMat(i,5);
    phi=phi2-phi1;
    xcom=rotvec(1,1);
    ycom=rotvec(2,1);
    
    % front layer
    [tension, regime, attach1, attalpha1, vec1, vec2, theta1, theta2, Loc, Lsvec, b, Fc]=kinematics(phi, xcom, ycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v);
    
    if (tension==1) % tension
        if (regime==1)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==2)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==3)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==4)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==5)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==6)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==7)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==8)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==9)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        end
    else % compression
        if (regime==1)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==2)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==3)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==4)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==5)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==6)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==7)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==8)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==9)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==10)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        elseif (regime==11)
            circlelattice(0,0,Rb1,boltr,0,rotangle,veccam1)
            circlelattice(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
        end
    end
end

%%%% draw front layer straps and fabcut
for i=1:jointnum
    phi1=LocMat(FullJointMat(i,1),1);
    phi2=LocMat(FullJointMat(i,2),1);
    veccam2=[LocMat(FullJointMat(i,2),2); LocMat(FullJointMat(i,2),3)];
    veccam1=[LocMat(FullJointMat(i,1),2); LocMat(FullJointMat(i,1),3)];
    difvec=veccam2-veccam1;
    rotangle=FullJointMat(i,3)-(pi/2)+phi1;
    rotvec=[cos(rotangle), -sin(rotangle); sin(rotangle), cos(rotangle)]\difvec;
    Rp1=FullJointMat(i,6);
    Rp2=FullJointMat(i,7);
    Rb1=Rp1-(t/2);
    Rb2=Rp2-(t/2);
    beta1=FullJointMat(i,4);
    beta2=FullJointMat(i,5);
    phi=phi2-phi1;
    xcom=rotvec(1,1);
    ycom=rotvec(2,1);
    
    % front layer
    [tension, regime, attach1, attalpha1, vec1, vec2, theta1, theta2, Loc, Lsvec, b, Fc]=kinematics(phi, xcom, ycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v);
    
    if (tension==1) % tension
        if (regime==1)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==2)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==3)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==4)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==5)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==6)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==7)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==8)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==9)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        end
    else % compression
        if (regime==1)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==2)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==3)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==4)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==5)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==6)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==7)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(vec1,vec2-vec1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==8)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==9)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==10)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(attach1,attalpha1-attach1,t,t/4,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        elseif (regime==11)
            fabcutlattice(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
            fabcutlattice(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
            straplattice(Loc,Lsvec,t,t/2,rotangle,veccam1)
        end
    end
end
function fabcutlattice(x,y,r,Z,dcut,beta,phi,rotangle,veccam1)
% draws fabricated cutout where x and y are center coordinates and r is
% the base circle radius. Z is the angle from
% the horiztontal to the normal joint direction. dcut is the
% fabricated cut by the strap that defines alpha. beta is define
% in the figure. phi is how much the cam rotates

inc=100;
alpha=acos((r-dcut)/r);
an1=Z-beta+phi;
an2=Z-beta-alpha+phi;

th = an1:-abs(an1-an2)/inc:an2;
xunit1 = r * cos(th) + x;
xunit = [xunit1, x+((r-dcut)*cos(Z-beta+phi))];
yunit1 = r * sin(th) + y;
yunit = [yunit1, y+((r-dcut)*sin(Z-beta+phi))];
zunit = zeros(size(yunit));

%rotate them back to the global coordinate system
xunitg=(xunit*cos(rotangle))-(yunit*sin(rotangle))+veccam1(1,1);
yunitg=(xunit*sin(rotangle))+(yunit*cos(rotangle))+veccam1(2,1);
xunit=xunitg;
yunit=yunitg;

patch(xunit, yunit, zunit, [.8 .8 .8], 'LineStyle', 'none');
function fabcutlattice2(x,y,r,Z,dcut,beta,phi,rotangle,veccam1)
% draws fabricated cutout where x and y are center coordinates and r is
% the base circle radius. Z is the angle from
% the horiztontal to the normal joint direction. dcut is the
% fabricated cut by the strap that defines alpha. beta is define
% in the figure. phi is how much the cam rotates

format long

hold on
inc=100;
alpha=acos((r-dcut)/r);
an1=Z-beta+phi;
an2=Z-beta-alpha+phi;

th = an1:-abs(an1-an2)/inc:an2;
xunit1 = r * cos(th) + x;
xunit = [xunit1, x+((r-dcut)*cos(Z-beta+phi))];
yunit1 = r * sin(th) + y;
yunit = [yunit1, y+((r-dcut)*sin(Z-beta+phi))];
zunit = zeros(size(yunit));

%rotate them back to the global coordinate system
xunitg=(xunit*cos(rotangle))-(yunit*sin(rotangle))+veccam1(1,1);
yunitg=(xunit*sin(rotangle))+(yunit*cos(rotangle))+veccam1(2,1);
xunit=xunitg;
yunit=yunitg;

patch(xunit, yunit, zunit, [1 1 1], 'LineStyle', 'none');

hold off
function [Wtot, tension, regime, bregime, bphi, bxcom, bycom, attach1, attalpha1, vec1, vec2, theta1, theta2, Loc, Lsvec, battach1, battalpha1, bvec1, bvec2, btheta1, btheta2, bLoc, bLsvec]=jointloads(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, O, Et, Ec, v, fric)
%This function calculates the load that needs to be imparted on cam 2 to
%move it to the positions xcom and ycom and the orientation phi.

bphi=-phi;
bxcom=-xcom;
bycom=ycom;

% front layer
[tension, regime, attach1, attalpha1, vec1, vec2, theta1, theta2, Loc, Lsvec, b, Fc]=kinematics(phi, xcom, ycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v);

% back layer
[btension, bregime, battach1, battalpha1, bvec1, bvec2, btheta1, btheta2, bLoc, bLsvec, bb, bFc]=kinematics(bphi, bxcom, bycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v);

if (tension==1) % tension scenario
    if (regime==1)
        [Wloadfront]=tensionregime1(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta2, vec1, vec2, Lsvec);
    elseif (regime==2)
        [Wloadfront]=tensionregime2(Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, vec1, vec2, Lsvec);
    elseif (regime==3)
        [Wloadfront]=tensionregime3(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, Loc, Lsvec);
    elseif (regime==4)
        [Wloadfront]=tensionregime4(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, attach1, attalpha1, Ec, Et, v, fric, theta2, vec1, vec2, Lsvec);
    elseif (regime==5)
        [Wloadfront]=tensionregime5(Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, attach1, attalpha1, Ec, Et, v, vec1, vec2, Lsvec);
    elseif (regime==6)
        [Wloadfront]=tensionregime6(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, attach1, attalpha1, Ec, Et, v, Loc, Lsvec);
    elseif (regime==7)
        [Wloadfront]=tensionregime7(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, theta2, vec1, vec2, Lsvec);
    elseif (regime==8)
        [Wloadfront]=tensionregime8(Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, vec1, vec2, Lsvec);
    elseif (regime==9)
        [Wloadfront]=tensionregime9(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, Loc, Lsvec);
    end
    if (bregime==1)
        [Wloadback]=tensionregime1(bphi, bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta2, bvec1, bvec2, bLsvec);
    elseif (bregime==2)
        [Wloadback]=tensionregime2(Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, bvec1, bvec2, bLsvec);
    elseif (bregime==3)
        [Wloadback]=tensionregime3(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, bLoc, bLsvec);
    elseif (bregime==4)
        [Wloadback]=tensionregime4(bphi, bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, battach1, battalpha1, Ec, Et, v, fric, btheta2, bvec1, bvec2, bLsvec);
    elseif (bregime==5)
        [Wloadback]=tensionregime5(Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, battach1, battalpha1, Ec, Et, v, bvec1, bvec2, bLsvec);
    elseif (bregime==6)
        [Wloadback]=tensionregime6(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, battach1, battalpha1, Ec, Et, v, bLoc, bLsvec);
    elseif (bregime==7)
        [Wloadback]=tensionregime7(bphi, bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta1, btheta2, bvec1, bvec2, bLsvec);
    elseif (bregime==8)
        [Wloadback]=tensionregime8(Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta1, bvec1, bvec2, bLsvec);
    elseif (bregime==9)
        [Wloadback]=tensionregime9(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta1, bLoc, bLsvec);
    end
    Wloadbackfinal=[-Wloadback(1,1); Wloadback(2,1); -Wloadback(3,1)];
    Wtot=(O/2)*(Wloadfront+Wloadbackfinal);
else % compression scenario
    if (regime==1)
        [Wloadfront, Ttanfront]=compressionregime1(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta2, vec1, vec2, Lsvec, Fc);
    elseif (regime==2)
        [Wloadfront, Ttanfront]=compressionregime2(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, vec1, vec2, Lsvec, Fc);
    elseif (regime==3)
        [Wloadfront, Ttanfront]=compressionregime3(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, Loc, Lsvec, Fc);
    elseif (regime==4)
        [Wloadfront, Ttanfront]=compressionregime4(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta2, vec1, vec2, Loc, Lsvec, Fc);
    elseif (regime==5)
        [Wloadfront, Ttanfront]=compressionregime5(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta2, attach1, vec1, vec2, Loc, Lsvec, Fc);
    elseif (regime==6)
        [Wloadfront, Ttanfront]=compressionregime6(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, theta2, vec1, vec2, Loc, Lsvec, Fc);
    elseif (regime==7)
        [Wloadfront, Ttanfront]=compressionregime7(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, vec1, vec2, Loc, Lsvec, Fc);
    elseif (regime==8)
        [Wloadfront, Ttanfront]=compressionregime8(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, vec1, Loc, Lsvec, Fc);
    elseif (regime==9)
        [Wloadfront, Ttanfront]=compressionregime9(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, Loc, Lsvec, Fc);
    elseif (regime==10)
        [Wloadfront, Ttanfront]=compressionregime10(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, attach1, attalpha1, Loc, Lsvec, Fc);
    elseif (regime==11)
        [Wloadfront, Ttanfront]=compressionregime11(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, Loc, Lsvec, Fc);
    end
    if (bregime==1)
        [Wloadback, Ttanb]=compressionregime1(bphi, bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta2, bvec1, bvec2, bLsvec, bFc);
    elseif (bregime==2)
        [Wloadback, Ttanb]=compressionregime2(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, bvec1, bvec2, bLsvec, bFc);
    elseif (bregime==3)
        [Wloadback, Ttanb]=compressionregime3(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, bLoc, bLsvec, bFc);
    elseif (bregime==4)
        [Wloadback, Ttanb]=compressionregime4(bphi, bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta2, bvec1, bvec2, bLoc, bLsvec, bFc);
    elseif (bregime==5)
        [Wloadback, Ttanb]=compressionregime5(bphi, bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta2, battach1, bvec1, bvec2, bLoc, bLsvec, bFc);
    elseif (bregime==6)
        [Wloadback, Ttanb]=compressionregime6(bphi, bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta1, btheta2, bvec1, bvec2, bLoc, bLsvec, bFc);
    elseif (bregime==7)
        [Wloadback, Ttanb]=compressionregime7(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta1, bvec1, bvec2, bLoc, bLsvec, bFc);
    elseif (bregime==8)
        [Wloadback, Ttanb]=compressionregime8(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta1, bvec1, bLoc, bLsvec, bFc);
    elseif (bregime==9)
        [Wloadback, Ttanb]=compressionregime9(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, btheta1, bLoc, bLsvec, bFc);
    elseif (bregime==10)
        [Wloadback, Ttanb]=compressionregime10(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, battach1, battalpha1, bLoc, bLsvec, bFc);
    elseif (bregime==11)
        [Wloadback, Ttanb]=compressionregime11(bxcom, bycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, bLoc, bLsvec, bFc);
    end
    Wloadbackfinal=[-Wloadback(1,1); Wloadback(2,1); -Wloadback(3,1)];
    % Determine if friction loads occur between the compressed cams
    unitcenter=[xcom, ycom]/sqrt(dotFast([xcom, ycom],[xcom, ycom]));
    tstar=t*(1-(Fc/(Ec*b*W)));
    fricloc=(sqrt((Rp1^2)-((b/2)^2))+(tstar/2))*unitcenter;
    dirvec=[fricloc(1,2), -fricloc(1,1)];
    unitdirvec=dirvec/sqrt(dot(dirvec,dirvec));
    Ttanback=-Ttanb;
    if(cutdec(Ttanfront,10)==0 && cutdec(Ttanback,10)==0)
        fricforce=[0,0];
    elseif(cutdec(Ttanfront,10)<0 && cutdec(Ttanback,10)<0)
        fricforce=-2*Fc*fric*unitdirvec;
    elseif(cutdec(Ttanfront,10)<0 && cutdec(Ttanback,10)==0)
        fricforce=-2*Fc*fric*unitdirvec;
    elseif(cutdec(Ttanfront,10)==0 && cutdec(Ttanback,10)<0)
        fricforce=-2*Fc*fric*unitdirvec;
    elseif(cutdec(Ttanfront,10)>0 && cutdec(Ttanback,10)>0)
        fricforce=2*Fc*fric*unitdirvec;
    elseif(cutdec(Ttanfront,10)>0 && cutdec(Ttanback,10)==0)
        fricforce=2*Fc*fric*unitdirvec;
    elseif(cutdec(Ttanfront,10)==0 && cutdec(Ttanback,10)>0)
        fricforce=2*Fc*fric*unitdirvec;
    else
        dirsign=Ttanfront+Ttanback;
        if(cutdec(dirsign,10)>0)
            fricforce=2*Fc*fric*unitdirvec;
        elseif(cutdec(dirsign,10)<0)
            fricforce=-2*Fc*fric*unitdirvec;
        else
            fricforce=[0,0];
        end
    end
    Wloadfric=[fricforce(1,1); fricforce(1,2); ((fricloc(1,1)*fricforce(1,2))-(fricloc(1,2)*fricforce(1,1)))];
    Wtot=(O/2)*(Wloadfront+Wloadbackfinal+Wloadfric);
end
function [tension, regime, attach1, attalpha1, vec1, vec2, theta1, theta2, Loc, Lsvec, b, Fc] = kinematics(phi, xcom, ycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v)
% This function determines the geometric compatiability of the strap within a single layer

Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
centers=[xcom, ycom];
D=sqrt(dotFast(centers,centers));
ancenters=acos(dotFast([1, 0],centers/D));
attach1=[(Rp1*cos(pi-beta1)), (Rp1*sin(pi-beta1))];
attalpha1=[(Rp1*cos(pi-beta1-alpha1)), (Rp1*sin(pi-beta1-alpha1))];
if (phi-beta2 < 0)
    anattach2=phi-beta2+(2*pi);
else
    anattach2=phi-beta2;
end
if (phi-beta2-alpha2 < 0)
    anattalpha2=phi-beta2-alpha2+(2*pi);
else
    anattalpha2=phi-beta2-alpha2;
end
attach2=[(Rp2*cos(anattach2)), (Rp2*sin(anattach2))];
attalpha2=[(Rp2*cos(anattalpha2)), (Rp2*sin(anattalpha2))];
vec1=centers+attach2;
magvec1=sqrt(dotFast(vec1,vec1));
angle1=acos(dotFast([1, 0],vec1/magvec1));
if(vec1(1,2)<0)
    angle1=((2*pi)-angle1);
    if(angle1>(3*pi/2))
        angle1=angle1-(2*pi);
    end
end
vec2=centers+attalpha2;
magvec2=sqrt(dotFast(vec2,vec2));
angle2=acos(dotFast([1, 0],vec2/magvec2));
if(vec2(1,2)<0)
    angle2=((2*pi)-angle2);
    if(angle2>(3*pi/2))
        angle2=angle2-(2*pi);
    end
end
slope1=attalpha1-attach1;
magslope1=sqrt(dotFast(slope1,slope1));
slopean=acos(dotFast([1, 0],slope1/magslope1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (cutdec(D-(Rp1+Rp2),10) >= 0)
    tension=1;
    b=0;
    Fc=0;
    
    attcen=centers-attach1;
    magattcen=sqrt(dotFast(attcen,attcen));
    halfan=asin(Rp2/magattcen);
    an=acos(dotFast([1, 0],attcen/magattcen));
    if((ycom-attach1(1,2))>=0)
        gamma=an-halfan;
    else
        gamma=((2*pi)-an-halfan);
        if(gamma>(3*pi/2))
            gamma=gamma-(2*pi);
        end
    end
    
    if(gamma>=slopean) % nothing on cam 1
        if((phi-beta2-alpha2)>(gamma-(pi/2))) % theta2 only
            regime=1;
            Ls=Rp2/tan(halfan);
            Lsvec=Ls*[cos(gamma), sin(gamma)];
            Loc=attach1;
            theta1=0;
            theta2=(phi-beta2-alpha2)-(gamma-(pi/2));
        elseif(((phi-beta2-alpha2)<=(gamma-(pi/2))) && ((phi-beta2-(alpha2/2))>(gamma-(pi/2)))) % no theta and stuck on corner of cam2
            regime=2;
            Lsvec=vec2-attach1;
            Loc=attach1;
            theta1=0;
            theta2=0;
        else % no theta and nothing on a corner
            regime=3;
            Lsvec=vec1-attach1;
            Loc=attach1;
            theta1=0;
            theta2=0;
        end
    else
        attcen=centers-attalpha1;
        magattcen=sqrt(dotFast(attcen,attcen));
        halfan=asin(Rp2/magattcen);
        an=acos(dotFast([1, 0],attcen/magattcen));
        if((ycom-attalpha1(1,2))>=0)
            gamma=an-halfan;
        else
            gamma=((2*pi)-an-halfan);
            if(gamma>(3*pi/2))
                gamma=gamma-(2*pi);
            end
        end
        if((gamma<slopean) && (gamma>=(slopean-(alpha2/2)))) % no theta 1 but stuck on corner of cam 1
            if((phi-beta2-alpha2)>(gamma-(pi/2))) % theta2 only
                regime=4;
                Ls=Rp2/tan(halfan);
                Lsvec=Ls*[cos(gamma), sin(gamma)];
                Loc=attalpha1;
                theta1=0;
                theta2=(phi-beta2-alpha2)-(gamma-(pi/2));
            elseif(((phi-beta2-alpha2)<=(gamma-(pi/2))) && ((phi-beta2-(alpha2/2))>(gamma-(pi/2)))) % no theta and stuck on corner of cam2
                regime=5;
                Lsvec=vec2-attalpha1;
                Loc=attalpha1;
                theta1=0;
                theta2=0;
            else
                Ltemp1=vec1-attalpha1;
                magLtemp1=sqrt(dotFast(Ltemp1,Ltemp1));
                anLtemp1=acos(dotFast([1, 0],Ltemp1/magLtemp1));
                Ltemp2=vec1-attach1;
                magLtemp2=sqrt(dotFast(Ltemp2,Ltemp2));
                anLtemp2=acos(dotFast([1, 0],Ltemp2/magLtemp2));
                if(anLtemp1<anLtemp2) % no theta but on corner of cam 1; cam 2 has a large negative phi
                    regime=6;
                    Lsvec=vec1-attalpha1;
                    Loc=attalpha1;
                    theta1=0;
                    theta2=0;
                else  % no theta and not on corner of cam 1; cam 2 has a large negative phi
                    regime=3;
                    Lsvec=vec1-attach1;
                    Loc=attach1;
                    theta1=0;
                    theta2=0;
                end
            end
        else % theta 1 exists
            h1=cutdec((D-(Rp1+Rp2)),10)/((Rp1+Rp2)/Rp1);
            h2=h1*(Rp2/Rp1);
            lam=acos(Rp1/(Rp1+h1));
            gamma=ancenters+lam-(pi/2);
            if((phi-beta2-alpha2)>(gamma-(pi/2))) % theta 1 and theta 2
                regime=7;
                Ls=(Rp1+Rp2)*tan(lam);
                Lsvec=Ls*[cos(gamma), sin(gamma)];
                Loc=Rp1*[cos(ancenters+lam), sin(ancenters+lam)];
                theta1=pi-ancenters-lam-alpha1-beta1;
                theta2=(phi-beta2-alpha2)-(gamma-(pi/2));
            else
                Ls8=sqrt((magvec2^2)-(Rp1^2));
                R1an8=(acos(Rp1/magvec2))+angle2;
                gamma8=R1an8-(pi/2);
                Ls9=sqrt((magvec1^2)-(Rp1^2));
                R1an9=(acos(Rp1/magvec1))+angle1;
                gamma9=R1an9-(pi/2);
                if((R1an8<(pi-beta1-alpha1)) && (R1an9<(pi-beta1-alpha1)))
                    if(gamma9>gamma8) % theta 1 but no theta 2 and on corner of cam 2
                        regime=8;
                        Lsvec=Ls8*[cos(gamma8),sin(gamma8)];
                        Loc=Rp1*[cos(R1an8), sin(R1an8)];
                        theta1=pi-alpha1-beta1-R1an8;
                        theta2=0;
                    else % theta 1 but no theta 2 and nothing on corner of cam 2
                        regime=9;
                        Lsvec=Ls9*[cos(gamma9),sin(gamma9)];
                        Loc=Rp1*[cos(R1an9), sin(R1an9)];
                        theta1=pi-alpha1-beta1-R1an9;
                        theta2=0;
                    end
                elseif(R1an8<(pi-beta1-alpha1))
                    regime=8;
                    Lsvec=Ls8*[cos(gamma8),sin(gamma8)];
                    Loc=Rp1*[cos(R1an8), sin(R1an8)];
                    theta1=pi-alpha1-beta1-R1an8;
                    theta2=0;
                elseif(R1an9<(pi-beta1-alpha1))
                    regime=9;
                    Lsvec=Ls9*[cos(gamma9),sin(gamma9)];
                    Loc=Rp1*[cos(R1an9), sin(R1an9)];
                    theta1=pi-alpha1-beta1-R1an9;
                    theta2=0;
                else
                    v11=vec1-attach1;
                    magv11=sqrt(dotFast(v11,v11));
                    anv11=acos(dotFast([1, 0],v11/magv11));
                    v12=vec1-attalpha1;
                    magv12=sqrt(dotFast(v12,v12));
                    anv12=acos(dotFast([1, 0],v12/magv12));
                    v21=vec2-attach1;
                    magv21=sqrt(dotFast(v21,v21));
                    anv21=acos(dotFast([1, 0],v21/magv21));
                    v22=vec2-attalpha1;
                    magv22=sqrt(dotFast(v22,v22));
                    anv22=acos(dotFast([1, 0],v22/magv22));
                    vangles=[anv11, anv12, anv21, anv22];
                    smallest=min(vangles);
                    if(anv11==smallest)
                        regime=3;
                        Lsvec=v11;
                        Loc=attach1;
                        theta1=0;
                        theta2=0;
                    elseif(anv12==smallest)
                        regime=6;
                        Lsvec=v12;
                        Loc=attalpha1;
                        theta1=0;
                        theta2=0;
                    elseif(anv21==smallest)
                        regime=2;
                        Lsvec=v21;
                        Loc=attach1;
                        theta1=0;
                        theta2=0;
                    elseif(anv22==smallest)
                        regime=5;
                        Lsvec=v22;
                        Loc=attalpha1;
                        theta1=0;
                        theta2=0;
                    end
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Compression
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    tension=0;
    
    Compress=abs(cutdec(D-(Rp1+Rp2),10));
    inc=dcut1/(1e5);
    Fc=inc;
    error=1e10;
    ct=0;
    while error > 1e-10
        ct=ct+1;
        b=4*sqrt((2*Fc*((1-(v^2))/Ec))/(pi*W*((1/Rb1)+(1/Rb2))));
        func=(((2*Fc)/(pi*W))*((1-(v^2))/Ec)*((2/3)+(log(8*Rb1/b))+(log(8*Rb2/b))))+((Fc*t)/(Ec*W*b))-Compress;
        binc=4*sqrt((2*(Fc+inc)*((1-(v^2))/Ec))/(pi*W*((1/Rb1)+(1/Rb2))));
        funcinc=(((2*(Fc+inc))/(pi*W))*((1-(v^2))/Ec)*((2/3)+(log(8*Rb1/binc))+(log(8*Rb2/binc))))+(((Fc+inc)*t)/(Ec*W*b))-Compress;
        deriv=(funcinc-func)/inc;
        newFc=Fc-(func/deriv);
        error=abs(Fc-newFc);
        Fc=newFc;
        if(ct>1e4)
            break
        end
    end
    halfanb=asin((b/2)/Rp1);
    angleb=ancenters-halfanb;
    vecb=Rp1*[cos(angleb), sin(angleb)];
    anglebback=ancenters+halfanb;
    vecbback=Rp1*[cos(anglebback), sin(anglebback)];
    
    if(angleb>=pi-beta1)
        aa=sqrt((D^2)+(Rp1^2)-(2*D*Rp1*cos(ancenters-(pi-beta1))));
        Ls=sqrt((aa^2)-(Rp2^2));
        curl=asin(Ls/aa);
        Zan=asin((Rp1*sin(ancenters-(pi-beta1)))/aa);
        Len=sqrt((D^2)+(Rp2^2)-(2*D*Rp2*cos(curl+Zan)));
        trian=((pi/2)-curl)+(pi-Zan-(ancenters-(pi-beta1)));
        if (trian<pi)
            sigan=asin((Ls*sin(trian))/Len);
            vecperp=Len*[cos(pi-beta1+sigan), sin(pi-beta1+sigan)];
        elseif (trian>pi)
            trian=(2*pi)-trian;
            sigan=asin((Ls*sin(trian))/Len);
            vecperp=Len*[cos(pi-beta1-sigan), sin(pi-beta1-sigan)];
        else
            vecperp=Len*[cos(pi-beta1), sin(pi-beta1)];
        end
        R2vecs=vecperp-centers;
        magR2vecs=sqrt(dotFast(R2vecs,R2vecs));
        anR2vecs=acos(dotFast([1, 0],R2vecs/magR2vecs));
        if(R2vecs(1,2)<0)
            anR2vecs=-anR2vecs;
        end
        slopevec1=vec1-attach1;
        magslopevec1=sqrt(dotFast(slopevec1,slopevec1));
        anslopevec1=acos(dotFast([1, 0],slopevec1/magslopevec1));
        if(slopevec1(1,2)<0)
            anslopevec1=(2*pi)-anslopevec1;
        end
        slopevec2=vec2-attach1;
        magslopevec2=sqrt(dotFast(slopevec2,slopevec2));
        anslopevec2=acos(dotFast([1, 0],slopevec2/magslopevec2));
        if(slopevec2(1,2)<0)
            anslopevec2=(2*pi)-anslopevec2;
        end
        if((phi-beta2-alpha2)>anR2vecs) % No theta 1 but theta 2. Strap stretched out.
            regime=1;
            Lsvec=vecperp-attach1;
            Loc=attach1;
            theta1=0;
            theta2=phi-beta2-alpha2-anR2vecs;
        elseif(anslopevec2<anslopevec1) %((phi-beta2-alpha2)<=anR2vecs) && ((phi-beta2-(alpha2/2))>anR2vecs)) % No theta 1 and no theta 2. Strap stretch and stuck on corner of cam 2
            regime=2;
            Lsvec=vec2-attach1;
            Loc=attach1;
            theta1=0;
            theta2=0;
        else % No theta 1 and no theta 2. Strap stretched but not stuck on any corners
            regime=3;
            Lsvec=vec1-attach1;
            Loc=attach1;
            theta1=0;
            theta2=0;
        end
    else
        R2vecs=vecb-centers;
        magR2vecs=sqrt(dotFast(R2vecs,R2vecs));
        anR2vecs=acos(dotFast([1, 0],R2vecs/magR2vecs));
        if(R2vecs(1,2)<0)
            anR2vecs=-anR2vecs;
        end
        if((angleb<pi-beta1) && (angleb>=pi-beta1-alpha1) && (anglebback>=pi-beta1)) % Front of crunch region after attach1 but before attalpha1 and back of region is before attach1
            regime=4; % only strap is in cruch region and theta2 only
            unitattach1=attach1/sqrt(dotFast(attach1,attach1));
            newlen=((b/2)/tan(halfanb))/cos(pi-beta1-ancenters);
            veccrunch=newlen*unitattach1;
            Lsvec=vecb-veccrunch;
            Loc=veccrunch;
            theta1=0;
            theta2=phi-beta2-alpha2-anR2vecs;
        elseif((angleb<pi-beta1) && (angleb>=pi-beta1-alpha1) && (anglebback<pi-beta1)) % Entire crunch region is between attach1 and attalpha1
            regime=5;
            Lsvec=vecb-vecbback;
            Loc=vecbback;
            attalpha1=vecbback; % Strap to the left starts at attach1 and goes to vecbback
            theta1=0;
            theta2=phi-beta2-alpha2-anR2vecs;
        elseif((angleb<pi-beta1-alpha1) && (anglebback>=pi-beta1)) % Crunch region after attalpha1 but back side is behind attach1
            regime=4; % only strap is in cruch region and theta2 only
            unitattach1=attach1/sqrt(dotFast(attach1,attach1));
            newlen=((b/2)/tan(halfanb))/cos(pi-beta1-ancenters);
            veccrunch=newlen*unitattach1;
            Lsvec=vecb-veccrunch;
            Loc=veccrunch;
            theta1=0;
            theta2=phi-beta2-alpha2-anR2vecs;
        elseif((angleb<pi-beta1-alpha1) && (anglebback<pi-beta1) && (anglebback>=pi-beta1-alpha1)) % Crunch region after attalpha1 and back side is between attach1 and attalpha1
            regime=5; % assume that cam 2 is wrapped with theta 2
            Lsvec=vecb-vecbback;
            Loc=vecbback;
            attalpha1=vecbback; % Strap to the left starts at attach1 and goes to vecbback
            theta1=0;
            theta2=phi-beta2-alpha2-anR2vecs;
        else % Front and back of crunch region is after attalpha1
            if((phi-beta2-alpha2)>(anR2vecs)) % Theta 2 and theta 1
                regime=6;
                Lsvec=vecb-vecbback; % the strap before connects from attach 1 to attalpha1 as usual
                Loc=vecbback;
                theta1=pi-beta1-alpha1-anglebback;
                theta2=phi-beta2-alpha2-anR2vecs;
            elseif(((phi-beta2-alpha2)<=(anR2vecs)) && ((phi-beta2)>(anR2vecs))) % No theta 2 but theta 1. Straight strap on right side
                regime=7;
                Lsvec=vecb-vecbback;
                Loc=vecbback;
                vec2=vecb;  % the strap after (right side) connects from vec1 to vecb
                theta1=pi-beta1-alpha1-anglebback;
                theta2=0;
            elseif(((phi-beta2)<=(anR2vecs)) && ((phi-beta2)>(anR2vecs-(2*asin((b/2)/Rp2))))) % No theta 2 but theta 1. Straight strap in crunched region from vecbback to where it attaches to cam 2
                regime=8;
                unitattach2=attach2/sqrt(dotFast(attach2,attach2));
                newlen=(D-((b/2)/tan(halfanb)))/cos(pi-ancenters+phi-beta2);
                veccrunch=centers+(newlen*unitattach2);
                Lsvec=veccrunch-vecbback; % No strap to the right in this regime
                Loc=vecbback;
                theta1=pi-beta1-alpha1-anglebback;
                theta2=0;
            elseif((phi-beta2)<=(anR2vecs-(2*asin((b/2)/Rp2))))
                vecstrapalpha=vec1-attalpha1;
                Lsalpha1=sqrt(dotFast(vecstrapalpha,vecstrapalpha));
                gammaalpha1=acos(dotFast([1, 0],vecstrapalpha/Lsalpha1));
                vecstrapbeta=vec1-attach1;
                Lsbeta1=sqrt(dotFast(vecstrapbeta,vecstrapbeta));
                gammabeta1=acos(dotFast([1, 0],vecstrapbeta/Lsbeta1));
                R1an=(acos(Rp1/magvec1))+angle1;
                if(R1an<(pi-beta1-alpha1)) % Strap lifts off cam 1 with theta1 and connects at vec1 on cam 2 rotated a large netative
                    regime=9;
                    Ls=sqrt((magvec1^2)-(Rp1^2));
                    gammatheta1=R1an-(pi/2);
                    Lsvec=Ls*[cos(gammatheta1), sin(gammatheta1)];
                    Loc=Rp1*[cos(R1an), sin(R1an)];
                    theta1=pi-beta1-alpha1-R1an;
                    theta2=0;
                elseif(gammabeta1<slopean) % Strap lifts off of corner of cam 1 with no theta1 and connects at vec1 on cam 2 rotated a large netative
                    regime=10;
                    Lsvec=vec1-attalpha1;
                    Loc=attalpha1;
                    theta1=0;
                    theta2=0;
                elseif(gammabeta1>=slopean) % Strap connects attach1 to vec1 on cam 2 rotated a large netative with no thetas
                    regime=11;
                    Lsvec=vec1-attach1;
                    Loc=attach1;
                    theta1=0;
                    theta2=0;
                end
            end
        end
    end
end
function [Wlatticetot]=latticeforce(FullJointMat, LocMat, jointnum, rm, cg, t, W, del, dcut1, dcut2, O, Et, Ec, v, fric)
% This function returns a vector containing all the 2D wrench vectors that
% must be imparted on all the non-grounded cams within a lattice so that
% they can be located and oriented in the positions contained in LocMat.

% This portion sums up all the loads on each cam imposed by their
% neighboring cams to which they are joined by straps
Wlattice=zeros(3*(rm-cg),1);
for j=1:(rm-cg)
    Wglobal=[0; 0; 0];
    for k=1:(2*jointnum)
        if(FullJointMat(k,2)==j)
            phi1=LocMat(FullJointMat(k,1),1);
            phi2=LocMat(FullJointMat(k,2),1);
            veccam1=[LocMat(FullJointMat(k,1),2); LocMat(FullJointMat(k,1),3)];
            veccam2=[LocMat(FullJointMat(k,2),2); LocMat(FullJointMat(k,2),3)];
            difvec=veccam2-veccam1;
            rotangle=FullJointMat(k,3)-(pi/2)+phi1;
            RotMat=[cos(rotangle), -sin(rotangle); sin(rotangle), cos(rotangle)];
            rotvec=RotMat\difvec;
            Rp1=FullJointMat(k,6);
            Rp2=FullJointMat(k,7);
            Rb1=Rp1-(t/2);
            Rb2=Rp2-(t/2);
            beta1=FullJointMat(k,4);
            beta2=FullJointMat(k,5);
            phi=phi2-phi1;
            xcom=rotvec(1,1);
            ycom=rotvec(2,1);
            [Wtot]=jointloads(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, O, Et, Ec, v, fric);
            
            n1x=cos(rotangle);
            n1y=sin(rotangle);
            n2x=cos(rotangle+(pi/2));
            n2y=sin(rotangle+(pi/2));
            transwrench=[n1x, n2x, 0; n1y, n2y, 0; ((veccam1(1,1)*n1y)-(veccam1(2,1)*n1x)), ((veccam1(1,1)*n2y)-(veccam1(2,1)*n2x)), 1];
            Wglobal=Wglobal+(transwrench*Wtot);
        end
    end
    Wlattice((3*(j-1))+1:(3*j),1)=Wglobal;
end

%%% This portion sums up additional collision forces on cams that bump into
%%% other cams to which they are not joined by straps. these forces prevent
%%% the cams from passing through each other as their lattice is deformed
% Woverlap=zeros(3*rm,1);
% [YoN, culpritvec]=overlap(LocMat);
% history=[0];
% for i=1:rm
%     for j=1:rm
%         gonogo=1;
%         for k=1:(2*jointnum)
%             if(FullJointMat(k,1)==i && FullJointMat(k,2)==j)
%                 gonogo=0;
%             end
%         end
%         if(gonogo==1 && culpritvec(i,j)~=0)
%             Compress=culpritvec(i,j);
%             inc=dcut1/(1e5);
%             Fc=inc;
%             error=1e10;
%             Rp1=LocMat(i,4);
%             Rp2=LocMat(j,4);
%             Rb1=Rp1-(t/2);
%             Rb2=Rp2-(t/2);
%             while (error > 1e-10)
%                 b=4*sqrt((2*Fc*((1-(v^2))/Ec))/(pi*W*((1/Rb1)+(1/Rb2))));
%                 func=(((2*Fc)/(pi*W))*((1-(v^2))/Ec)*((2/3)+(log(8*Rb1/b))+(log(8*Rb2/b))))+((Fc*t)/(Ec*W*b))-Compress;
%                 binc=4*sqrt((2*(Fc+inc)*((1-(v^2))/Ec))/(pi*W*((1/Rb1)+(1/Rb2))));
%                 funcinc=(((2*(Fc+inc))/(pi*W))*((1-(v^2))/Ec)*((2/3)+(log(8*Rb1/binc))+(log(8*Rb2/binc))))+(((Fc+inc)*t)/(Ec*W*b))-Compress;
%                 deriv=cutdec((funcinc-func),10)/inc;
%                 newFc=Fc-(func/deriv);
%                 error=abs(Fc-newFc);
%                 Fc=newFc;
%             end
%             veccami=[LocMat(i,2); LocMat(i,3)];
%             veccamj=[LocMat(j,2); LocMat(j,3)];
%             vecij=veccamj-veccami;
%             vecji=veccami-veccamj;
%             unitvecij=vecij/sqrt(dot(vecij,vecij));
%             unitvecji=vecji/sqrt(dot(vecji,vecji));
%             Fij=Fc*unitvecij;
%             Fji=Fc*unitvecji;
%             Woverlap((3*(i-1))+1:(3*i),1)=Woverlap((3*(i-1))+1:(3*i),1)+[Fji(1,1); Fji(2,1); (veccami(1,1)*Fji(2,1))-(veccami(2,1)*Fji(1,1))];
%             Woverlap((3*(j-1))+1:(3*j),1)=Woverlap((3*(j-1))+1:(3*j),1)+[Fij(1,1); Fij(2,1); (veccamj(1,1)*Fij(2,1))-(veccamj(2,1)*Fij(1,1))];
%         end
%     end
% end
% Wcontact=Woverlap(1:(3*(rm-cg)),1);

Wlatticetot=Wlattice; %+Wcontact;
function [YoN, culpritvec]=overlap(LocMat)
% This function determines if any circular cams in the lattice overlap and
% returns a matrix that stores information pertaining to how much they overlap

[rm, cm]= size(LocMat);
culpritvec=zeros(rm, rm);
dismat=zeros(rm, rm);
YoN=0;
for i=1:rm
    for j=1:rm
        pos1=[LocMat(i,2), LocMat(i,3)];
        pos2=[LocMat(j,2), LocMat(j,3)];
        dismat(i,j)=sqrt(dot((pos1-pos2),(pos1-pos2)));
        if((i~=j) && cutdec((dismat(i,j)-LocMat(i,4)-LocMat(j,4)),10)<0)
            YoN=1;
            culpritvec(i,j)=cutdec(abs(dismat(i,j)-LocMat(i,4)-LocMat(j,4)),10);
        end
    end
end
function h = plaincircle(x,y,r)
% draws circle where x and y are center coordinates and r is the radius
inc=250;
th = 0:pi/inc:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
% zunit = zeros(size(yunit));
h = plot(xunit, yunit, 'k');
function straplattice(Loc,Lsvec,t,extra,rotangle,veccam1)
% This funtion draws a rectangular strap of thickness t that begins at Loc
% and points in the direction and length of Lsvec. It also adds an extra
% length to both sides of length extra.

Ls=sqrt(dot(Lsvec,Lsvec));
unitLs=Lsvec/Ls;
perpLs=[unitLs(1,2), -unitLs(1,1)];
strapx=[(Loc(1,1)+((t/2)*perpLs(1,1))-(extra*unitLs(1,1))), (Loc(1,1)-((t/2)*perpLs(1,1))-(extra*unitLs(1,1))), (Loc(1,1)-((t/2)*perpLs(1,1))+((Ls+extra)*unitLs(1,1))), (Loc(1,1)+((t/2)*perpLs(1,1))+((Ls+extra)*unitLs(1,1)))];
strapy=[(Loc(1,2)+((t/2)*perpLs(1,2))-(extra*unitLs(1,2))), (Loc(1,2)-((t/2)*perpLs(1,2))-(extra*unitLs(1,2))), (Loc(1,2)-((t/2)*perpLs(1,2))+((Ls+extra)*unitLs(1,2))), (Loc(1,2)+((t/2)*perpLs(1,2))+((Ls+extra)*unitLs(1,2)))];
strapz = zeros(size(strapy));

%rotate them back to the global coordinate system
xunitg=(strapx*cos(rotangle))-(strapy*sin(rotangle))+veccam1(1,1);
yunitg=(strapx*sin(rotangle))+(strapy*cos(rotangle))+veccam1(2,1);
strapx=xunitg;
strapy=yunitg;

patch(strapx, strapy, strapz, [0.6 0.6 0.6], 'LineStyle', 'none');
function [Wload]=tensionregime1(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta2, vec1, vec2, Lsvec)
% This function calculates the loads on a layer in the tension scenario
% within regime 1

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp2*sin(alpha2/2))-(Rp2*theta2);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=(2*Rp2*sin(alpha2/2))+(Rp2*((exp(fric*theta2)-1)/fric))+(Li*exp(fric*theta2));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
R2vec=vec2-[xcom, ycom];
perpvec=-[R2vec(1,2), -R2vec(1,1)];
unitperp=(perpvec)/sqrt(dotFast((perpvec),(perpvec)));
forcecorner2=Talpha2*unitperp;
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

% sum of reaction normal forces from the strap on the cam
dd=phi-beta2-alpha2;
sd=sin(dd);
cd=cos(dd);
sdt=sin(dd-real(theta2));
cdt=cos(dd-real(theta2));
Ncon=Talpha2*exp(fric*dd)*(2/(1+(fric^2)));
NFx=Ncon*(((sd/(2*exp(fric*dd)))-(fric*cd/(2*exp(fric*dd))))-((sdt/(2*exp(fric*(dd-real(theta2)))))-(fric*cdt/(2*exp(fric*(dd-real(theta2)))))));
NFy=Ncon*(((-fric*sd/(2*exp(fric*dd)))-(cd/(2*exp(fric*dd))))-((-fric*sdt/(2*exp(fric*(dd-real(theta2)))))-(cdt/(2*exp(fric*(dd-real(theta2)))))));
NMz=(xcom*NFy)-(ycom*NFx);
Wforcenormal=[NFx; NFy; NMz];

% sum of reaction tangential friction forces from the strap on the cam
fFx=-fric*NFy;
fFy=fric*NFx;
fMz=(fFy*xcom)-(fFx*ycom)+(Rb2*Talpha2*(exp(fric*real(theta2))-1));
Wforcetan=[fFx; fFy; fMz];

Wforce=Wforceattach+Wforcecorner+Wforcenormal+Wforcetan;

% % moment to resist the strap bent around the cams
% Momenttheta2=(E*W*(t^3))/(12*Rp2);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     rotbeta1=(2*pi)-acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% else
%     rotbeta1=acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% end
% if(Talpha2>0)
%     Momentbeta1=LivingHingeK*(rotbeta1);
%     Momentalpha2=LivingHingeK*(alpha2/2);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentbeta1=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta2+Momentbeta1+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload]=tensionregime2(Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, vec1, vec2, Lsvec)
% This function calculates the loads on a layer in the tension scenario
% within regime 2

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls+(2*Rp2*sin(alpha2/2))-Lo,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=Lo;
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
if(Ls==0)
    perpLoc=[Locforce2(1,2), -Locforce2(1,1)];
    forcecorner2=Talpha2*(perpLoc/sqrt(dotFast(perpLoc,perpLoc)));
else
    forcecorner2=Talpha2*(Lsvec/Ls);
end
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

Wforce=Wforceattach+Wforcecorner;

% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     rotbeta1=(2*pi)-acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% else
%     rotbeta1=acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% end
% slopevec=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
% anbetween=acos(dotFast(slopevec,Lsvec/Ls));
% if(Talpha2>0)
%     Momentbeta1=LivingHingeK*(rotbeta1);
%     Momentalpha2=LivingHingeK*(anbetween);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentbeta1=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momentbeta1+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload]=tensionregime3(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, Loc, Lsvec)
% This function calculates the loads on a layer in the tension scenario
% within regime 3

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Lo,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=Lo;
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% force to resist tension in the strap
Locforce=Loc+Lsvec;
if(Ls==0)
    perpLoc=[Locforce(1,2), -Locforce(1,1)];
    force=Talpha2*(perpLoc/sqrt(dotFast(perpLoc,perpLoc)));
else
    force=Talpha2*(Lsvec/Ls);
end
Wforce=[force(1,1); force(1,2); ((Locforce(1,1)*force(1,2))-(Locforce(1,2)*force(1,1)))]; % tension force from strap

% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Lsvec(1,2)<0)
%     rotbeta1=(2*pi)-acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% else
%     rotbeta1=acos(dotFast([1, 0],Lsvec/Ls))-((pi/2)-beta1);
% end
% vec1=Loc+Lsvec;
% R2vec=vec1-[xcom, ycom];
% perpvec=[R2vec(1,2), -R2vec(1,1)];
% unitperp=perpvec/sqrt(dotFast(perpvec,perpvec));
% if(unitperp(1,2)>0)
%     anperp=(2*pi)-acos(dotFast([1,0],unitperp));
% else
%     anperp=acos(dotFast([1,0],unitperp));
% end
% unitLs=-Lsvec/Ls;
% if(unitLs(1,2)>0)
%     anLs=(2*pi)-acos(dotFast([1,0],unitLs));
% else
%     anLs=acos(dotFast([1,0],unitLs));
% end
% anbetween=anLs-anperp;
%
% if(Talpha2>0)
%     Momentbeta1=LivingHingeK*(rotbeta1);
%     Momentbeta2=LivingHingeK*(anbetween);
% else
%     Momentbeta1=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momentbeta1+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload]=tensionregime4(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, attach1, attalpha1, Ec, Et, v, fric, theta2, vec1, vec2, Lsvec)
% This function calculates the loads on a layer in the tension scenario
% within regime 4

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp2*sin(alpha2/2))-(Rp2*theta2);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls+(2*Rp1*sin(alpha1/2))-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=(2*Rp2*sin(alpha2/2))+(Rp2*((exp(fric*theta2)-1)/fric))+(Li*exp(fric*theta2));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
R2vec=vec2-[xcom, ycom];
perpvec=-[R2vec(1,2), -R2vec(1,1)];
unitperp=(perpvec)/sqrt(dotFast((perpvec),(perpvec)));
forcecorner2=Talpha2*unitperp;
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

% sum of reaction normal forces from the strap on the cam
dd=phi-beta2-alpha2;
sd=sin(dd);
cd=cos(dd);
sdt=sin(dd-real(theta2));
cdt=cos(dd-real(theta2));
Ncon=Talpha2*exp(fric*dd)*(2/(1+(fric^2)));
NFx=Ncon*(((sd/(2*exp(fric*dd)))-(fric*cd/(2*exp(fric*dd))))-((sdt/(2*exp(fric*(dd-real(theta2)))))-(fric*cdt/(2*exp(fric*(dd-real(theta2)))))));
NFy=Ncon*(((-fric*sd/(2*exp(fric*dd)))-(cd/(2*exp(fric*dd))))-((-fric*sdt/(2*exp(fric*(dd-real(theta2)))))-(cdt/(2*exp(fric*(dd-real(theta2)))))));
NMz=(xcom*NFy)-(ycom*NFx);
Wforcenormal=[NFx; NFy; NMz];

% sum of reaction tangential friction forces from the strap on the cam
fFx=-fric*NFy;
fFy=fric*NFx;
fMz=(fFy*xcom)-(fFx*ycom)+(Rb2*Talpha2*(exp(fric*real(theta2))-1));
Wforcetan=[fFx; fFy; fMz];

Wforce=Wforceattach+Wforcecorner+Wforcenormal+Wforcetan;

% % moment to resist the strap bent around the cams
% Momenttheta2=(E*W*(t^3))/(12*Rp2);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% slopevec=(attalpha1-attach1)/sqrt(dotFast((attalpha1-attach1),(attalpha1-attach1)));
% rotalpha1=acos(dotFast([1, 0],slopevec))-acos(dotFast([1, 0],Lsvec/Ls));
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*-rotalpha1;
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentalpha2=LivingHingeK*(alpha2/2);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta2+Momentalpha1+Momentbeta1+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload]=tensionregime5(Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, attach1, attalpha1, Ec, Et, v, vec1, vec2, Lsvec)
% This function calculates the loads on a layer in the tension scenario
% within regime 5

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls+(2*Rp2*sin(alpha2/2))+(2*Rp1*sin(alpha1/2))-Lo,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=Lo;
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
if(Ls==0)
    perpLoc=[Locforce2(1,2), -Locforce2(1,1)];
    forcecorner2=Talpha2*(perpLoc/sqrt(dotFast(perpLoc,perpLoc)));
else
    forcecorner2=Talpha2*(Lsvec/Ls);
end
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

Wforce=Wforceattach+Wforcecorner;

% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% slopevec1=(attalpha1-attach1)/sqrt(dotFast((attalpha1-attach1),(attalpha1-attach1)));
% rotalpha1=acos(dotFast([1, 0],slopevec1))-acos(dotFast([1, 0],Lsvec/Ls));
% slopevec2=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
% anbetween=acos(dotFast(slopevec2,Lsvec/Ls));
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-rotalpha1);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentalpha2=LivingHingeK*(anbetween);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momentalpha1+Momentbeta1+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload]=tensionregime6(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, attach1, attalpha1, Ec, Et, v, Loc, Lsvec)
% This function calculates the loads on a layer in the tension scenario
% within regime 6

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls+(2*Rp1*sin(alpha1/2))-Lo,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=Lo;
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% force to resist tension in the strap
Locforce=Loc+Lsvec;
if(Ls==0)
    perpLoc=[Locforce(1,2), -Locforce(1,1)];
    force=Talpha2*(perpLoc/sqrt(dotFast(perpLoc,perpLoc)));
else
    force=Talpha2*(Lsvec/Ls);
end
Wforce=[force(1,1); force(1,2); ((Locforce(1,1)*force(1,2))-(Locforce(1,2)*force(1,1)))]; % tension force from strap

% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% slopevec1=(attalpha1-attach1)/sqrt(dotFast((attalpha1-attach1),(attalpha1-attach1)));
% rotalpha1=acos(dotFast([1, 0],slopevec1))-acos(dotFast([1, 0],Lsvec/Ls));
% vec1=Loc+Lsvec;
% R2vec=vec1-[xcom, ycom];
% perpvec=[R2vec(1,2), -R2vec(1,1)];
% unitperp=perpvec/sqrt(dotFast(perpvec,perpvec));
% if(unitperp(1,2)>0)
%     anperp=(2*pi)-acos(dotFast([1,0],unitperp));
% else
%     anperp=acos(dotFast([1,0],unitperp));
% end
% unitLs=-Lsvec/Ls;
% if(unitLs(1,2)>0)
%     anLs=(2*pi)-acos(dotFast([1,0],unitLs));
% else
%     anLs=acos(dotFast([1,0],unitLs));
% end
% anbetween=anLs-anperp;
%
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-rotalpha1);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentalpha2=LivingHingeK*(anbetween);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentalpha2=0;
% end
%
% WMoment=[0; 0; (Momentalpha1+Momentbeta1+Momentalpha2)];

Wload=Wforce; %+WMoment;
function [Wload]=tensionregime7(phi, xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, theta2, vec1, vec2, Lsvec)
% This function calculates the loads on a layer in the tension scenario
% within regime 7

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp1*sin(alpha1/2))-(2*Rp2*sin(alpha2/2))-(Rp1*theta1)-(Rp2*theta2);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch);
    denominat=(2*Rp2*sin(alpha2/2))+(2*Rp1*sin(alpha1/2)*exp(fric*(theta2-theta1)))+(Rp2*((exp(fric*theta2)-1)/fric))+(Rp1*((exp(fric*theta1)-1)/fric)*exp(fric*(theta2-theta1)))+(Li*exp(fric*theta2));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
R2vec=vec2-[xcom, ycom];
perpvec=-[R2vec(1,2), -R2vec(1,1)];
unitperp=(perpvec)/sqrt(dotFast((perpvec),(perpvec)));
forcecorner2=Talpha2*unitperp;
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

% sum of reaction normal forces from the strap on the cam
dd=phi-beta2-alpha2;
sd=sin(dd);
cd=cos(dd);
sdt=sin(dd-real(theta2));
cdt=cos(dd-real(theta2));
Ncon=Talpha2*exp(fric*dd)*(2/(1+(fric^2)));
NFx=Ncon*(((sd/(2*exp(fric*dd)))-(fric*cd/(2*exp(fric*dd))))-((sdt/(2*exp(fric*(dd-real(theta2)))))-(fric*cdt/(2*exp(fric*(dd-real(theta2)))))));
NFy=Ncon*(((-fric*sd/(2*exp(fric*dd)))-(cd/(2*exp(fric*dd))))-((-fric*sdt/(2*exp(fric*(dd-real(theta2)))))-(cdt/(2*exp(fric*(dd-real(theta2)))))));
NMz=(xcom*NFy)-(ycom*NFx);
Wforcenormal=[NFx; NFy; NMz];

% sum of reaction tangential friction forces from the strap on the cam
fFx=-fric*NFy;
fFy=fric*NFx;
fMz=(fFy*xcom)-(fFx*ycom)+(Rb2*Talpha2*(exp(fric*real(theta2))-1));
Wforcetan=[fFx; fFy; fMz];

Wforce=Wforceattach+Wforcecorner+Wforcenormal+Wforcetan;

% % moment to resist the strap bent around the cams
% Momenttheta1=-(E*W*(t^3))/(12*Rp1);
% Momenttheta2=(E*W*(t^3))/(12*Rp2);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-alpha1/2);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentalpha2=LivingHingeK*(alpha2/2);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta1+Momenttheta2+Momentalpha1+Momentbeta1+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload]=tensionregime8(Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, vec1, vec2, Lsvec)
% This function calculates the loads on a layer in the tension scenario
% within regime 8

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp1*sin(alpha1/2))-(Rp1*theta1);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls+(2*Rp2*sin(alpha2/2))-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch)*(exp(fric*theta1));
    denominat=(2*Rp1*sin(alpha1/2))+(Rp1*((exp(fric*theta1)-1)/fric))+(Li*exp(fric*theta1));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% force to resist tension in the strap where it attaches to the cam
Locforce1=vec1;
dirforce=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
forceattach=Talpha2*dirforce;
Wforceattach=[forceattach(1,1); forceattach(1,2); ((Locforce1(1,1)*forceattach(1,2))-(Locforce1(1,2)*forceattach(1,1)))];

% reaction force at corner from the strap
Locforce2=vec2;
dirforce=(vec2-vec1)/sqrt(dotFast((vec2-vec1),(vec2-vec1)));
forcecorner1=Talpha2*dirforce;
if(Ls==0)
    perpLoc=[Locforce2(1,2), -Locforce2(1,1)];
    forcecorner2=Talpha2*(perpLoc/sqrt(dotFast(perpLoc,perpLoc)));
else
    forcecorner2=Talpha2*(Lsvec/Ls);
end
forcecorner=forcecorner1+forcecorner2;
Wforcecorner=[forcecorner(1,1); forcecorner(1,2); ((Locforce2(1,1)*forcecorner(1,2))-(Locforce2(1,2)*forcecorner(1,1)))];

Wforce=Wforceattach+Wforcecorner;

% % moment to resist the strap bent around the cams
% Momenttheta1=-(E*W*(t^3))/(12*Rp1);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% slopevec1=(vec1-vec2)/sqrt(dotFast((vec1-vec2),(vec1-vec2)));
% anbetween=acos(dotFast(slopevec1,Lsvec/Ls));
%
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-alpha1/2);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentalpha2=LivingHingeK*(anbetween);
%     Momentbeta2=LivingHingeK*(alpha2/2);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentalpha2=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta1+Momentalpha1+Momentbeta1+Momentalpha2+Momentbeta2)];

Wload=Wforce; %+WMoment;
function [Wload]=tensionregime9(xcom, ycom, Rp1, Rp2, t, W, del, dcut1, dcut2, beta1, beta2, Ec, Et, v, fric, theta1, Loc, Lsvec)
% This function calculates the loads on a layer in the tension scenario
% within regime 9

E=(Ec+Et)/2;
Rb1=Rp1-(t/2);
Rb2=Rp2-(t/2);
alpha1=acos((Rb1-dcut1)/Rb1);
alpha2=acos((Rb2-dcut2)/Rb2);
L=(Rp1*((pi/2)-alpha1-beta1))+(2*Rp1*sin(alpha1/2))+(Rp2*((pi/2)-alpha2-beta2))+(2*Rp2*sin(alpha2/2));
Lo=L-del;
Li=Lo-(2*Rp1*sin(alpha1/2))-(Rp1*theta1);
Ls=sqrt(dotFast(Lsvec,Lsvec));
stretch=cutdec(Ls-Li,10);
if (stretch > 0)
    numerat=(Et*t*W)*(stretch)*(exp(fric*theta1));
    denominat=(2*Rp1*sin(alpha1/2))+(Rp1*((exp(fric*theta1)-1)/fric))+(Li*exp(fric*theta1));
    Talpha2=(numerat/denominat);
else
    Talpha2=0;
end

% force to resist tension in the strap
Locforce=Loc+Lsvec;
if(Ls==0)
    perpLoc=[Locforce(1,2), -Locforce(1,1)];
    force=Talpha2*(perpLoc/sqrt(dotFast(perpLoc,perpLoc)));
else
    force=Talpha2*(Lsvec/Ls);
end
Wforce=[force(1,1); force(1,2); ((Locforce(1,1)*force(1,2))-(Locforce(1,2)*force(1,1)))]; % tension force from strap

% % moment to resist the strap bent around the cams
% Momenttheta1=-(E*W*(t^3))/(12*Rp1);
%
% % living hinge moment to resist the strap bent around sharp corners
% Con=1;
% LivingHingeK=Con*((E*W*(t^2))/(2*(1+v)))*((1/3)-((0.21*t/W)*(1-((t^4)/(12*(W^4))))));
% vec1=Loc+Lsvec;
% R2vec=vec1-[xcom, ycom];
% perpvec=[R2vec(1,2), -R2vec(1,1)];
% unitperp=perpvec/sqrt(dotFast(perpvec,perpvec));
% if(unitperp(1,2)>0)
%     anperp=(2*pi)-acos(dotFast([1,0],unitperp));
% else
%     anperp=acos(dotFast([1,0],unitperp));
% end
% unitLs=-Lsvec/Ls;
% if(unitLs(1,2)>0)
%     anLs=(2*pi)-acos(dotFast([1,0],unitLs));
% else
%     anLs=acos(dotFast([1,0],unitLs));
% end
% anbetween=anLs-anperp;
%
% if(Talpha2>0)
%     Momentalpha1=LivingHingeK*(-alpha1/2);
%     Momentbeta1=LivingHingeK*(-alpha1/2);
%     Momentbeta2=LivingHingeK*(anbetween);
% else
%     Momentalpha1=0;
%     Momentbeta1=0;
%     Momentbeta2=0;
% end
%
% WMoment=[0; 0; (Momenttheta1+Momentalpha1+Momentbeta1+Momentbeta2)];

Wload=Wforce; %+WMoment;
function wrapstrapconnectlattice(x,y,r,t,Z,beta,C,phi,rotangle,veccam1)
% draws front layer strap sector where x and y are center coordinates and r is
% the base circle radius and t is the strap thickness. Z is the angle
% from the horizontal to the perpendicular line of the joint. beta is defined
% in the figure, C is how many strap thickness lengths the strap attaches to
% the cams, and phi is how much the cam rotates

inc=100;
connectprt=(C*t/r);

% connector sector points
an1=Z-beta+connectprt+phi;
an2=Z-beta+phi;
thc = an1:-abs(an1-an2)/inc:an2;
xunit1c = (r+t)*cos(thc);
xunitc = [xunit1c, 0]+x;
yunit1c = (r+t)*sin(thc);
yunitc = [yunit1c, 0]+y;
zunitc=zeros(size(yunitc));

%rotate them back to the global coordinate system
xunitg=(xunitc*cos(rotangle))-(yunitc*sin(rotangle))+veccam1(1,1);
yunitg=(xunitc*sin(rotangle))+(yunitc*cos(rotangle))+veccam1(2,1);
xunitc=xunitg;
yunitc=yunitg;

patch(xunitc, yunitc, zunitc, [0.6 0.6 0.6], 'LineStyle', 'none');
function wrapstrapthetalattice(x,y,r,t,Z,dcut,beta,theta,phi,rotangle,veccam1)
% draws back layer strap sector where x and y are center coordinates and r is
% the base circle radius and t is the strap thickness. Z is the angle
% from the horizontal to the perpendicular line of the joint. dcut is the
% fabricated cut by the strap that defines alpha. beta and theta are defined
% in the figure, and phi is how much the cam rotates

inc=100;
alpha=acos((r-dcut)/r);

% theta sector points
an3=Z-beta-alpha+phi;
an4=Z-beta-alpha-theta+phi;
ththe = an3:-abs(an3-an4)/inc:an4;
xunit1the = (r+t)*cos(ththe);
xunitthe = [xunit1the, 0]+x;
yunit1the = (r+t)*sin(ththe);
yunitthe = [yunit1the, 0]+y;
zunitthe=zeros(size(yunitthe));

%rotate them back to the global coordinate system
xunitg=(xunitthe*cos(rotangle))-(yunitthe*sin(rotangle))+veccam1(1,1);
yunitg=(xunitthe*sin(rotangle))+(yunitthe*cos(rotangle))+veccam1(2,1);
xunitthe=xunitg;
yunitthe=yunitg;

patch(xunitthe, yunitthe, zunitthe, [0.6 0.6 0.6], 'LineStyle', 'none');
function drawlayerundeformed(FullJointMat,LocMat,jointnum, rm, C, boltr, t, del, dcut1, dcut2)

RelativeLoc=FullJointMat(:,1:2);
for i=1:rm
    for j=1:(2*jointnum)
        if(FullJointMat(j,1)==i)
            an1=FullJointMat(j,3)+(pi/2)-FullJointMat(j,4);
            vecattach1=FullJointMat(j,6)*[cos(an1), sin(an1)];
            vecstrap=(FullJointMat(j,8)-del)*[sin(an1), -cos(an1)];
            vecattach2=FullJointMat(j,7)*[cos(an1), sin(an1)];
            RelativeLoc(j,3:4)=vecattach1+vecstrap+vecattach2;
        end
    end
end

UpdatedLocMat=zeros(rm,4);
UpdatedLocMat(:,4)=LocMat(:,4);
UpdatedLocMat(1,:)=LocMat(1,:);

kk=1;
history=zeros(1,rm);
history(1,1)=kk;
subholder=0;
ccc=0;
while(ismember(0,history)==1)
    cc=0;
    subhistory=0;
    for j=1:(2*jointnum)
        if(FullJointMat(j,1)==kk)
            cc=cc+1;
            UpdatedLocMat(FullJointMat(j,2),1)=UpdatedLocMat(kk,1)+FullJointMat(j,5)-FullJointMat(j,4);
            Relloc=[cos(UpdatedLocMat(kk,1)), -sin(UpdatedLocMat(kk,1)); sin(UpdatedLocMat(kk,1)), cos(UpdatedLocMat(kk,1))]*[RelativeLoc(j,3); RelativeLoc(j,4)];
            UpdatedLocMat(FullJointMat(j,2),2)=UpdatedLocMat(kk,2)+Relloc(1,1);
            UpdatedLocMat(FullJointMat(j,2),3)=UpdatedLocMat(kk,3)+Relloc(2,1);
            subhistory(1,cc)=FullJointMat(j,2);
        end
    end
    ccc=ccc+1;
    subholder(ccc,1:cc)=subhistory;
    for ee=1:ccc
        for hh=1:cc
            if(ismember(subholder(ee,hh),history)==0)
                kk=subholder(ee,hh);
                history(1,kk)=kk;
            end
        end
    end
end

LocMat=UpdatedLocMat;

for i=1:jointnum
    phi1=LocMat(FullJointMat(i,1),1);
    phi2=LocMat(FullJointMat(i,2),1);
    veccam2=[LocMat(FullJointMat(i,2),2); LocMat(FullJointMat(i,2),3)];
    veccam1=[LocMat(FullJointMat(i,1),2); LocMat(FullJointMat(i,1),3)];
    difvec=veccam2-veccam1;
    rotangle=FullJointMat(i,3)-(pi/2)+phi1;
    rotvec=[cos(rotangle), -sin(rotangle); sin(rotangle), cos(rotangle)]\difvec;
    Rp1=FullJointMat(i,6);
    Rp2=FullJointMat(i,7);
    Rb1=Rp1-(t/2);
    Rb2=Rp2-(t/2);
    beta1=FullJointMat(i,4);
    beta2=FullJointMat(i,5);
    phi=phi2-phi1;
    xcom=rotvec(1,1);
    ycom=rotvec(2,1);
    
    veclocstr1=Rp1*[cos(pi-beta1), sin(pi-beta1)];
    vecper=[veclocstr1(1,2), -veclocstr1(1,1)];
    unitstrapvec=(vecper)/sqrt(dot(vecper,vecper));
    Lsvec=(FullJointMat(i,8)-del)*unitstrapvec;
    
    wrapstrapconnectlattice(0,0,Rb1,t,pi,beta1,C,0,rotangle,veccam1)
    wrapstrapconnectlattice(xcom,ycom,Rb2,t,2*pi,beta2,C,phi,rotangle,veccam1)
    straplattice(veclocstr1,Lsvec,t,0,rotangle,veccam1)
end
for i=1:jointnum
    phi1=LocMat(FullJointMat(i,1),1);
    phi2=LocMat(FullJointMat(i,2),1);
    veccam2=[LocMat(FullJointMat(i,2),2); LocMat(FullJointMat(i,2),3)];
    veccam1=[LocMat(FullJointMat(i,1),2); LocMat(FullJointMat(i,1),3)];
    difvec=veccam2-veccam1;
    rotangle=FullJointMat(i,3)-(pi/2)+phi1;
    rotvec=[cos(rotangle), -sin(rotangle); sin(rotangle), cos(rotangle)]\difvec;
    Rp1=FullJointMat(i,6);
    Rp2=FullJointMat(i,7);
    Rb1=Rp1-(t/2);
    Rb2=Rp2-(t/2);
    beta1=FullJointMat(i,4);
    beta2=FullJointMat(i,5);
    phi=phi2-phi1;
    xcom=rotvec(1,1);
    ycom=rotvec(2,1);
    
    circlelattice2(0,0,Rb1,boltr,0,rotangle,veccam1)
    circlelattice2(xcom,ycom,Rb2,boltr,phi,rotangle,veccam1)
end
for i=1:jointnum
    phi1=LocMat(FullJointMat(i,1),1);
    phi2=LocMat(FullJointMat(i,2),1);
    veccam2=[LocMat(FullJointMat(i,2),2); LocMat(FullJointMat(i,2),3)];
    veccam1=[LocMat(FullJointMat(i,1),2); LocMat(FullJointMat(i,1),3)];
    difvec=veccam2-veccam1;
    rotangle=FullJointMat(i,3)-(pi/2)+phi1;
    rotvec=[cos(rotangle), -sin(rotangle); sin(rotangle), cos(rotangle)]\difvec;
    Rp1=FullJointMat(i,6);
    Rp2=FullJointMat(i,7);
    Rb1=Rp1-(t/2);
    Rb2=Rp2-(t/2);
    beta1=FullJointMat(i,4);
    beta2=FullJointMat(i,5);
    phi=phi2-phi1;
    xcom=rotvec(1,1);
    ycom=rotvec(2,1);
    
    fabcutlattice2(0,0,Rb1,pi,dcut1,beta1,0,rotangle,veccam1)
    fabcutlattice2(xcom,ycom,Rb2,2*pi,dcut2,beta2,phi,rotangle,veccam1)
end
function [PitchLoc]=PitchLocation(phi, xcom, ycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v)
%This function calculates the load that needs to be imparted on cam 2 to
%move it to the positions xcom and ycom and the orientation phi.

bphi=-phi;
bxcom=-xcom;
bycom=ycom;

% front layer
[tension, regime, attach1, attalpha1, vec1, vec2, theta1, theta2, Loc, Lsvec, b, Fc]=kinematics(phi, xcom, ycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v);

% back layer
[btension, bregime, battach1, battalpha1, bvec1, bvec2, btheta1, btheta2, bLoc, bLsvec, bb, bFc]=kinematics(bphi, bxcom, bycom, Rp1, Rp2, t, W, dcut1, dcut2, beta1, beta2, Ec, v);

Ls=sqrt(dotFast(Lsvec,Lsvec));
if(cutdec(Ls,10)==0)
    perp=[Loc(1,2), -Loc(1,1)];
    unitLsvec=perp/sqrt(dotFast(perp,perp));
else
    unitLsvec=Lsvec/sqrt(dotFast(Lsvec,Lsvec));
end

bLs=sqrt(dotFast(bLsvec,bLsvec));
if(cutdec(bLs,10)==0)
    bperp=[bLoc(1,2), -bLoc(1,1)];
    unitbLsvec=bperp/sqrt(dotFast(bperp,bperp));
else
    unitbLsvec=bLsvec/sqrt(dotFast(bLsvec,bLsvec));
end
bbLoc=[-bLoc(1,1), bLoc(1,2)];
unitbbLsvec=[-unitbLsvec(1,1), unitbLsvec(1,2)];

if(cutdec(acos(dotFast(unitLsvec,unitbbLsvec)),10)==0 || cutdec(acos(dotFast(unitLsvec,unitbbLsvec)),10)==pi)
    if(tension==1)
        PitchLoc=Loc;
    elseif(regime==4 && Ls<(b/2))
        PitchLoc=Loc;
    elseif(bregime==4 && bLs<(bb/2))
        PitchLoc=bLoc;
    elseif(regime==8 && Ls<(b/2))
        PitchLoc=Loc+Lsvec;
    elseif(bregime==8 && bLs<(bb/2))
        PitchLoc=bLoc+bLsvec;
    else
        PitchLoc=bbLoc+((bb/2)*unitbbLsvec);
    end
else
    tt=(bbLoc(1,1)-Loc(1,1))/(unitLsvec(1,1)-unitbbLsvec(1,1));
    PitchLoc=Loc+(unitLsvec*tt);
end

%% Misc. custom functions
function grayImage = grayscale(rgbImage)
redChannel = rgbImage(:, :, 1);
greenChannel = rgbImage(:, :, 2);
blueChannel = rgbImage(:, :, 3);
% Do the weighted sum.
grayImage = .299*double(redChannel) + ...
    .587*double(greenChannel) + ...
    .114*double(blueChannel);
grayImage = uint8(grayImage);
function varargout = cinput(arg1)
%GINPUT Graphical input from mouse.
%   [X,Y] = CINPUT(N) gets N points from the current axes and returns
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse.  Data points are entered by pressing
%   a mouse button or any key on the keyboard except carriage return,
%   which terminates the input before N points are entered.
%
%   [X,Y] = CINPUT gathers an unlimited number of points until the
%   return key is pressed.
%
%   [X,Y,BUTTON] = CINPUT(N) returns a third result, BUTTON, that
%   contains a vector of integers specifying which mouse button was
%   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
%   was used.
%
%   Examples:
%       [x,y] = cinput;
%
%       [x,y] = cinput(5);
%
%       [x, y, button] = cinput(1);
%
%   See also GTEXT, WAITFORBUTTONPRESS.

%   Copyright 1984-2013 The MathWorks, Inc.

x = []; y = []; userInputs = [];

if ~matlab.ui.internal.isFigureShowEnabled
    error(message('MATLAB:hg:NoDisplayNoFigureSupport', 'ginput'))
end

% Check inputs
if nargin == 0
    how_many = -1;
else
    how_many = arg1;
    
    if ~isPositiveScalarIntegerNumber(how_many)
        error(message('MATLAB:ginput:NeedPositiveInt'))
    end
    if how_many == 0
        % If input argument is equal to zero points,
        % give a warning and return empty for the outputs.
        warning (message('MATLAB:ginput:InputArgumentZero'));
    end
end

% Get figure
fig = gcf;
figure(fig);

% Make sure the figure has an axes
gca(fig);

% Setup the figure to disable interactive modes and activate pointers.
initialState = setupFcn(fig);

% onCleanup object to restore everything to original state in event of
% completion, closing of figure errors or ctrl+c.
c = onCleanup(@() restoreFcn(initialState));

drawnow
char = 0;

while how_many ~= 0
    try
        mode = waitForUserInput(fig);
    catch %#ok<CTCH>
        cleanup(c);
        if(ishghandle(fig))
            error(message('MATLAB:ginput:Interrupted'));
        else
            error(message('MATLAB:ginput:FigureDeletionPause'));
        end
    end
    
    
    % Make sure figure has not been closed
    checkFigureAvailable();
    
    if (isCorrectFigure(fig))
        switch mode
            case 'key'
                char = get(fig, 'CurrentCharacter');
                curUserInput = abs(get(fig, 'CurrentCharacter'));
            case 'mouse'
                curUserInput = get(fig, 'SelectionType');
                if strcmp(curUserInput,'open')
                    curUserInput = 1;
                elseif strcmp(curUserInput,'normal')
                    curUserInput = 1;
                elseif strcmp(curUserInput,'extend')
                    curUserInput = 2;
                elseif strcmp(curUserInput,'alt')
                    curUserInput = 3;
                else
                    error(message('MATLAB:ginput:InvalidSelection'))
                end
        end
        axes_handle = get(fig,'CurrentAxes');
        pt = get(axes_handle, 'CurrentPoint');
        
        how_many = how_many - 1;
        
        if(char == 13)
            % if the return key was pressed, char will == 13,
            % and that's our signal to break out of here whether
            % or not we have collected all the requested data
            % points.
            % If this was an early breakout, don't include
            % the <Return> key info in the return arrays.
            % We will no longer count it if it's the last input.
            break;
        end
        
        x = [x;pt(1,1)]; %#ok<AGROW>
        y = [y;pt(1,2)]; %#ok<AGROW>
        userInputs = [userInputs;curUserInput]; %#ok<AGROW>
    end
end

% Cleanup and Restore
cleanup(c);

if nargout == 1
    varargout{1} = [x y];
else
    varargout{1} = x;
end
if nargout > 1
    varargout{2} = y;
end
if nargout > 2
    varargout{3} = userInputs;
end
function valid = isPositiveScalarIntegerNumber(how_many)

valid = ~ischar(how_many) && ...            % is numeric
    isscalar(how_many) && ...           % is scalar
    (fix(how_many) == how_many) && ...  % is integer in value
    how_many >= 0;                      % is positive
function mode = waitForUserInput(fig)
waitfor(fig,'UserData')
% Extract mode to determine if key or mouse was used
mode = get(fig,'UserData');
if ischar(mode)
    ud = strsplit(mode, '_');
    mode = ud{1};
end% Reset user data to prepare for next trigger
set(fig,'UserData',[])
function initialState = setupFcn(fig)

% Store Figure Handle.
initialState.figureHandle = fig;

% Suspend figure functions
initialState.uisuspendState = uisuspend(fig);

% Disabling ^C for edit menu so the only ^C is for interrupting the function
initialState.AcceleratorMenu = findall(fig,'Type','uimenu','Accelerator','C');
set(initialState.AcceleratorMenu,'Accelerator','');

% Extract user data
initialState.PreviousUserData = get(fig,'UserData');

% Set callbacks to distinguish from key and mouse triggers
% Using random numbers to distinguish things like double clicks.
set(fig, 'WindowButtondownFcn', @(~,~)set(fig, 'UserData', ['mouse_' num2str(rand)]));
set(fig, 'WindowKeyPressFcn', @(~,~) set(fig, 'UserData', ['key_' num2str(rand)]));

% Disable Plottools Buttons
initialState.toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn')];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','off');
end

%Setup empty pointer
P = ones(16)+1;
P(1,:) = 1; P(16,:) = 1;
P(:,1) = 1; P(:,16) = 1;
P(1:4,8:9) = 1; P(13:16,8:9) = 1;
P(8:9,1:4) = 1; P(8:9,13:16) = 1;
P(5:12,5:12) = NaN; % Create a transparent region in the center
cdata = P;
hotspot = [9, 9];
set(gcf,'Pointer','custom','PointerShapeCData',cdata,'PointerShapeHotSpot',hotspot)

% Get the initial Figure Units
initialState.fig_units = get(fig,'Units');
function restoreFcn(initialState)
if ishghandle(initialState.figureHandle)
    
    % Figure Units
    set(initialState.figureHandle,'Units',initialState.fig_units);
    
    % Reset user data
    set(initialState.figureHandle,'UserData',initialState.PreviousUserData)
    
    % Enable Ctrl+c
    set(initialState.AcceleratorMenu,'Accelerator','C');
    
    % Plottools Icons
    if ~isempty(initialState.toolbar) && ~isempty(initialState.ptButtons)
        set (initialState.ptButtons(1),'Enable',initialState.ptState{1});
        set (initialState.ptButtons(2),'Enable',initialState.ptState{2});
    end
    
    % UISUSPEND
    uirestore(initialState.uisuspendState);
end
function updateCrossHair(fig, crossHair)
% update cross hair for figure.
gap = 3; % 3 pixel view port between the crosshairs
cp = hgconvertunits(fig, [fig.CurrentPoint 0 0], fig.Units, 'pixels', fig);
cp = cp(1:2);
figPos = hgconvertunits(fig, fig.Position, fig.Units, 'pixels', fig.Parent);
figWidth = figPos(3);
figHeight = figPos(4);

% Early return if point is outside the figure
if cp(1) < gap || cp(2) < gap || cp(1)>figWidth-gap || cp(2)>figHeight-gap
    return
end

set(crossHair, 'Visible', 'on');
thickness = 1; % 1 Pixel thin lines.
set(crossHair(1), 'Position', [0 cp(2) cp(1)-gap thickness]);
set(crossHair(2), 'Position', [cp(1)+gap cp(2) figWidth-cp(1)-gap thickness]);
set(crossHair(3), 'Position', [cp(1) 0 thickness cp(2)-gap]);
set(crossHair(4), 'Position', [cp(1) cp(2)+gap thickness figHeight-cp(2)-gap]);
function checkFigureAvailable()
% See if root has children
figchildren = allchild(0);
if isempty(figchildren)
    error(message('MATLAB:ginput:FigureUnavailable'));
end
function valid = isCorrectFigure(fig)
% g467403 - ginput failed to discern clicks/keypressed on the figure it was
% registered to if the figure's handleVisibility was set to 'callback'
figchildren = allchild(0);
% Select figure at top
ptr_fig = figchildren(1);
% Check if they match
valid = isequal(ptr_fig,fig);
function cleanup(c)
if isvalid(c)
    delete(c);
end
