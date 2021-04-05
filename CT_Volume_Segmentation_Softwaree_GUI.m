function varargout = CT_Volume_Segmentation_Softwaree_GUI(varargin)
% This file is developed by HKUST ECE FYP 2020-2021 YW03a-20 - CT Image Segmentation
% Student Lee Man View Raymond (mvrlee@connect.ust.hk)
% Student Au Pui Sze (psauab@connect.ust.hk)
% Supervised by Professor Yu Weichuan and his assistant Zhang Xuechen, 

% CT_VOLUME_SEGMENTATION_SOFTWAREE_GUI MATLAB code for CT_Volume_Segmentation_Softwaree_GUI.fig
%      CT_VOLUME_SEGMENTATION_SOFTWAREE_GUI, by itself, creates a new CT_VOLUME_SEGMENTATION_SOFTWAREE_GUI or raises the existing
%      singleton*.
%
%      H = CT_VOLUME_SEGMENTATION_SOFTWAREE_GUI returns the handle to a new CT_VOLUME_SEGMENTATION_SOFTWAREE_GUI or the handle to
%      the existing singleton*.
%
%      CT_VOLUME_SEGMENTATION_SOFTWAREE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CT_VOLUME_SEGMENTATION_SOFTWAREE_GUI.M with the given input arguments.
%
%      CT_VOLUME_SEGMENTATION_SOFTWAREE_GUI('Property','Value',...) creates a new CT_VOLUME_SEGMENTATION_SOFTWAREE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CT_Volume_Segmentation_Softwaree_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CT_Volume_Segmentation_Softwaree_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CT_Volume_Segmentation_Softwaree_GUI

% Last Modified by GUIDE v2.5 17-Mar-2021 14:10:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CT_Volume_Segmentation_Softwaree_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CT_Volume_Segmentation_Softwaree_GUI_OutputFcn, ...
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

% --- Executes just before CT_Volume_Segmentation_Softwaree_GUI is made visible.
function CT_Volume_Segmentation_Softwaree_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CT_Volume_Segmentation_Softwaree_GUI (see VARARGIN)

% Create variables under the handles structure to save parameters for
% sharing across different callback functions
handles.selected_input_path = '';
handles.selected_input = [];
handles.segmented_input = [];
handles.previous_slide_no = 0;
handles.previous_iter_num = 40;
handles.total_no_of_slides = 0;
handles.time_interval = 100;
handles.phi = [];
handles.segmented_flag = false;
set(handles.PlayToggleBtn, 'String', 'Play');
% Choose default command line output for CT_Volume_Segmentation_Softwaree_GUI
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CT_Volume_Segmentation_Softwaree_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Default function generated function by GUIDE
% --- Outputs from this function are returned to the command line.
function varargout = CT_Volume_Segmentation_Softwaree_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Default function generated function by GUIDE
% --- Executes during object creation, after setting all properties.
function InputEditTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputEditTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in BrowseBtn.
function BrowseBtn_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% allow user to select *.nii files as the segmentation input
[file,path] = uigetfile('*.nii','Select input for CT Volume Segmentation');
% Do nothing if the browsing action is cancelled
if isequal(file,0)
else
    % concatenate the file path and file as one string and display it in the
    % input edit text box
    selected_input_path = strcat(path,file);
    set(handles.InputEditTextBox, 'String', selected_input_path);
    % read the input as a 3d array
    selected_input = niftiread(selected_input_path);
    total_no_of_slides = size(selected_input,3);
    % display the total number of slides in the total slide number static
    % text box
    total_no_of_slides_str = sprintf('/%d', total_no_of_slides);
    set(handles.TotalSlideNumStaticTextBox, 'String', total_no_of_slides_str);
    % use the middle slide number as the default displaying slide
    middle_slide_no = round(total_no_of_slides / 2);
    set(handles.CurrentSlideNumberEditTextBox, 'String', middle_slide_no);
    % set the max property of the slider bar so that it wont display slides
    % out of range of the input
    set(handles.SlideNumberSlider, 'Max', total_no_of_slides);
    set(handles.SlideNumberSlider, 'Value', middle_slide_no);
    % update the parameters in the handles structure
    handles.previous_slide_no = middle_slide_no;
    handles.selected_input_path = selected_input_path;
    handles.selected_input = selected_input;
    handles.total_no_of_slides = total_no_of_slides;
    handles.segmented_flag = false;
    handles.phi = [];
    % display the middle slide 
    imshow(squeeze(selected_input(:,:,middle_slide_no,:)), [double(min(handles.selected_input(:))), double(max(handles.selected_input(:)))], 'Parent', handles.DisplayPanel);
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SegmentBtn.
function SegmentBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% do nothing if the input is empty
if isequal(handles.selected_input_path, '')
else
    totaltime = 0;
    t0 = cputime;

    % reading image
    IMG = handles.selected_input;
    IMG = double(IMG);

    % pre-processing
    % IMG = 255*mat2gray(IMG);  % Image normalization
    IMG = medfilt3(IMG);  % ok
    % IMG = histeq(IMG);    % histogram equalization, blank, not work

    % delete black slices
    [x,y,z] = size(IMG);
    boolX = false(1,x);
    boolY = false(1,y);
    boolZ = false(1,z);
    boolinX = false;
    boolinY = false;
    boolinZ = false;
    EDGEWIDTH = 5;
    for i = 1:x
        if boolinX == false
            if all(all(IMG(i+EDGEWIDTH,:,:) == 0))
                boolX(i) = true;
            elseif (i < EDGEWIDTH)
                boolinX = true;
            elseif ~all(all(IMG(i-EDGEWIDTH,:,:) == 0))
                boolinX = true;
            end
        elseif (i > EDGEWIDTH)
            if all(all(IMG(i-EDGEWIDTH,:,:) == 0))
                boolX(i) = true;
            end
        end
    end
    for i = 1:y
        if boolinY == false
            if all(all(IMG(:,i+EDGEWIDTH,:) == 0))
                boolY(i) = true;
            elseif (i < EDGEWIDTH)
                boolinY = true;
            elseif ~all(all(IMG(:,i-EDGEWIDTH,:) == 0))
                boolinY = true;
            end
        elseif (i > EDGEWIDTH)
            if all(all(IMG(:,i-EDGEWIDTH,:) == 0))
                boolY(i) = true;
            end
        end
    end
    for i = 1:z
        if boolinZ == false
            if all(all(IMG(:,:,i+EDGEWIDTH) == 0))
                boolZ(i) = true;
            elseif (i < EDGEWIDTH)
                boolinZ = true;
            elseif ~(all(all(IMG(:,:,i-EDGEWIDTH) == 0)))
                boolinZ = true;
            end
        elseif (i > EDGEWIDTH)
            if all(all(IMG(:,:,i-EDGEWIDTH) == 0))
                boolZ(i) = true;
            end
        end
    end
    IMG(boolX,:,:) = [];
    IMG(:,boolY,:) = [];
    IMG(:,:,boolZ) = [];
        
    handles.segmented_input = IMG;
    
    % constant initialize
    ITR_INNER = 10;
    ITR_OUTER = str2double(get(handles.IterNumEditTextBox, 'String'));
    EPSILON = 1;
    SIGMA = 4;
    ONE_K = imgaussfilt3(IMG,SIGMA);
    K = fspecial('gaussian',2*ceil(2*SIGMA)+1,SIGMA);
    LAMBDA = 1;
    COEF_AL = 0.001*255^2;
    COEF_DRLSE = 1;
    TIMESTEP = 0.1;
    
    % initialize level set function phi and bias field b
    [x,y,z] = size(IMG);
    phi = ones(size(IMG));
    b = phi;
    phi(ceil(x/2-25):ceil(x/2+25),ceil(y/2-25):ceil(y/2+25),ceil(z/2-25):ceil(z/2+25)) = -1;
    m = UpdateM(phi,EPSILON);    

    terminal_output = get(handles.TerminalEditTextBox, 'String');
    if isempty(terminal_output)
        terminal_output = {};
    end
    str = 'Starting the segmentation process...';
    terminal_output = [str; terminal_output];
    set(handles.TerminalEditTextBox, 'String', terminal_output);
    drawnow
    
    for i = 1:ITR_OUTER
        terminal_output = get(handles.TerminalEditTextBox, 'String');
        if isempty(terminal_output)
            terminal_output = {};
        end
        % reflect segmentation progress in the terminal
        str = sprintf('Current iteration number = %d/%d', i, ITR_OUTER);
        terminal_output = [str; terminal_output];
        set(handles.TerminalEditTextBox, 'String', terminal_output);
        % allow the user interface to update 
        drawnow
        % update (b*K) and (b^2*K)
        bk = convn(b,K,'same');
        b2k = convn(b.^2,K,'same');
        % energy minimization with respect to c
        c = UpdateC(bk,b2k,IMG,m);
        % energy minimization with respect to phi
        phi = UpdatePHI(bk,b2k,IMG,m,c,phi,ONE_K,ITR_INNER,EPSILON,LAMBDA,COEF_AL,COEF_DRLSE,TIMESTEP);
        % update m
        m = UpdateM(phi,EPSILON);
        % energy minimization with respect to b
        b = UpdateB(IMG,m,c,K);
    end

    % visualize the result in UI
    % figure;
    % imshow3D(IMG,phi); 
    % update the parameters in the handles structure
    handles.segmented_flag = true;
    handles.phi = phi;
    total_no_of_slides = size(IMG,3);
    handles.total_no_of_slides = total_no_of_slides;
    total_no_of_slides_str = sprintf('/%d', total_no_of_slides);
    set(handles.TotalSlideNumStaticTextBox, 'String', total_no_of_slides_str);
    middle_slide_no = round(total_no_of_slides / 2);
    set(handles.CurrentSlideNumberEditTextBox, 'String', middle_slide_no);
    set(handles.SlideNumberSlider, 'Max', total_no_of_slides);
    set(handles.SlideNumberSlider, 'Value', middle_slide_no);
    guidata(hObject,handles);
    
    val = round(str2double(get(handles.CurrentSlideNumberEditTextBox, 'String')));
    imshow(squeeze(handles.segmented_input(:,:,val,:)), [double(min(handles.segmented_input(:))), double(max(handles.segmented_input(:)))], 'Parent', handles.DisplayPanel);
    % display the level set function on top of the input
    hold on;
    contour(squeeze(phi(:,:,val)),[0 0],'r')
    hold off;
    
    % record the time used for the segmentation and show it in the terminal
    t1 = cputime;
    totaltime = totaltime + t1 - t0;
    terminal_output = get(handles.TerminalEditTextBox, 'String');
    if isempty(terminal_output)
        terminal_output = {};
    end
    str = sprintf('The segmentation process is finished, total time used (s) = %d', totaltime);
    terminal_output = [str; terminal_output];
    set(handles.TerminalEditTextBox, 'String', terminal_output);
    % totaltime
end

% energy minimization with respect to c
function c = UpdateC(bk,b2k,IMG,m)
for i = 1:2
    term1 = bk.*IMG.*m(:,:,:,i);
    term2 = b2k.*m(:,:,:,i);
    c(i) = sum(term1,'all')/sum(term2,'all');
end

% energy minimization with respect to phi
function phi = UpdatePHI(bk,b2k,IMG,m,c,phi,ONE_K,ITR_INNER,EPSILON,LAMBDA,COEF_AL,COEF_DRLSE,TIMESTEP)
e = zeros(size(m));
for i = 1:2
    term1 = IMG.^2.*ONE_K;
    term2 = 2.*c(i).*IMG.*bk;
    term3 = c(i).^2*b2k;
    e(:,:,:,i) = term1-term2+term3;
end
for i = 1:ITR_INNER
    DD = DiracDelta(phi,EPSILON);
    CVT = Curvature(phi);
    phi = NeumannBoundCond(phi);
    globalterm = CVT-LAMBDA.*((c(1)-IMG).^2-(c(2)-IMG).^2);
    term1 = -LAMBDA.*(e(:,:,:,1)-e(:,:,:,2));
    term2 = COEF_AL.*DD.*CVT;
    term3 = COEF_DRLSE.*(4.*del2(phi)-CVT);
    phi = phi+TIMESTEP*globalterm;
    phi = phi+TIMESTEP*(term1+term2+term3);
end

% energy minimization with respect to b
function b = UpdateB(IMG,m,c,K)
J1 = zeros(size(IMG));
J2 = J1;
for i = 1:2
    J1 = J1+c(i).*m(:,:,:,i);
    J2 = J2+c(i).^2.*m(:,:,:,i);
end
b = convn(IMG.*J1,K,'same')./convn(J2,K,'same');

function m = UpdateM(phi,EPSILON)
m(:,:,:,1) = Heaviside(phi,EPSILON);
m(:,:,:,2) = 1-Heaviside(phi,EPSILON);

function phi = NeumannBoundCond(phi)
PADSIZE = 5;
[x,y,z] = size(phi);
phi = padarray(phi,[PADSIZE,PADSIZE,PADSIZE]);
phi([1 x],[1 y],[1 z]) = phi([3 x-2],[3 y-2],[3 z-2]);  
phi([1 x],2:end-1,2:end-1) = phi([3 x-2],2:end-1,2:end-1);          
phi(2:end-1,[1 y],2:end-1) = phi(2:end-1,[3 y-2],2:end-1);
phi(2:end-1,2:end-1,[1 z]) = phi(2:end-1,2:end-1,[3 z-2]);
phi = phi(PADSIZE+1:end-PADSIZE,PADSIZE+1:end-PADSIZE,PADSIZE+1:end-PADSIZE);

function CVT = Curvature(phi)
[phi_x,phi_y,phi_z] = gradient(phi);
length = sqrt(phi_x.^2+phi_y.^2+phi_z.^2+1e-10);
Fx = phi_x./length;
Fy = phi_y./length;
Fz = phi_z./length;
CVT = divergence(Fx,Fy,Fz);

function m = Heaviside(phi,EPSILON)    
m = 0.5*(1+(2/pi)*atan(phi./EPSILON));

function DD = DiracDelta(phi,EPSILON)
DD = (EPSILON/pi)./(EPSILON^2.+phi.^2);

% --- Executes on slider movement.
function SlideNumberSlider_Callback(hObject, eventdata, handles)
% hObject    handle to SlideNumberSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Do nothing if the input is empty
if isequal(handles.selected_input_path, '')
else
    % get the slider value, round off to the nearest integer, and update
    % the current slide number edit text box
    slide_no = round(get(handles.SlideNumberSlider, 'Value'));
    set(handles.CurrentSlideNumberEditTextBox, 'String', num2str(slide_no));
    if ~handles.segmented_flag
        imshow(squeeze(handles.selected_input(:,:,slide_no,:)), [double(min(handles.selected_input(:))), double(max(handles.selected_input(:)))], 'Parent', handles.DisplayPanel);
    else
        imshow(squeeze(handles.segmented_input(:,:,slide_no,:)), [double(min(handles.segmented_input(:))), double(max(handles.segmented_input(:)))], 'Parent', handles.DisplayPanel);
    end
    % display the level set function on top of the input
    if (handles.segmented_flag)
        hold on;
        contour(squeeze(handles.phi(:,:,slide_no)),[0 0],'r')
        hold off;
    end
    % update the parameters in the handles structure
    handles.previous_slide_no = slide_no;
    guidata(hObject,handles);
end

% Default function generated by GUIDE
% --- Executes during object creation, after setting all properties.
function SlideNumberSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlideNumberSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function CurrentSlideNumberEditTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentSlideNumberEditTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrentSlideNumberEditTextBox as text
%        str2double(get(hObject,'String')) returns contents of CurrentSlideNumberEditTextBox as a double

% do nothing if the input is empty
if isequal(handles.selected_input_path, '')
else
    % get the value of the current slide number edit text box
    val = round(str2double(get(handles.CurrentSlideNumberEditTextBox, 'String')));
    set(handles.CurrentSlideNumberEditTextBox, 'String', num2str(val));
    % validate user input, must be a number within the range of the input,
    % update the silder bar
    if ~isnan(val) & val >= get(handles.SlideNumberSlider, 'Min') & val <= get(handles.SlideNumberSlider, 'Max')
        set(handles.SlideNumberSlider, 'Value', val);
        handles.previous_slide_no = val;
        guidata(hObject, handles);
        if ~handles.segmented_flag
            imshow(squeeze(handles.selected_input(:,:,val,:)), [double(min(handles.selected_input(:))), double(max(handles.selected_input(:)))], 'Parent', handles.DisplayPanel);
        else
            imshow(squeeze(handles.segmented_input(:,:,val,:)), [double(min(handles.segmented_input(:))), double(max(handles.segmented_input(:)))], 'Parent', handles.DisplayPanel);
        end
        if (handles.segmented_flag)
            hold on;
            contour(squeeze(handles.phi(:,:,val)),[0 0],'r')
            hold off;
        end
    else
        % unsupported input handling, display the previous input
        set(handles.CurrentSlideNumberEditTextBox, 'String', num2str(handles.previous_slide_no)); 
    end
end

% Default function generated by GUIDE
% --- Executes during object creation, after setting all properties.
function CurrentSlideNumberEditTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentSlideNumberEditTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in PlayToggleBtn.
function PlayToggleBtn_Callback(hObject, eventdata, handles)
% hObject    handle to PlayToggleBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlayToggleBtn

% ensure that the play button does not toggle if the input is empty
if isequal(handles.selected_input_path, '')
    set(hObject, 'Value', 0);
else
    slide_no = str2double(get(handles.CurrentSlideNumberEditTextBox, 'String'));
    guidata(hObject, handles);
    % toggle the button
    if get(hObject, 'Value')
        set(handles.PlayToggleBtn, 'String', 'Pause');
    else
        set(handles.PlayToggleBtn, 'String', 'Play');
    end
    while get(hObject,'Value')
        % limit the range of playing to the range of the input
        if slide_no + 1 <= handles.total_no_of_slides
            slide_no = slide_no + 1;
        else
            set(hObject, 'Value', 0);
            set(handles.PlayToggleBtn, 'String', 'Play');
        end
        % update the slider bar and the current slide number edit text box
        % simultaneously while playing the slides
        set(handles.CurrentSlideNumberEditTextBox, 'String', num2str(slide_no));
        set(handles.SlideNumberSlider, 'Value', slide_no);
        if ~handles.segmented_flag
            imshow(squeeze(handles.selected_input(:,:,slide_no,:)), [double(min(handles.selected_input(:))), double(max(handles.selected_input(:)))], 'Parent', handles.DisplayPanel);
        else
            imshow(squeeze(handles.segmented_input(:,:,slide_no,:)), [double(min(handles.segmented_input(:))), double(max(handles.segmented_input(:)))], 'Parent', handles.DisplayPanel);
        end
        if (handles.segmented_flag)
            hold on;
            contour(squeeze(handles.phi(:,:,slide_no)),[0 0],'r')
            hold off;
        end
        % pause with the user defined time interval to allow observable
        % update of the slides
        pause(handles.time_interval/1000);
    end
end

% Executes when a mouse scroll is detected, displays different slides
% according to the mouse scroll movement
function MouseScroll_Callback(hObject, eventdata, handles)
if isequal(handles.selected_input_path, '')
else
    scroll_count = eventdata.VerticalScrollCount;
    slide_no = str2double(get(handles.CurrentSlideNumberEditTextBox, 'String'));
    slide_no = slide_no - scroll_count;
    % limit the range of the mouse scroll to the range of the input
    if slide_no < 1
        slide_no = 1;
    elseif slide_no > handles.total_no_of_slides
        slide_no = handles.total_no_of_slides;
    end
    % update the parameters in the handles structure
    handles.previous_slide_no = slide_no;
    % update the slider bar, the current slide number edit text box and the
    % display panel simultaneously during mouse scroll
    set(handles.SlideNumberSlider, 'Value', slide_no);
    set(handles.CurrentSlideNumberEditTextBox, 'String', num2str(slide_no));
    if ~handles.segmented_flag
        imshow(squeeze(handles.selected_input(:,:,slide_no,:)), [double(min(handles.selected_input(:))), double(max(handles.selected_input(:)))], 'Parent', handles.DisplayPanel);
    else
        imshow(squeeze(handles.segmented_input(:,:,slide_no,:)), [double(min(handles.segmented_input(:))), double(max(handles.segmented_input(:)))], 'Parent', handles.DisplayPanel);
    end
    if (handles.segmented_flag)
        hold on;
        contour(squeeze(handles.phi(:,:,slide_no)),[0 0],'r')
        hold off;
    end
    guidata(hObject, handles);
end

% Executes when a key press on the keyboard is detected
function KeyPress_Callback(hObject, eventdata, handles)
if isequal(handles.selected_input_path, '')
else
    current_object = gco;
    slide_no = str2double(get(handles.CurrentSlideNumberEditTextBox, 'String'));
    switch eventdata.Key
        case {'uparrow', 'leftarrow'}
            % make sure the user did not select the slider bar as the
            % current object, otherwise the key press callback will
            % contradict with the slider bar behaviour
            if ~strcmp(current_object.Tag, 'SlideNumberSlider')  
                % limit the range to the range of input
                if slide_no - 1 < 1
                    slide_no = 1;
                else
                    slide_no = slide_no - 1;
                end
                % update the slider bar, the current slide number edit
                % text box and the display panel simultaneously during key
                % press
                set(handles.SlideNumberSlider, 'Value', slide_no);
                set(handles.CurrentSlideNumberEditTextBox, 'String', num2str(slide_no));
                if ~handles.segmented_flag
                    imshow(squeeze(handles.selected_input(:,:,slide_no,:)), [double(min(handles.selected_input(:))), double(max(handles.selected_input(:)))], 'Parent', handles.DisplayPanel);
                else
                    imshow(squeeze(handles.segmented_input(:,:,slide_no,:)), [double(min(handles.segmented_input(:))), double(max(handles.segmented_input(:)))], 'Parent', handles.DisplayPanel);
                end
                if (handles.segmented_flag)
                    hold on;
                    contour(squeeze(handles.phi(:,:,slide_no)),[0 0],'r')
                    hold off;
                end
                % update the parameters in the handles structure
                handles.previous_slide_no = slide_no;
                guidata(hObject,handles);
            end
        % similar behaviour with the above case
        case {'downarrow', 'rightarrow'}
            if ~strcmp(current_object.Tag, 'SlideNumberSlider')   
                if slide_no + 1 >= handles.total_no_of_slides
                    slide_no = handles.total_no_of_slides;
                else
                    slide_no = slide_no + 1;
                end
                set(handles.SlideNumberSlider, 'Value', slide_no);
                set(handles.CurrentSlideNumberEditTextBox, 'String', num2str(slide_no));
                if ~handles.segmented_flag
                    imshow(squeeze(handles.selected_input(:,:,slide_no,:)), [double(min(handles.selected_input(:))), double(max(handles.selected_input(:)))], 'Parent', handles.DisplayPanel);
                else
                    imshow(squeeze(handles.segmented_input(:,:,slide_no,:)), [double(min(handles.segmented_input(:))), double(max(handles.segmented_input(:)))], 'Parent', handles.DisplayPanel);
                end
                if (handles.segmented_flag)
                    hold on;
                    contour(squeeze(handles.phi(:,:,slide_no)),[0 0],'r')
                    hold off;
                end
                handles.previous_slide_no = slide_no;
                guidata(hObject,handles);
            end
    end
end

function TimeIntervalEditTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to TimeIntervalEditTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeIntervalEditTextBox as text
%        str2double(get(hObject,'String')) returns contents of TimeIntervalEditTextBox as a double
val = str2double(get(handles.TimeIntervalEditTextBox, 'String'));
% validate user input, must be a number
if ~isnan(val)
    handles.time_interval = val;
    guidata(hObject, handles);
else
    set(handles.TimeIntervalEditTextBox, 'String', handles.time_interval);
end

% Default function generated by GUIDE
% --- Executes during object creation, after setting all properties.
function TimeIntervalEditTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeIntervalEditTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Default function generated by GUIDE
% --- Executes during object creation, after setting all properties.
function TerminalEditTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TerminalEditTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function IterNumEditTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to IterNumEditTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IterNumEditTextBox as text
%        str2double(get(hObject,'String')) returns contents of IterNumEditTextBox as a double
val = round(str2double(get(handles.IterNumEditTextBox, 'String')));
set(handles.IterNumEditTextBox, 'String', num2str(val));
% validate user input, must be a number
if ~isnan(val)
    handles.previous_iter_num = val;
    guidata(hObject, handles);
else
    set(handles.IterNumEditTextBox, 'String', num2str(handles.previous_iter_num)); 
end

% Default function generated by GUIDE
% --- Executes during object creation, after setting all properties.
function IterNumEditTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IterNumEditTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
