function varargout = CellPopSim_GUI(varargin)
% CELLPOPSIM_GUI MATLAB code for CellPopSim_GUI.fig
%      CELLPOPSIM_GUI, by itself, creates a new CELLPOPSIM_GUI or raises the existing
%      singleton*.
%
%      H = CELLPOPSIM_GUI returns the handle to a new CELLPOPSIM_GUI or the handle to
%      the existing singleton*.
%
%      CELLPOPSIM_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLPOPSIM_GUI.M with the given input arguments.
%
%      CELLPOPSIM_GUI('Property','Value',...) creates a new CELLPOPSIM_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellPopSim_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellPopSim_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellPopSim_GUI

% Last Modified by GUIDE v2.5 21-Nov-2017 14:10:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellPopSim_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CellPopSim_GUI_OutputFcn, ...
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


% --- Executes just before CellPopSim_GUI is made visible.
function CellPopSim_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellPopSim_GUI (see VARARGIN)

% Choose default command line output for CellPopSim_GUI
handles.output = hObject;

data_controller = varargin{1};
handles.data_controller = data_controller;

set(handles.model_chooser,'String',data_controller.models);

set(handles.plot_type_chooser,'String',{'N(t)','generation(t)'});

handles.Tmax = 10;
set(handles.Tmax_edit,'String',num2str(handles.Tmax));

handles.dt = data_controller.dt;
set(handles.dt_edit,'String',num2str(handles.dt));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CellPopSim_GUI wait for user response (see UIRESUME)
% uiwait(handles.CellPopSim);


% --- Outputs from this function are returned to the command line.
function varargout = CellPopSim_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in test_button.
function test_button_Callback(hObject, eventdata, handles)
% hObject    handle to test_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;

% % % % % % dc.set_kidney_development;
% % % % % % Nini = 1;
% % % % % % %
% % % % % % Tmax = 18*24;
% % % % % % %
% % % % % % [cell_numbers, gen_numbers] = dc.simulate(Nini,Tmax);
% % % % % % figure;
% % % % % % for k=1:numel(dc.cell_types)
% % % % % %      semilogy(dc.t/24,cell_numbers(k,:),'linewidth',2);
% % % % % %      hold on;
% % % % % % end
% % % % % % plot(dc.t/24,sum(cell_numbers,1),'r:','linewidth',2);
% % % % % % %%%%%%%%%%%%%
% % % % % % dpc = [11.5 12 12.5 13.5 14 14.5 15.5 16.5 17.5];
% % % % % % nvox = [420862 1170317 1833695 4488798 6602763 9235365 17267749 24330325 36193224];
% % % % % % um = 1e-6;
% % % % % % voxvol = 4.17*um^3; % um^3
% % % % % % R = 6*um; % 7 microns radius of cell
% % % % % % v = (4/3*pi*R^3)/voxvol; % volume of cell in voxels
% % % % % % semilogy(dpc,nvox/v,'ko-','markersize',12,'linewidth',2);
% % % % % % hold on;
% % % % % % %%%%%%%%%%%%%%
% % % % % % % tau = log(2)/(19.7/24) % tau for Tc = 20 hours
% % % % % % tau = 1/0.8441;
% % % % % % model = exp((dc.t/24)/tau)*(3/13.5);
% % % % % % semilogy(dc.t/24,model,'b.-');
% % % % % % %%%%%%%%%%%%%
% % % % % % hold off;
% % % % % % grid on;
% % % % % % xlabel('days past inception');
% % % % % % ylabel('cell N(t)');
% % % % % % legend([dc.cell_types 'total' 'Lefevre 2017' '2^{t/T_c} for CM']);
% % % % % % 
% % % % % % figure;
% % % % % % for k=1:numel(dc.cell_types)
% % % % % %      plot(dc.t/24,gen_numbers(k,:),'linewidth',2);
% % % % % %      hold on;
% % % % % % end
% % % % % % hold off;
% % % % % % grid on;

model_index = get(handles.model_chooser,'Value');
strings = get(handles.model_chooser,'String');
model = char(strings(model_index));

Nini = 1;

Tmax = str2double(get(handles.Tmax_edit,'String'));

switch model
    case 'cancer'
        dc.set_cancer; %Tmax = 6.5; % days
    case 'cancer_delayed'
        dc.set_cancer_delayed; %Tmax = 13; % days        
    case 'kidney_development'
        dc.set_kidney_development; %Tmax = 11; % days
    case 'disappearance'
        dc.set_disappearance; %Tmax = 13; % days       
end

dc.simulate(Nini,24*Tmax);

if isempty(dc.cell_numbers)
    % cla(handles.main_plot,'reset');
    return;
end

mode_index = get(handles.plot_type_chooser,'Value');
strings = get(handles.plot_type_chooser,'String');
mode = char(strings(mode_index));
%
visualize(handles,dc.cell_numbers,dc.gen_numbers,mode);

% --------------------------------------------------------------------
function Tmax_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Tmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tmax_edit as text
%        str2double(get(hObject,'String')) returns contents of Tmax_edit as a double
value = fix(str2double(get(hObject,'String')));
if isnumeric(value) && value > 0
    handles.Tmax = value;
    guidata(hObject,handles);
else
    value = handles.Tmax;
    set(hObject,'String',num2str(value));
end
uiresume(handles.CellPopSim);

% --- Executes during object creation, after setting all properties.
function Tmax_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in model_chooser.
function model_chooser_Callback(hObject, eventdata, handles)
% hObject    handle to model_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_chooser contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_chooser


% --- Executes during object creation, after setting all properties.
function model_chooser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cancel_simulation.
function cancel_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
% dc.interrupt_simulations;
notify(dc,'interrupt_simulations');

% --------------------------------------------------------------------
function save_current_plot_as_fig_Callback(hObject, eventdata, handles)
% hObject    handle to save_current_plot_as_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uiputfile('*.fig','Save figure');
if pathname == 0 %if the user pressed cancelled, then we exit this callback 
 return 
end 
saveName = fullfile(pathname,filename);

f_new = figure;
ax_new = copyobj(handles.main_plot,f_new);
set(ax_new,'Position','default')
saveas(f_new,saveName,'fig');
close(f_new);

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function visualize(handles,cell_numbers,gen_numbers,mode)
dc = handles.data_controller;
cla(handles.main_plot,'reset');
%
switch mode
    case 'N(t)'
        for k=1:numel(dc.cell_types)
             semilogy(handles.main_plot,dc.t/24,cell_numbers(k,:),'linewidth',2);
             hold on;
        end
        plot(handles.main_plot,dc.t/24,sum(cell_numbers,1),'r:','linewidth',2);
        hold off;                
        legend(handles.main_plot,[dc.cell_types 'total'],'fontsize',14);
        grid(handles.main_plot,'on');
        xlabel(handles.main_plot,'days past inception','fontsize',14);
        ylabel(handles.main_plot,'cell N(t)','fontsize',14);
    case 'generation(t)'
        for k=1:numel(dc.cell_types)            
             gens = gen_numbers(k,:);
             range = gens>0;
             plot(handles.main_plot,dc.t(range)/24,gens(range),'linewidth',2);                          
             hold on;
        end
        hold off;
        legend(handles.main_plot,dc.cell_types,'fontsize',14);
        grid(handles.main_plot,'on');
        xlabel(handles.main_plot,'days past inception','fontsize',14);
        ylabel(handles.main_plot,'cell generation(t)','fontsize',14);
end
% --------------------------------------------------------------------


% --- Executes on selection change in plot_type_chooser.
function plot_type_chooser_Callback(hObject, eventdata, handles)
% hObject    handle to plot_type_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_type_chooser contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_type_chooser
dc = handles.data_controller;
if isempty(dc.cell_numbers), return, end
mode_index = get(handles.plot_type_chooser,'Value');
strings = get(handles.plot_type_chooser,'String');
mode = char(strings(mode_index));
visualize(handles,dc.cell_numbers,dc.gen_numbers,mode);

% --- Executes during object creation, after setting all properties.
function plot_type_chooser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dt_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt_edit as text
%        str2double(get(hObject,'String')) returns contents of dt_edit as a double
value = str2double(get(hObject,'String'));
if isnumeric(value) && value > 1/12
    handles.dt = value;
    guidata(hObject,handles);
else
    value = handles.dt;
    set(hObject,'String',num2str(value));
end
dc = handles.data_controller;
dc.dt = handles.dt;
uiresume(handles.CellPopSim);


% --- Executes during object creation, after setting all properties.
function dt_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
