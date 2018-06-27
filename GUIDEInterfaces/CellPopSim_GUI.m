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

% Last Modified by GUIDE v2.5 20-Jun-2018 16:49:15

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

set(handles.plot_type_chooser,'String',{'N(t)','generation(t)','total #'});

plot(handles.graph_pane,handles.data_controller.G, ... 
                    'Layout','layered', ...
                    'EdgeLabel',handles.data_controller.G.Edges.Weight);
set(handles.graph_pane, 'xticklabel', [], 'yticklabel', []);
set(handles.graph_pane, 'xtick', [], 'ytick', []);

handles.Tmax = 10;
set(handles.Tmax_edit,'String',num2str(handles.Tmax));

handles.dt = data_controller.dt;
set(handles.dt_edit,'String',num2str(handles.dt));

set(handles.edges_table,'ColumnName',{'S','T','weight'});
set(handles.edges_table,'Data',[data_controller.S data_controller.T num2cell(data_controller.weights)]);

set(handles.phenotypes_table,'RowName',data_controller.G.Nodes.Name);
set(handles.phenotypes_table,'ColumnName',data_controller.phenotype_properties_names');
set(handles.phenotypes_table,'Data',data_controller.phenotype_properties);

set(handles.logY,'Value',1);

set(handles.show_experimental_curve,'Value',0);

set(handles.main_plot, 'xticklabel', [], 'yticklabel', []);
set(handles.main_plot, 'xtick', [], 'ytick', []);
%%%

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

Tmax = str2double(get(handles.Tmax_edit,'String'));
dc.Tmax = Tmax;
dc.simulate(24*Tmax);
handles.data_controller = dc;

% Update handles structure
guidata(hObject, handles);

if isempty(dc.cell_numbers)
    % cla(handles.main_plot,'reset');
    return;
end

mode_index = get(handles.plot_type_chooser,'Value');
strings = get(handles.plot_type_chooser,'String');
mode = char(strings(mode_index));
%
visualize(handles,mode);

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
dc = handles.data_controller;
dc.Tmax = handles.Tmax;
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


% --- Executes on button press in cancel_simulation.
function cancel_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
notify(dc,'interrupt_simulations');

% --------------------------------------------------------------------
function save_current_plot_as_fig_Callback(hObject, eventdata, handles)
% hObject    handle to save_current_plot_as_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
root_dir = pwd;
cd(dc.DefaultDirectory);

[filename,pathname] = uiputfile('*.fig','Save figure');
if pathname == 0 %if the user pressed cancelled, then we exit this callback 
    cd(root_dir);
 return 
end 

saveName = fullfile(pathname,filename);

f_new = figure;
ax_new = copyobj(handles.main_plot,f_new);
set(ax_new,'Position','default');
legend(ax_new,'-DynamicLegend');
saveas(f_new,saveName,'fig');
close(f_new);
%
cd(root_dir);

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function visualize(handles,mode)
dc = handles.data_controller;
cell_numbers = dc.cell_numbers*dc.experimental_curve_scale;
gen_numbers = dc.gen_numbers;
cla(handles.main_plot,'reset');
%
if isempty(cell_numbers), return, end;
%
is_logY = get(handles.logY,'Value');
%
show_experimental_curve = ... 
    get(handles.show_experimental_curve,'Value') && ~isempty(dc.experimental_curve);
%
switch mode
    case 'N(t)'
        for k=1:size(dc.phenotype_properties,1)
             if is_logY
                semilogy(handles.main_plot,dc.t/24,cell_numbers(k,:),'linewidth',2);
             else
                plot(handles.main_plot,dc.t/24,cell_numbers(k,:),'linewidth',2);                 
             end
             hold on;
        end
        %
        if is_logY
            semilogy(handles.main_plot,dc.t/24,sum(cell_numbers,1),'r:','linewidth',2);
        else
            plot(handles.main_plot,dc.t/24,sum(cell_numbers,1),'r:','linewidth',2);
        end
        %
        if show_experimental_curve
            te = (dc.experimental_curve(:,1) - dc.experimental_curve_shift)/24;
            Ne = dc.experimental_curve(:,2);
            hold on;
            if is_logY
                semilogy(handles.main_plot,te,Ne,'ko-','linewidth',2,'markersize',10);
            else
                plot(handles.main_plot,te,Ne,'ko-','linewidth',2,'markersize',10);
            end            
        end
        %
        hold off; 
        if show_experimental_curve
            legend(handles.main_plot, ...
                [dc.G.Nodes.Name' 'total' dc.experimental_curve_name],'fontsize',14);
        else
            legend(handles.main_plot,[dc.G.Nodes.Name' 'total'],'fontsize',14);
        end
        grid(handles.main_plot,'on');
        xlabel(handles.main_plot,'days past inception','fontsize',12);
        ylabel(handles.main_plot,'cell N(t)','fontsize',12);
        
    case 'generation(t)'
        for k=1:size(dc.phenotype_properties,1)
             gens = gen_numbers(k,:);
             range = gens>0;
             plot(handles.main_plot,dc.t(range)/24,gens(range),'linewidth',2);                          
             hold on;
        end
        hold off;
        legend(handles.main_plot,dc.G.Nodes.Name','fontsize',14);
        grid(handles.main_plot,'on');
        xlabel(handles.main_plot,'days past inception','fontsize',12);
        ylabel(handles.main_plot,'cell generation(t)','fontsize',12);
        
    case 'total #'
        Y = sum(cell_numbers,2);
        Y = Y/sum(Y)*100;
        for k=1:size(dc.phenotype_properties,1)
             plot(handles.main_plot,[k k],[0 Y(k)],'linewidth',20);                          
             hold on;
        end
        hold off;
        axis([0 numel(Y)+1 0 100]);
        grid(handles.main_plot,'on');
        set(handles.main_plot,'XGrid','off')
        legend(handles.main_plot,dc.G.Nodes.Name','fontsize',14);
        set(handles.main_plot,'xticklabel',[]);
        ylabel(handles.main_plot,'total cell number %','fontsize',12);
end
% --------------------------------------------------------------------


% --- Executes on selection change in plot_type_chooser.
function plot_type_chooser_Callback(hObject, eventdata, handles)
% hObject    handle to plot_type_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_type_chooser contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_type_chooser
mode_index = get(handles.plot_type_chooser,'Value');
strings = get(handles.plot_type_chooser,'String');
mode = char(strings(mode_index));
visualize(handles,mode);

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
if isnumeric(value) && value >= 1/12
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


% --- Executes on button press in update_model.
function update_model_Callback(hObject, eventdata, handles)
% hObject    handle to update_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
dc.phenotype_properties = get(handles.phenotypes_table,'Data');
%
D = get(handles.edges_table,'Data');
dc.weights = cell2mat(D(:,3));
dc.G = digraph(dc.S,dc.T,dc.weights);
plot(handles.graph_pane,handles.data_controller.G, ... 
                    'Layout','layered', ...
                    'EdgeLabel',handles.data_controller.G.Edges.Weight);
set(handles.graph_pane, 'xticklabel', [], 'yticklabel', []);
set(handles.graph_pane, 'xtick', [], 'ytick', []);

% --------------------------------------------------------------------
function load_model_Callback(hObject, eventdata, handles)
% hObject    handle to load_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
[fname, fpath] = uigetfile('*.xml','Load Model..',dc.DefaultDirectory);
if fpath == 0; return; end
filespec = fullfile(fpath,fname);
try
    dc.load_model(filespec);
    plot(handles.graph_pane,dc.G, ... 
                    'Layout','layered', ...
                    'EdgeLabel',dc.G.Edges.Weight);    
    set(handles.graph_pane, 'xticklabel', [], 'yticklabel', []);
    set(handles.graph_pane, 'xtick', [], 'ytick', []);                
    
    set(handles.edges_table,'Data',[dc.S dc.T num2cell(dc.weights)]);
    set(handles.edges_table,'ColumnName',{'S','T','weight'});
    
    set(handles.phenotypes_table,'RowName',dc.G.Nodes.Name);
    set(handles.phenotypes_table,'ColumnName',dc.phenotype_properties_names');
    set(handles.phenotypes_table,'Data',dc.phenotype_properties); 
    %
    handles.dt = dc.dt;
    handles.Tmax = dc.Tmax;
    set(handles.dt_edit,'String',num2str(handles.dt));
    set(handles.Tmax_edit,'String',num2str(handles.Tmax));
    %
    % show experimental curve if present
    is_logY = get(handles.logY,'Value');
       if ~isempty(dc.experimental_curve)
           if ischar(dc.experimental_curve)
               curve = eval(dc.experimental_curve);
                te = curve(:,1)/24;
                Ne = curve(:,2);          
                dc.experimental_curve = curve;
           else
            te = dc.experimental_curve(:,1)/24;
            Ne = dc.experimental_curve(:,2);
           end
            if is_logY
                semilogy(handles.main_plot,te,Ne,'ko-','linewidth',2,'markersize',10);
            else
                plot(handles.main_plot,te,Ne,'ko-','linewidth',2,'markersize',10);
            end
            grid(handles.main_plot,'on');
            xlabel(handles.main_plot,'days past inception','fontsize',12);
            ylabel(handles.main_plot,'cell N(t)','fontsize',12);
           legend(handles.main_plot,{dc.experimental_curve_name},'fontsize',14);           
            set(handles.show_experimental_curve,'Value',1);
       else
           cla(handles.main_plot,'reset');
           set(handles.main_plot, 'xticklabel', [], 'yticklabel', []);
           set(handles.main_plot, 'xtick', [], 'ytick', []);
       end
        
       set(handles.CellPopSim,'Name',['CellPopSim: ' filespec]);
catch
    errordlg('Error while trying to load model');
end


% --------------------------------------------------------------------
function save_model_Callback(hObject, eventdata, handles)
% hObject    handle to save_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
[fname, fpath] = uiputfile('*.xml','Save Model as..',dc.DefaultDirectory);
if fpath == 0; return; end
filespec = fullfile(fpath,fname);
try
    dc.save_model(filespec);
catch
    errordlg('Error while trying to save model');
end


% --- Executes on button press in logY.
function logY_Callback(hObject, eventdata, handles)
% hObject    handle to logY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logY
plot_type_chooser_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function load_experimental_curve_Callback(hObject, eventdata, handles)
% hObject    handle to load_experimental_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
[fname, fpath] = uigetfile({'*.xls;*.xlsx;*.csv'},'Load experimental curve..',dc.DefaultDirectory);
try
    filespec = fullfile(fpath,fname);
    [NUM,TXT,RAW]=xlsread(filespec);
    dc.experimental_curve = NUM;
    str = strrep(fname,'.xlsx',' ');
        str = strrep(str,'.xls',' ');
            str = strrep(str,'.csv',' ');
                str = strrep(str,'_',' ');
    dc.experimental_curve_name = str;
    disp(fname);
catch err
    disp(err.message)
end
set(handles.show_experimental_curve,'Value',1);
plot_type_chooser_Callback(hObject, eventdata, handles);






% --- Executes on button press in show_experimental_curve.
function show_experimental_curve_Callback(hObject, eventdata, handles)
% hObject    handle to show_experimental_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_experimental_curve
plot_type_chooser_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function load_state_Callback(hObject, eventdata, handles)
% hObject    handle to load_state (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
             [fname, fpath] = uigetfile({'*.mat'},'select state file',dc.DefaultDirectory);
             if fpath == 0, return, end;
             filespec = fullfile(fpath,fname);
             try
                dc.load_state(filespec);
                    mode_index = get(handles.plot_type_chooser,'Value');
                    strings = get(handles.plot_type_chooser,'String');
                    mode = char(strings(mode_index));  
                    visualize(handles,mode);                    
                        plot(handles.graph_pane,handles.data_controller.G, ... 
                                            'Layout','layered', ...
                                            'EdgeLabel',handles.data_controller.G.Edges.Weight);
                        set(handles.graph_pane, 'xticklabel', [], 'yticklabel', []);
                        set(handles.graph_pane, 'xtick', [], 'ytick', []);
             catch
                errordlg('Error while trying to load state');
             end

% --------------------------------------------------------------------
function save_state_Callback(hObject, eventdata, handles)
% hObject    handle to save_state (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
            [fname, fpath] = uiputfile('*.mat','save state as..',[dc.DefaultDirectory filesep 'CellPopSim_state']);
            if fpath == 0; return; end
            filespec = fullfile(fpath,fname);
            try
                dc.save_state(filespec);
            catch
                errordlg('Error while trying to save state');
            end







