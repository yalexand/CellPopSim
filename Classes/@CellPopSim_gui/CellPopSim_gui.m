classdef CellPopSim_gui
           
    properties
        
        window;
        data_controller;
        
    end
    
    methods
        %     
        function obj = CellPopSim_gui()
                    
            obj.data_controller = CellPopSim_data_controller;
            
            obj.window = CellPopSim_GUI(obj.data_controller);
            
            handles = guidata(obj.window); 
            handles.window = obj.window;
            handles.use_popup = true;            
            
            guidata(obj.window,handles);
            
            close all;             
            set(obj.window,'Visible','on');
            set(obj.window,'CloseRequestFcn',@obj.close_request_fcn);                        
            % WinOnTop(obj.window,true);
            
            if isdeployed
            waitfor(obj.window);
            end
            
        end
        %        
        function close_request_fcn(obj,~,~)            
            handles = guidata(obj.window);
            obj.data_controller = [];
            % actually close window
            delete(handles.window);            
        end
        %               
    end
    
end
