%--------------------------------------------------------------------                      
        function set_disappearance(obj,~,~)            
            obj.cell_types = {'T1','T2','T3'}; % three progressive clones
            obj.cell_cycle_phases_durations = [1,6,5; ...  % slow
                                               2,6,5; ...  % 
                                               3,6,5];  % 
            %
            obj.cell_cycle_G1_duration_std = [2,2,2]; % hours
            obj.cell_cycle_G2M_duration_std = [2,2,2];
            %
            obj.cell_cycle_Progression_threshold_uncertainty = [0.1,0.1,0.1]; % generations   
            obj.cell_cycle_Progression_threshold = [8,6,4];% generations   
            %
            obj.cell_exit_probability = [0,0,0.7]; 
            %
            obj.Progression_probability = [1,1,1];
            %
            obj.dt = 1;        
        end        