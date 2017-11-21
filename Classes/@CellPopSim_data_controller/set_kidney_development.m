% corrected
       function set_kidney_development(obj,~,~)            
            %                                                    
            % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4903905/
            % "..early cell cycles are exceptionally fast: 
            % nuclei undergo 13 mitotic divisions in 2–2.5 hours, whereas an average tissue culture cell 
            % takes 8–24 hours to go through one cell cycle"
            % (that gives Tc about 0.2h = 2.5h/13) 
            %
            obj.cell_types = {'jump-start','carrier','CM','RV'}; % fast/slow "Cap Mesenchyme", and "Renal Ventricular"
            obj.cell_cycle_phases_durations = [0.01,0.2,0.01; ...
                                               5,6,4; ...      % 
                                               10,6,4; ...
                                               900,900,900];
            %
            obj.cell_cycle_G1_duration_std = [0.001,1,1,1]; % hours
            obj.cell_cycle_G2M_duration_std = [0.001,1,1,1];
            %
            obj.cell_cycle_Progression_threshold_uncertainty = [1,.2,1,1]; % generations   
            %
            obj.cell_cycle_Progression_threshold = [13,32-13,39.3-32,100];% 18Oct            
            %
            obj.cell_exit_probability = [0,0.5,0,0]; % probability to exit after mitosis
            %
            obj.Progression_probability = [1,1,1,1];
            %
            %obj.dt = 2;                                    
        end
