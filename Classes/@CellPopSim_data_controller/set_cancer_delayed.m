%--------------------------------------------------------------------                      
        function set_cancer_delayed(obj,~,~)            
            % new layout
            names = { 'T1' 'T2'; ...
                      'T2' 'T3'; ...
                      'T3' 'T4' };
            S = names(:,1);
            T = names(:,2);
            %
            weights = [1 1 1]';
            %
            obj.G = digraph(S,T,weights);            
            %
            obj.cell_types = obj.G.Nodes.Name';
            %
            ph_T1 = phenotype('T1');
                ph_T1.cell_cycle_duration = 12;
                ph_T1.cell_cycle_duration_std = 2;
                ph_T1.cell_cycle_Progression_threshold = 8;                   
                ph_T1.cell_cycle_Progression_threshold_uncertainty = 0.1;
                ph_T1.apoptosis_probability = 0;
                ph_T1.quescience_probability = 0; 
            %    
            ph_T2 = phenotype('T2');
                ph_T2.cell_cycle_duration = 13;
                ph_T2.cell_cycle_duration_std = 2;
                ph_T2.cell_cycle_Progression_threshold = 12;                   
                ph_T2.cell_cycle_Progression_threshold_uncertainty = 0.1;        
                ph_T2.apoptosis_probability = 0.5;
                ph_T2.quescience_probability = 0; 
            %                
            ph_T3 = phenotype('T3');
                ph_T3.cell_cycle_duration = 14;
                ph_T3.cell_cycle_duration_std = 2;
                ph_T3.cell_cycle_Progression_threshold = 4;                   
                ph_T3.cell_cycle_Progression_threshold_uncertainty = 0.1;        
                ph_T3.apoptosis_probability = 0;
                ph_T3.quescience_probability = 0; 
            %                
            ph_T4 = phenotype('T4');
                ph_T4.cell_cycle_duration = 15;
                ph_T4.cell_cycle_duration_std = 2;
                ph_T4.cell_cycle_Progression_threshold = 80;                   
                ph_T4.cell_cycle_Progression_threshold_uncertainty = 0.1;        
                ph_T4.apoptosis_probability = 0.2;
                ph_T4.quescience_probability = 0; 
            %                            
            obj.phenotypes{1} = ph_T1;
            obj.phenotypes{2} = ph_T2;
            obj.phenotypes{3} = ph_T3;
            obj.phenotypes{4} = ph_T4;                              
            %
            % obj.cell_types = {'T1','T2','T3','T4'}; % three progressive clones
            obj.cell_cycle_phases_durations = [1,6,5; ...   %
                                               2,6,5; ...   % 
                                               3,6,5; ...   % 
                                               4,6,5];      %
            %
            % the code below is functional for now - more or less repeats setups
            % above
            obj.cell_cycle_G1_duration_std = [2,2,2,2]; % hours
            obj.cell_cycle_G2M_duration_std = [2,2,2,2];
            %
            obj.cell_cycle_Progression_threshold_uncertainty = [0.1,0.1,0.1,0.1]; % generations   
            obj.cell_cycle_Progression_threshold = [8,12,4,80];% generations   
            %
            obj.cell_exit_probability = [0,0.5,0.2,0]; 
            %
            obj.Progression_probability = [1,1,1,1];
        end        