
% God's version
classdef CellPopSim_data_controller < handle
           
    properties
        
        dt = 0.25;      % 15 minutes
        t = [];         % global time
                
        cells_data = [];        
        
        cell_types = [];
        %
        cell_cycle_phases_durations = [];
        cell_cycle_G1_duration_std = [];
        cell_cycle_G2M_duration_std = [];
        cell_cycle_Progression_threshold = [];                   
        cell_cycle_Progression_threshold_uncertainty = [];
        %
        cell_exit_probability = [];
        %
        Progression_probability = [];        
        %
        phenotypes = [];
        %
        
        cell_buffer = [];
        crnt_buffer = [];
        next_buffer = [];
        
        models = {'cancer','cancer_delayed','kidney_development','disappearance'};
        
        simulations_interrupt_flag  = false;
        
        cell_numbers = [];
        gen_numbers = [];
        
        G = []; % graph
        
    end
%--------------------------------------------------------------------                  
    events
        interrupt_simulations;
    end        
%--------------------------------------------------------------------              
    methods
%--------------------------------------------------------------------      
        function obj = CellPopSim_data_controller(varargin) 
            
            addlistener(obj,'interrupt_simulations',@obj.on_interrupt_simulations);
            obj.create_hematopoietic_graph;
            %obj.create_simple_branching_graph;
        end
%--------------------------------------------------------------------                      
        function on_interrupt_simulations(obj,~,~)
            obj.simulations_interrupt_flag = true;
        end        
%--------------------------------------------------------------------                      
        function close_request_fcn(obj,~,~)
            obj.cell_buffer = [];
            obj.crnt_buffer = [];
            obj.next_buffer = [];
            clear all;
        end
%--------------------------------------------------------------------
% 1 birth time
% 2 G2Mout time
% 3 type
% 4 generation
function Cell = get_ini_cell(obj,~,~)
         Cell = zeros(1,4);
         Cell = obj.set_type_index(Cell,1);
         Cell = obj.set_birth_and_G2Mout_times(Cell,0);
         Cell = obj.set_generation(Cell,1);
end
%--------------------------------------------------------------------
function out_rec = set_generation(obj,rec,value_to_set,~)
    out_rec = rec;
    out_rec(4) = value_to_set;
end
%--------------------------------------------------------------------
function out_rec = set_type_index(obj,rec,value_to_set,~)
    out_rec = rec;
    out_rec(3) = value_to_set;
end    
%--------------------------------------------------------------------    
function out_rec = set_birth_and_G2Mout_times(obj,rec,birth_time,~)
    out_rec = rec;
    type_index = rec(3);
    out_rec(1) = birth_time;
            TG1 = obj.cell_cycle_phases_durations(type_index,1);
            TS = obj.cell_cycle_phases_durations(type_index,2);
            TG2M = obj.cell_cycle_phases_durations(type_index,3);    
    TG1_std = obj.cell_cycle_G1_duration_std(type_index);
    TG2M_std = obj.cell_cycle_G2M_duration_std(type_index);
    out_rec(2) = birth_time + TS + normrnd(TG1,TG1_std) + normrnd(TG2M,TG2M_std);
end
%
%--------------------------------------------------------------------
function ret = get_generation(obj,rec,~)
    ret = rec(4);
end
%--------------------------------------------------------------------
function ret = get_type_index(obj,rec,~)
    ret = rec(3);
end
%--------------------------------------------------------------------    
function [t1,t2] = get_birth_and_G2Mout_times(obj,rec,~)
    t1 = rec(1);
    t2 = rec(2);    
end
%--------------------------------------------------------------------    
function [c1,c2] = update_cell(obj,in_cell,~)
    c1 = in_cell;
    c2 = in_cell;
    %
    prev_gen = obj.get_generation(in_cell);
    [~,t2] = obj.get_birth_and_G2Mout_times(in_cell);
    %
    c1 = obj.set_generation(c1,prev_gen+1);
    c2 = obj.set_generation(c2,prev_gen+1);    
    %
    % check if type is to be changed and redo settings if needed - c1
    c1 = obj.Progress(c1);
    c2 = obj.Progress(c2);
    %        
    if rand < obj.cell_exit_probability(obj.get_type_index(c1))
        c1 = [];
    else
        c1 = obj.set_birth_and_G2Mout_times(c1,t2);
    end
    %
    if rand  < obj.cell_exit_probability(obj.get_type_index(c2))
        c2 = [];
    else
        c2 = obj.set_birth_and_G2Mout_times(c2,t2);
    end
    %
end
%--------------------------------------------------------------------    
        function out_cell = Progress(obj,in_cell,~)
            %
            out_cell = in_cell;
            type_ind = obj.get_type_index(out_cell);  
            %
            Ntypes = numel(obj.cell_types);
            if type_ind == Ntypes, return, end;

            g0 = obj.cell_cycle_Progression_threshold(type_ind);
            g_scale = obj.cell_cycle_Progression_threshold_uncertainty(type_ind);
            g = obj.get_generation(in_cell);
            p = Generalised_logistic_function(double(g),double(g0),double(g_scale),1);
            %
            P = obj.Progression_probability(type_ind);
            %
            if p*P > rand
                next_type_ind = min(type_ind+1,Ntypes);                
                out_cell = obj.set_type_index(out_cell,next_type_ind);
                out_cell = obj.set_generation(out_cell,1); % local counter
            end
        end  
%--------------------------------------------------------------------    
        function simulate(obj,Nini,Tmax,~)
            %
            obj.cell_buffer = [];
            obj.crnt_buffer = [];
            obj.next_buffer = [];                        
            %
            obj.cell_buffer = zeros(2^22,4,'single');
            obj.crnt_buffer = zeros(2^21,4,'single');
            obj.next_buffer = zeros(2^21,4,'single');                        
            %
            ini_cells = [];
            for k=1:Nini
                Cell = obj.get_ini_cell();
                ini_cells = [ ini_cells; Cell];
            end
            %
            obj.crnt_buffer(1:Nini,:) = ini_cells;
            obj.cell_buffer(1:Nini,:) = ini_cells;            
            %            
            n_cells_crnt = Nini;
            Ncells = Nini;
            %            
            tic            
            hw = waitbar(0,'Creating cell pool - please wait');   
            g=0;
            while true
                g=g+1;
                %
                cnt = 0; % counter for newly added cells
                for k=1:n_cells_crnt
                    [c1,c2] = obj.update_cell(obj.crnt_buffer(k,:));
                    if ~isempty(c1)
                        cnt = cnt+1;
                        obj.next_buffer(cnt,:)=c1;
                    end
                    if ~isempty(c2)
                        cnt = cnt+1;
                        obj.next_buffer(cnt,:)=c2;
                    end                    
                end
                % save current buffer to cells
                if 1~=g
                    obj.cell_buffer(Ncells+1:Ncells+n_cells_crnt,:) = obj.crnt_buffer(1:n_cells_crnt,:);
                    Ncells = Ncells + n_cells_crnt;                    
                end
                % replace current buffer by the next buffer
                n_cells_crnt = cnt;
                obj.crnt_buffer(1:n_cells_crnt,:) = obj.next_buffer(1:n_cells_crnt,:);
                %
                Tmax_cur = min(squeeze(obj.crnt_buffer(1:n_cells_crnt,2)));
                %
                disp([g Tmax_cur Ncells]);
                %
                if Tmax_cur > Tmax, break, end
                %                                
                if ~isempty(hw), waitbar(Tmax_cur/Tmax,hw); drawnow, end
                %
                if obj.simulations_interrupt_flag
                    if ~isempty(hw), delete(hw), drawnow; end            
                    obj.cell_buffer = [];
                    obj.crnt_buffer = [];
                    obj.next_buffer = [];
                    obj.simulations_interrupt_flag = false;
                    obj.cell_numbers = [];
                    obj.gen_numbers = [];                    
                    return;
                end
            end 
            if ~isempty(hw), delete(hw), drawnow; end            
            %
            lapsed_time = toc
            %            
            N_tpoints = floor(Tmax/obj.dt);
            obj.t = (0:N_tpoints-1)*obj.dt;
            obj.cell_numbers = zeros(numel(obj.cell_types),N_tpoints);
            obj.gen_numbers = zeros(numel(obj.cell_types),N_tpoints);
            %
            hw = waitbar(0,'Harvesting time dependencies, please wait'); 
            inds = 1:N_tpoints;
            for k=1:Ncells
                % run through all cells updating every time point
                [t1,t2] = obj.get_birth_and_G2Mout_times(obj.cell_buffer(k,:));
                type_index = obj.get_type_index(obj.cell_buffer(k,:));
                gen = obj.get_generation(obj.cell_buffer(k,:));
                z = (obj.t>t1).*(obj.t<t2);
                z=z.*inds;
                z=z(0~=z);
                for m=1:numel(z)
                    obj.cell_numbers(type_index,z(m)) = obj.cell_numbers(type_index,z(m))+1;
                    obj.gen_numbers(type_index,z(m)) = obj.gen_numbers(type_index,z(m))+gen;
                end
                %
                if ~isempty(hw), waitbar(k/Ncells,hw); drawnow, end;
                %
                if obj.simulations_interrupt_flag
                    if ~isempty(hw), delete(hw), drawnow; end            
                    obj.cell_buffer = [];
                    obj.crnt_buffer = [];
                    obj.next_buffer = [];
                    obj.simulations_interrupt_flag = false;
                    obj.cell_numbers = [];
                    obj.gen_numbers = [];
                    return;
                end                
                %
            end
            if ~isempty(hw), delete(hw), drawnow; end;                        
            
            obj.gen_numbers = obj.gen_numbers./obj.cell_numbers;
            obj.gen_numbers(isnan(obj.gen_numbers))=0;
            %
        end
%--------------------------------------------------------------------    
        function create_graph(obj,~,~)
            if isempty(obj.cell_types), return, end
            %
            N = numel(obj.cell_types);            
            % for now linear only            
            S_range = 1:N-1;
            T_range = 2:N;      
            weights = ones(size(S_range))';
            obj.G = digraph(obj.cell_types(S_range),obj.cell_types(T_range),weights);            
        end
%--------------------------------------------------------------------
        function create_hematopoietic_graph(obj,~,~) %default graph            
            [~,TXT,~] = xlsread('hematopoiesis_edges.xls');
            S = TXT(:,1);
            T = TXT(:,2);
            obj.G = digraph(S,T,ones(size(T)));
        end
%--------------------------------------------------------------------            
        function create_simple_branching_graph(obj,~,~)
            names = { 'T1' 'T2'; ...
                      'T2' 'T3';
                      'T2' 'T4' };
            S = names(:,1);
            T = names(:,2);
            %
            weights = [1 0.4 0.6]';
            %
            obj.G = digraph(S,T,weights);
            %
            obj.cell_types = obj.G.Nodes.Name;
            %
            ph_T1 = phenotype('T1');
                ph_T1.cell_cycle_duration = 6;
                ph_T1.cell_cycle_duration_std = 1;
                ph_T1.cell_cycle_Progression_threshold = 12;                   
                ph_T1.cell_cycle_Progression_threshold_uncertainty = 0.5;        
                ph_T1.apoptosis_probability = 0;
                ph_T1.quescience_probability = 0; 
            %    
            ph_T2 = phenotype('T2');
                ph_T2.cell_cycle_duration = 6;
                ph_T2.cell_cycle_duration_std = 1;
                ph_T2.cell_cycle_Progression_threshold = 12;                   
                ph_T2.cell_cycle_Progression_threshold_uncertainty = 0.5;        
                ph_T2.apoptosis_probability = 0;
                ph_T2.quescience_probability = 0; 
            %                
            ph_T3 = phenotype('T3');
                ph_T3.cell_cycle_duration = 6;
                ph_T3.cell_cycle_duration_std = 1;
                ph_T3.cell_cycle_Progression_threshold = 12;                   
                ph_T3.cell_cycle_Progression_threshold_uncertainty = 0.5;        
                ph_T3.apoptosis_probability = 0;
                ph_T3.quescience_probability = 0; 
            %                
            ph_T4 = phenotype('T4');
                ph_T4.cell_cycle_duration = 6;
                ph_T4.cell_cycle_duration_std = 1;
                ph_T4.cell_cycle_Progression_threshold = 12;                   
                ph_T4.cell_cycle_Progression_threshold_uncertainty = 0.5;        
                ph_T4.apoptosis_probability = 0;
                ph_T4.quescience_probability = 0; 
            %                            
            obj.phenotypes{1} = ph_T1;
            obj.phenotypes{2} = ph_T2;
            obj.phenotypes{3} = ph_T3;
            obj.phenotypes{4} = ph_T4;
            %
            % TEST CODE TO GET DIFFERENTIATION PROBABILITIES FROM
            % CURRENT PHENOTYPE
%             for k=1:numel(obj.cell_types)
%                 name = obj.cell_types{k};
%                 current_node = findnode(obj.G,name);
%                 succ = successors(obj.G,current_node); % immediate_successors
%                 for s=1:numel(succ)
%                     edge = findedge(obj.G,current_node,succ(s));
%                     % weights between current node and successors
%                     obj.G.Edges.Weight(edge)
%                 end
%             end
            %
        end
%--------------------------------------------------------------------
    end  
end
