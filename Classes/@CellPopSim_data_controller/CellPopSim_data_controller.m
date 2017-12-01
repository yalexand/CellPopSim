
% God's version
classdef CellPopSim_data_controller < handle
           
    properties
        
        dt = 0.25;      % 15 minutes
        
        Tmax = 6;       % days
        
        t = [];         % global time
                
        cells_data = [];        
                
        cell_buffer = [];
        crnt_buffer = [];
        next_buffer = [];
                                        
        G = []; % graph        
        S = [];
        T = [];
        weights = [];

        phenotype_properties_names = { ...
        'Tc', ...
        'Tc_std', ...
        'Gt', ...
        'Gt_std', ...
        'A', ...
        'Q'}
        
        phenotype_properties = [];
                        
        simulations_interrupt_flag  = false;
        
        cell_numbers = [];
        gen_numbers = [];
                
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
        end
%--------------------------------------------------------------------                      
        function save_model(obj,fullfilename,~)
            settings = [];
            %
            settings.dt = obj.dt;
            settings.Tmax = obj.Tmax;            
            settings.S = obj.S;
            settings.T = obj.T;
            settings.weights = obj.weights;
            settings.phenotype_properties_names = obj.phenotype_properties_names;
            settings.phenotype_properties = obj.phenotype_properties;            
            %
            xml_write(fullfilename, settings);            
        end        
%--------------------------------------------------------------------                      
        function load_model(obj,fullfilename,~)
            if ~exist(fullfilename,'file'), return, end;
            [ settings, ~ ] = xml_read(fullfilename);    
            %
            obj.dt = settings.dt;
            obj.Tmax = settings.Tmax;            
            obj.S = settings.S';
            obj.T = settings.T';
            obj.weights = settings.weights;
            obj.phenotype_properties_names = settings.phenotype_properties_names;
            obj.phenotype_properties = settings.phenotype_properties;
            %
            obj.G = digraph(obj.S,obj.T,obj.weights);
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
     Tc = obj.phenotype_properties(type_index,1);
     Tc_std = obj.phenotype_properties(type_index,2);
     out_rec(2) = birth_time + normrnd(Tc,Tc_std);
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
    pA = obj.phenotype_properties(obj.get_type_index(c1),5); % exit probability
    if rand < pA
        c1 = [];
    else
        c1 = obj.set_birth_and_G2Mout_times(c1,t2);
    end
    %
    pA = obj.phenotype_properties(obj.get_type_index(c2),5); % exit probability    
    if rand < pA
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
              Ntypes = size(obj.phenotype_properties,1);
              curnt = obj.get_type_index(out_cell);
              if curnt == Ntypes, return, end;
              %
              g0 = obj.phenotype_properties(curnt,3);
              g_scale = obj.phenotype_properties(curnt,4);
              g = obj.get_generation(in_cell);
              p = Generalised_logistic_function(double(g),double(g0),double(g_scale),1);
              if p > rand              
                  %
                  succ = successors(obj.G,curnt); % immediate_successors
                  N_succ = numel(succ);
                  if 0==N_succ
                      return;
                  elseif 1==N_succ
                      nxt = succ;
                  else
                      w = zeros(1,N_succ); % define weights
                      for s=1:N_succ
                        edge = findedge(obj.G,curnt,succ(s));
                        % weights between current node and successors
                        w(s) = obj.G.Edges.Weight(edge);
                      end
                      nxt = curnt + roulette_wheel(w,1);
                  end
                  %
                  out_cell = obj.set_type_index(out_cell,nxt);
                  out_cell = obj.set_generation(out_cell,1);
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
            obj.cell_numbers = zeros(size(obj.phenotype_properties,1),N_tpoints);
            obj.gen_numbers = zeros(size(obj.phenotype_properties,1),N_tpoints);
            
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
        function create_hematopoietic_graph(obj,~,~) %default graph            
            [~,TXT,~] = xlsread('hematopoiesis_edges.xls');
            obj.S = TXT(:,1);
            obj.T = TXT(:,2);
            obj.weights = ones(size(obj.T));
            obj.G = digraph(obj.S,obj.T,obj.weights);
            %
            obj.phenotype_properties = zeros(numel(obj.G.Nodes.Name),numel(obj.phenotype_properties_names));
        end
%--------------------------------------------------------------------            
    end  
end
