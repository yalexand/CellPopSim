classdef phenotype % < handle           
    properties               
        name = [];
        cell_cycle_duration = [];
        cell_cycle_duration_std = [];
        cell_cycle_Progression_threshold = [];                   
        cell_cycle_Progression_threshold_uncertainty = [];
        %
        apoptosis_probability = [];        
        quescience_probability = []; 
    end
    
    methods        
    %--------------------------------------------------------------------      
        function obj = phenotype(varargin) 
            obj.name = varargin{1};
        end
    %--------------------------------------------------------------------  
    end
    
end
