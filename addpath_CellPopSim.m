function addpath_CellPopSim()
        
    if ~isdeployed

            thisdir = fileparts( mfilename( 'fullpath' ) );
            addpath( thisdir,...
                    [thisdir filesep 'Classes'],...
                    [thisdir filesep 'GUIDEInterfaces'],...                    
                    [thisdir filesep 'HelperFunctions']); 
    end
            
end