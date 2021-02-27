classdef BootstrapPsuedoRangeGroupSet < PsuedoRangeGroupSet
    %BOOTSTRAPPSUEDORANGEGROUPSET Class pertaining to bootstrapped sets of 
    %PseudoRangeGroups
    
    properties
        psuedoRangeGroupSetParent;
    end
    
    methods(Access = public)
        function obj = BootstrapPsuedoRangeGroupSet(psuedoRangeGroup, number_of_samples, number_per_group)
            
            indices = randi([1,psuedoRangeGroup.n_measurements], number_of_samples, number_per_group);
            
            [mrhos,xsats,ysats,zsats,bsats] = psuedoRangeGroup.get_data();
            
            psuedoRangeGroups = PsuedoRangeGroup.empty(size(indices,1),0);
            for i = 1:size(indices,1)
                prg = PsuedoRangeGroup(...
                    mrhos(indices(i,:)), ...
                    xsats(indices(i,:)), ...
                    ysats(indices(i,:)), ...
                    zsats(indices(i,:)), ...
                    bsats(indices(i,:)));
                psuedoRangeGroups(i) = prg;
            end
            
            obj@PsuedoRangeGroupSet(psuedoRangeGroups);
            obj.psuedoRangeGroupSetParent = psuedoRangeGroup;
        end
    end
end
