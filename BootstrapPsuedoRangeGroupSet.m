classdef BootstrapPsuedoRangeGroupSet < PsuedoRangeGroupSet
    %BOOTSTRAPPSUEDORANGEGROUPSET Class pertaining to bootstrapped sets of 
    %PseudoRangeGroups
    
    properties
        psuedoRangeGroupSetParent;
    end
    
    methods(Access = public)
        function obj = BootstrapPsuedoRangeGroupSet(psuedoRangeGroup, number_of_samples, number_per_group)
            
            %indices = randi([1,psuedoRangeGroup.n_measurements], number_of_samples, number_per_group);
            indices = zeros(number_of_samples, number_per_group);
            for i = 1:number_of_samples
                while true
                    row = randi([1,psuedoRangeGroup.n_measurements], 1, number_per_group);
                    if length(unique(row)) > 6
                        break;
                    end
                end
                indices(i, :) = row;
            end
            
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
