classdef RANCOPsuedoRangeGroupSet < PsuedoRangeGroupSet
    %RANCOPSUEDORANGEGROUPSET Class pertaining to RANCO sets of 
    %PseudoRangeGroups
    
    properties
        psuedoRangeGroupSetParent;
    end
    
    methods(Access = public)
        function obj = RANCOPsuedoRangeGroupSet(psuedoRangeGroup, number_per_group)
            
            % compute the set of indices that make all combinations of the
            % psuedo ranges taken number_per_group at a time
            indices = nchoosek(1:psuedoRangeGroup.n_measurements, number_per_group);
            
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
