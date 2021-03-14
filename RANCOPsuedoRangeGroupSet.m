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
        function sat_errors = satellite_errors(obj)
            % computes the individual satellite RANCO errors

            n_subsets = length(obj.psuedoRangeGroupSet);
            n_sats = obj.psuedoRangeGroupSetParent.n_measurements;

            [xrs, ns, Gs] = obj.solve_each_newton_raphson();

            sat_errors = zeros(n_subsets, n_sats);

            for i = 1:n_subsets
                xr = xrs(:,i);
                n = ns(:,i);
                G = squeeze(Gs(i,:,:));

                LOSs = obj.psuedoRangeGroupSetParent.geometry_matrix(xr);

                sat_errors(i,:) = LOSs / G * n;
            end
        end
        function threshold = threshold(obj)
            % computes the individual satellite RANCO thresholds

            n_subsets = length(obj.psuedoRangeGroupSet);
            n_sats = obj.psuedoRangeGroupSetParent.n_measurements;
            % assumed GPS estimated covariance
            W =  (diag(6*ones(obj.psuedoRangeGroupSet(1).n_measurements,1))).^2;
            sig = 0;

            [xrs, ~, Gs] = obj.solve_each_newton_raphson();

            threshold = zeros(n_subsets, n_sats);
            for i = 1:n_subsets
                xr = xrs(:,i);
                G = squeeze(Gs(i,:,:));
                LOSs = obj.psuedoRangeGroupSetParent.geometry_matrix(xr);
                for j = 1:n_sats
                    g = LOSs(j,:)';
                    threshold(i,j) = sqrt(g' / (G'/W*G) * g + sig);
                end
            end
        end
    end
end
