classdef PsuedoRangeGroup
    %PSUEDORANGEGROUP Methods pertaining to PsuedoRangeGroups
    
    properties(Access = protected)
        mrhos
        xsats
        ysats
        zsats
        bsats
        n_measurements
        data_matrix
        tolerance = 1e-3;
        max_iter = 50;
    end
    
    methods(Access = public)
        function obj = PsuedoRangeGroup(rhos, xsats, ysats, zsats, bsats)
            % Constructs a PsuedoRangeGroup Class from input vectors
            % containing the measured psuedoranges, and satellite states
            
            obj.n_measurements = size(rhos, 1);
            
            assert(size(xsats, 1) == obj.n_measurements, "Inputs must each have same length.");
            assert(size(ysats, 1) == obj.n_measurements, "Inputs must each have same length.");
            assert(size(zsats, 1) == obj.n_measurements, "Inputs must each have same length.");
            assert(size(bsats, 1) == obj.n_measurements, "Inputs must each have same length.");
            
            obj.mrhos = rhos;
            obj.xsats = xsats;
            obj.ysats = ysats;
            obj.zsats = zsats;
            obj.bsats = bsats;
        end
        function xr = solve_newton_raphson(obj)
            % Computes the position solution with the Newton Raphson
            % method.
            
            xr = zeros(4,1);
            
            i = 0;
            while true
                
                drho = obj.mrhos - obj.predict_psuedoranges(xr);
                
                G = obj.geometry_matrix(xr);
                
                del_xb = G \ drho;
                
                xr = xr + del_xb;
                
                if norm(del_xb) < obj.tolerance
                    break;
                end
                
                if i == obj.max_iter
                    error("Reached Max Iteration");
                end
                
                i = i+1;
            end
        end
    end
    
    methods(Access = protected)
        function G = geometry_matrix(obj, x0)
            % Constructs the geometry matrix give an initial guess x0 and
            % the available satellite state data.
            
            G = ones(obj.n_measurements, 4);
            
            dx = obj.xsats - x0(1);
            dy = obj.ysats - x0(2);
            dz = obj.zsats - x0(3);
            
            G(:,1) = -dx;
            G(:,2) = -dy;
            G(:,3) = -dz;
            
            rho0s = sqrt(sum(G(:,1:3).^2, 2));
            
            G(:,1:3) = G(:,1:3) ./ rho0s;
        end
        
        function rhos = predict_psuedoranges(obj, x0)
            % Computes the predicted psuedoranges given an initial guess x0
            % and the available satellite state data.
            
            dx = obj.xsats - x0(1);
            dy = obj.ysats - x0(2);
            dz = obj.zsats - x0(3);
            b = x0(4);
            
            rhos = sqrt(sum(dx.^2 + dy.^2 + dz.^2, 2)) + b - obj.bsats;
        end  
    end
end
