classdef PsuedoRangeGroup < handle
    %PSUEDORANGEGROUP Methods pertaining to PsuedoRangeGroups
    
    properties(Access = protected)
        mrhos
        xsats
        ysats
        zsats
        bsats
        
        tolerance = 1e-6;
        max_iter = 50;
    end
    
    properties(GetAccess = public, SetAccess = protected)
        n_measurements
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
        function [xr, n, G] = solve_newton_raphson(obj)
            % Computes the position solution with the Newton Raphson
            % method.
            
            if length(obj.mrhos) < 4
                xr = nan(4,1);
                n = nan(4,1);
                G = nan;
                return;
            end
            
            xr = zeros(4,1);
            
            i = 0;
            while true
                
                drho = obj.mrhos - obj.predict_psuedoranges(xr);
                
                G = obj.geometry_matrix(xr);
                
                del_xb = G \ drho;
                
                if norm(del_xb) < obj.tolerance
                    break;
                end
                
                xr = xr + del_xb;

                if i == obj.max_iter
                    error("Reached Max Iteration");
                end
                
                i = i+1;
            end
            n = obj.mrhos - obj.predict_psuedoranges(xr);
        end
        function [mrhos,xsats,ysats,zsats,bsats] = get_data(obj)
            % get data associated with this object
            mrhos = obj.mrhos;
            xsats = obj.xsats;
            ysats = obj.ysats;
            zsats = obj.zsats;
            bsats = obj.bsats;
        end
    end
    
    methods(Access = public)
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
        
        function [S,PDOP,TDOP,GDOP] = DOPcalcs(obj, x0)
            % calculates dilution of precision values and sigma matrix given geometry matrix
            % S = sigma matrix
            % PDOP = position dilution of precision
            % TDOP = time dilution of precision
            % GDOP = geometric dilution of precision
            
            % geometry matrix
            G = obj.geometry_matrix(x0);
            
            % forming the matrix
            DOPmatrix = inv(G'*G);
            
            % elements along the diagonal
            XDOP2 = DOPmatrix(1,1);
            YDOP2 = DOPmatrix(2,2);
            ZDOP2 = DOPmatrix(3,3);
            TDOP2 = DOPmatrix(4,4);
            
            % DOP values
            PDOP = sqrt(XDOP2+YDOP2+ZDOP2);
            TDOP = sqrt(TDOP2);
            GDOP = sqrt((PDOP^2)+(TDOP^2));
            
            % user range error sigma
            sigma_URE = 6.7/3;
            
            % sigma matrix with sigmas along diagonal
            S = (sigma_URE^2).*DOPmatrix;
            
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
