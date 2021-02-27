classdef PsuedoRangeGroupSet < handle
    %BOOTSTRAPPSUEDORANGEGROUPSET Class pertaining to sets of
    %PseudoRangeGroups
    
    properties
        psuedoRangeGroupSet;
    end
    
    methods(Access = public)
        function obj = PsuedoRangeGroupSet(psuedoRangeGroupArray)
            obj.psuedoRangeGroupSet = psuedoRangeGroupArray;
        end
        function xr_mean = mean_newton_raphson_solution(obj)
            n = length(obj.psuedoRangeGroupSet);
            xrs = zeros(4,n);
            for i = 1:n
                xrs(:,i) = obj.psuedoRangeGroupSet(i).solve_newton_raphson();
            end
            xr_mean = mean(xrs,2);
        end
        function xr_cov = covariance_newton_raphson_solution(obj)
            n = length(obj.psuedoRangeGroupSet);
            xrs = zeros(n,4);
            for i = 1:n
                xrs(i,:) = obj.psuedoRangeGroupSet(i).solve_newton_raphson()';
            end
            xr_cov = cov(xrs);
        end
    end
end

