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
        function [xrs, ns, Gs] = solve_each_newton_raphson(obj)
            N = length(obj.psuedoRangeGroupSet);
            xrs = zeros(4,N);
            ns = zeros(obj.psuedoRangeGroupSet(1).n_measurements,N);
            [~,~,G0] = obj.psuedoRangeGroupSet(1).solve_newton_raphson();
            Gs = zeros(N, size(G0,1), size(G0,2));
            for i = 1:N
                [xr, n, G] = obj.psuedoRangeGroupSet(i).solve_newton_raphson();
                xrs(:,i) = xr;
                ns(:,i) = n;
                Gs(i,:,:) = G;
            end
        end
        function xr_mean = mean_newton_raphson_solution(obj)
            xrs = obj.solve_each_newton_raphson();
            xr_mean = mean(xrs,2);
        end
        function xr_cov = covariance_newton_raphson_solution(obj)
            xrs = obj.solve_each_newton_raphson();
            xr_cov = cov(xrs');
        end
        function Gs = geometry_matrices(obj)
            number_per_group = size(obj.psuedoRangeGroupSet(1).geometry_matrix([0;0;0;0]),1);
            n = length(obj.psuedoRangeGroupSet);
            Gs = zeros(n, number_per_group, 4);
            for i = 1:n
                x0 = obj.psuedoRangeGroupSet(i).solve_newton_raphson();
                Gs(i, :, :) = obj.psuedoRangeGroupSet(i).geometry_matrix(x0);
            end
        end
    end
end

